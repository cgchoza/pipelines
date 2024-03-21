#Copyright (C) 2022 Inter-University Institute for Data Intensive Astronomy
#See processMeerKAT.py for license details.

import sys
import os
import shutil

import config_parser
import bookkeeping
from crosscal_scripts.config import CONFIG_PATH
from casarecipes.almapolhelpers import xyamb
import numpy as np

from casatasks import *
logfile=casalog.logfile()
from casatools import msmetadata
msmd = msmetadata()

import logging
from time import gmtime
logging.Formatter.converter = gmtime
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s: %(message)s", level=logging.INFO)

def qu_polfield(polfield, visname):
    """
    Given the pol source name and the reference frequency, returns the fractional Q and U
    calculated from Perley & Butler 2013
    """

    msmd.open(visname)
    meanfreq = msmd.meanfreq(0, unit='GHz')
    msmd.done()

    if polfield in ["3c286", "1328+307", "1331+305", "J1331+3030"]:
        #f_coeff=[1.2515,-0.4605,-0.1715,0.0336]    # coefficients for model Stokes I spectrum from Perley and Butler 2013
        perley_frac = np.array([0.086, 0.098, 0.101, 0.106, 0.112, 0.115, 0.119, 0.121, 0.123])
        perley_f = np.array([1.02, 1.47, 1.87, 2.57, 3.57, 4.89, 6.68, 8.43, 11.3])
        pa_polcal = np.array([33.0]*8 + [34.0])
    elif polfield in ["3C138", "0518+165", "0521+166", "J0521+1638"]:
        #f_coeff=[1.0332,-0.5608,-0.1197,0.041]    # coefficients for model Stokes I spectrum from Perley and Butler 2013
        perley_frac = np.array([0.056,0.075,0.084,0.09,0.104,0.107,0.10])
        perley_f = np.array([1050,1450,1640,1950,2450,2950,3250])
        pa_polcal = np.array([-14.0,-11.0,-10.0,-10.0,-10.0,-10.0,-10.0])
    elif polfield in ["3c48", "0134+329", "0137+331", "J0137+3309"]:
        perley_frac = np.array([0.003, 0.005, 0.007])
        perley_f = np.array([1050,1450,1640])
        pa_polcal = np.array([25, 140, -5])
    elif polfield == "J1130-1449": #Manual model from Russ Taylor, taken from MeerKAT polarisation calibrator project
        perley_frac = np.array([0.038, 0.050, 0.056])
        perley_f = np.array([1050, 1450, 1640])
        pa_polcal = np.array([145, 66, 45])
    else:
        # This should never happen.
        raise ValueError("Invalid polarization field. Exiting.")

    p = np.polyfit(perley_f, perley_frac, deg=2)
    p = np.poly1d(p)

    pa = np.polyfit(perley_f, pa_polcal, deg=2)
    pa = np.poly1d(pa)

    #polpoly = np.poly1d(polcoeffs)
    #polval = polpoly(meanfreq)

    # BEWARE: Stokes I coeffs are in log-log space, so care must be taken while converting to linear.
    # They are in np.log10 space, not np.log space!
    # Coefficients are from Perley-Butler 2013 (An Accurate Flux Density Scale from 1 to 50 GHz)
    #stokesIpoly = np.poly1d(stokesIcoeffs)
    #stokesIval = stokesIpoly(np.log10(meanfreq))
    ## in Jy
    #stokesIval = np.power(10, stokesIval)

    q = p(meanfreq) * np.cos(2*np.deg2rad(pa(meanfreq)))
    u = p(meanfreq) * np.sin(2*np.deg2rad(pa(meanfreq)))

    return q, u

def do_cross_cal(visname, fields, calfiles, referenceant, caldir,
        minbaselines, standard):

    polfield = bookkeeping.polfield_name(visname)
    print(fields)
    print(polfield)
    if polfield == '':
        polfield = fields.secondaryfield
    else:
        polqu = qu_polfield(polfield, visname)

    if not os.path.isdir(caldir):
        os.makedirs(caldir)
    elif not os.path.isdir(caldir+'_round1'):
        os.rename(caldir,caldir+'_round1')
        os.makedirs(caldir)

    base = visname.replace('.ms', '')
    dtempfile   = os.path.join(caldir, base+'.dtempcal')
    xy0ambpfile = os.path.join(caldir, base+'.xyambcal')
    xy0pfile    = os.path.join(caldir, base+'.xycal')
    xpfile      = os.path.join(caldir, base+'.xfcal')
    kcrossfile = os.path.join(caldir, base+'.Kcross0')


    logger.info(" starting bandpass -> %s" % calfiles.bpassfile)
    bandpass(vis=visname, caltable = calfiles.bpassfile,
            field = fields.bpassfield, refant = referenceant,
            minblperant = minbaselines, solnorm = False,  solint = '10min',
            combine = 'scan', bandtype = 'B', fillgaps = 8,
            parang = False, append = False)
    bookkeeping.check_file(calfiles.bpassfile)
    flagdata(vis=calfiles.bpassfile, datacolumn='CPARAM', mode='rflag', timedevscale=5.0, freqdevscale=5.0, action='apply')

    logger.info("starting \'Dflls\' polcal -> %s"  % calfiles.dpolfile)
    polcal(vis=visname, caltable = calfiles.dpolfile, field = fields.bpassfield,
            refant = '', solint = 'inf', combine = 'scan',
            poltype = 'Dflls', preavg= 200.0,
            gaintable = [calfiles.bpassfile],
            gainfield = [fields.bpassfield],
            append = False)
    bookkeeping.check_file(calfiles.dpolfile)
    flagdata(vis=calfiles.dpolfile, datacolumn='CPARAM', mode='rflag', timedevscale=5.0, freqdevscale=5.0, action='apply')

    logger.info(" starting gain calibration\n -> %s" % calfiles.gainfile)
    gaincal(vis=visname, caltable = calfiles.gainfile,
            field = fields.gainfields, refant = referenceant,
            minblperant = minbaselines, solnorm = False,  gaintype = 'T',
            solint = 'inf', combine = '', calmode='ap',
            gaintable=[calfiles.bpassfile,calfiles.dpolfile],
            gainfield=[fields.bpassfield,fields.bpassfield],
            parang = False, append = False)
    bookkeeping.check_file(calfiles.gainfile)

    print("The polfield is: ", polfield)
    if polfield != fields.secondaryfield:
        logger.info(" starting pol calibrator gain calibration\n -> %s" % calfiles.gainfile)
        gaincal(vis=visname, caltable = calfiles.gainfile,
                field = polfield, refant = referenceant,
                minblperant = minbaselines, solnorm = False,  gaintype = 'T',
                solint = 'inf', combine = '', calmode='ap',
                gaintable=[calfiles.bpassfile,calfiles.dpolfile],
                gainfield=[fields.bpassfield,fields.bpassfield],
                parang = False, append = True)
        bookkeeping.check_file(calfiles.gainfile)

    flagdata(vis=calfiles.gainfile, datacolumn='CPARAM', mode='rflag', timedevscale=5.0, freqdevscale=5.0, action='apply')
    # Only run fluxscale if bootstrapping
    if len(fields.gainfields.split(',')) > 1:
        logger.info(" starting fluxscale -> %s", calfiles.fluxfile)
        rmtables(os.path.join(caldir, calfiles.fluxfile))
        fluxscale(vis=visname, caltable=calfiles.gainfile,
                reference=[fields.fluxfield], transfer='',
                fluxtable=calfiles.fluxfile, append=False, display=False,
                listfile = os.path.join(caldir,'fluxscale_xy_yx.txt'))
        bookkeeping.check_file(calfiles.fluxfile)


    # # Best scan to calibrate cross-hands will be where the polarization signal is 
    # # # minimum in XX and YY (i.e., maximum in XY and YX); find the scan using the
    # # # gain calibration for the phase/polarization calibrator
    # # # This code taken from ALMA pipeline
    # tb.open(calfiles.gainfile)
    # scans = tb.getcol('SCAN_NUMBER')
    # gains = np.squeeze(tb.getcol('CPARAM'))
    # tb.close()
    # scan_list = np.array(list(set(scans)))
    # ratios = np.zeros(len(scan_list))
    # for si, s in enumerate(scan_list):
    #         filt = scans == s
    #         ratio = np.sqrt(np.average(np.power(np.abs(gains[0,filt])/np.abs(gains[1,filt])-1.0,2.)))
    #         ratios[si] = ratio

    # best_scan_index = np.argmin(ratios)
    # best_scan = scan_list[best_scan_index]
    # print(f"Scan with highest expected X-Y signal: {best_scan}")

    # # Kcross calibration
    # gaincal(vis=visname, caltable=kcrossfile, refant=referenceant, solint='inf', gaintype='KCROSS', scan=str(best_scan), smodel=[1, 0, 1, 0], calmode='ap', 
    #         minblperant=1, refantmode='strict', parang=True)
        

    if polfield == fields.secondaryfield:
        # Cannot resolve XY ambiguity so write into final file directly
        xyfile = xy0pfile
    else:
        xyfile = xy0ambpfile

    logger.info("\n Starting x-y phase calibration\n -> %s" % xy0ambpfile)
    gaincal(vis=visname, caltable = xyfile, field = polfield,
            refant = referenceant, solint = 'inf', combine = 'scan',
            gaintype = 'XYf+QU', minblperant = minbaselines,
            preavg = 120.0,
            gaintable = [calfiles.bpassfile, calfiles.dpolfile, calfiles.gainfile],
            gainfield = [fields.bpassfield, fields.bpassfield, polfield],
            append = False)
    bookkeeping.check_file(xyfile)

    if polfield != fields.secondaryfield:
        logger.info("\n Check for x-y phase ambiguity.")
        logger.info("Polarization qu is {0}".format(polqu))
        S = xyamb(xytab=xy0ambpfile, qu=polqu, xyout = xy0pfile)
        #logger.info("smodel = {0}".format(S))
        flagdata(vis=xy0pfile, datacolumn='CPARAM', mode='rflag', timedevscale=5.0, freqdevscale=5.0, action='apply')
    else:
        flagdata(vis=xyfile, datacolumn='CPARAM', mode='rflag', timedevscale=5.0, freqdevscale=5.0, action='apply')

################### RUN IT DOWN HERE ###################

taskvals,config = config_parser.parse_config(filename=CONFIG_PATH)
visname = config['data']['vis'].strip("'")

calfiles, caldir = bookkeeping.bookkeeping(visname)
fields = bookkeeping.get_field_ids(config['fields'])

minbaselines = taskvals['crosscal']['minbaselines']
standard = taskvals['crosscal']['standard']
refant = taskvals['crosscal']['refant']

do_cross_cal(visname, fields, calfiles, f"'{refant}'", caldir, minbaselines, standard)


