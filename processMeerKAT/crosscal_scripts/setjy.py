#Copyright (C) 2022 Inter-University Institute for Data Intensive Astronomy
#See processMeerKAT.py for license details.

import os, sys, shutil

import bookkeeping
import config_parser
from crosscal_scripts.config import CONFIG_PATH
import numpy as np
import logging
from time import gmtime
logging.Formatter.converter = gmtime

from casatasks import *
logfile=casalog.logfile()
casalog.setlogfile('logs/{SLURM_JOB_NAME}-{SLURM_JOB_ID}.casa'.format(**os.environ))
from casatools import msmetadata
import casampi
msmd = msmetadata()

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s: %(message)s", level=logging.INFO)

def linfit(xInput, xDataList, yDataList):
    """

    """
    y_predict = np.poly1d(np.polyfit(xDataList, yDataList, 1))
    yPredict = y_predict(xInput)
    return yPredict


def do_setjy(visname, spw, fields, standard, dopol=False):

    delmod(vis=visname) #if this isn't called, setjy job completes but has exit code 1; clearcal(vis=visname) also works

    fluxlist = ["J0408-6545", "0408-6545", ""]

    msmd.open(visname)
    fnames = fields.fluxfield.split(",")
    for fname in fnames:
        if fname.isdigit():
            fname = msmd.namesforfields(int(fname))

    do_manual = False
    for ff in fluxlist:
        if ff in fnames:
            setjyname = ff
            do_manual = True
            break
        else:
            setjyname = fields.fluxfield.split(",")[0]

    if do_manual:
        smodel = [17.066, 0.0, 0.0, 0.0]
        spwMeanFreq = msmd.meanfreq(0, unit='GHz')
        spix = [-1.179]    # check to see if this is a problem
        reffreq = f"{spwMeanFreq}GHz"

        logger.info("Using manual flux density scale - ")
        logger.info("Flux model: %s ", smodel)
        logger.info("Spix: %s", spix)
        logger.info("Ref freq %s", reffreq)

        setjy(vis=visname,field=setjyname,scalebychan=True,standard="manual",fluxdensity=smodel,spix=spix,reffreq=reffreq)
    else:
        setjy(vis=visname, field=setjyname, spw=spw, scalebychan=True, standard=standard)

    fieldnames = msmd.fieldnames()

    if dopol:
        # Check if 3C286 exists in the data
        is3C286 = False
        try:
            calibrator_3C286 = list(set(["3C286", "1328+307", "1331+305", "J1331+3030"]).intersection(set(fieldnames)))[0]
        except IndexError:
            calibrator_3C286 = []

        if len(calibrator_3C286):
            is3C286 = True
            id3C286 = str(msmd.fieldsforname(calibrator_3C286)[0])

        if is3C286:
            logger.info("Detected calibrator name(s):  %s" % calibrator_3C286)
            logger.info("Flux and spectral index taken/calculated from:  https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/fdscale")
            logger.info("Estimating polarization index and position angle of polarized emission from linear fit based on: Perley & Butler 2013 (https://ui.adsabs.harvard.edu/abs/2013ApJS..204...19P/abstract)")
            # central freq of spw
            spwMeanFreq = msmd.meanfreq(0, unit='GHz')
            freqList = np.array([1.02, 1.47, 1.87, 2.57, 3.57, 4.89, 6.68, 8.43, 11.3])

            # fractional linear polarisation
            fracPolList = np.array([0.086, 0.098, 0.101, 0.106, 0.112, 0.115, 0.119, 0.121, 0.123])
            polindex = linfit(spwMeanFreq, freqList, fracPolList)
            logger.info("Predicted polindex at frequency %s: %s", spwMeanFreq, polindex)
            # position angle of polarized intensity
            polPositionAngleList = np.array([33.0]*8 + [34.0])
            polangle = linfit(spwMeanFreq, freqList, polPositionAngleList)
            logger.info("Predicted pol angle at frequency %s: %s", spwMeanFreq, polangle)

            reffreq = f"{spwMeanFreq}GHz"
            logger.info("Ref freq %s", reffreq)
            setjy(vis=visname,
                field=id3C286,
                scalebychan=True,
                standard="manual",
                fluxdensity=[-14.6, 0.0, 0.0, 0.0],
                #spix=-0.52, # between 1465MHz and 1565MHz
                reffreq=reffreq,
                polindex=[polindex],
                polangle=[polangle],
                rotmeas=0)


        # Check if 3C138 exists in the data
        is3C138 = False
        try:
            calibrator_3C138 = list(set(["3C138", "0518+165", "0521+166", "J0521+1638"]).intersection(set(fieldnames)))[0]
        except IndexError:
            calibrator_3C138 = []

        if len(calibrator_3C138):
            is3C138 = True
            id3C138 = str(msmd.fieldsforname(calibrator_3C138)[0])

        if is3C138:
            logger.info("Detected calibrator name(s):  %s" % calibrator_3C138)
            logger.info("Flux and spectral index taken/calculated from:  https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/fdscale")
            logger.info("Estimating polarization index and position angle of polarized emission from linear fit based on: Perley & Butler 2013 (https://ui.adsabs.harvard.edu/abs/2013ApJS..204...19P/abstract)")
            # central freq of spw
            spwMeanFreq = msmd.meanfreq(0, unit='GHz')
            freqList = np.array([1.05, 1.45, 1.64, 1.95])
            # fractional linear polarisation
            fracPolList = [0.056, 0.075, 0.084, 0.09]
            polindex = linfit(spwMeanFreq, freqList, fracPolList)
            logger.info("Predicted polindex at frequency %s: %s", spwMeanFreq, polindex)
            # position angle of polarized intensity
            polPositionAngleList = [-14, -11, -10, -10]
            polangle = linfit(spwMeanFreq, freqList, polPositionAngleList)
            logger.info("Predicted pol angle at frequency %s: %s", spwMeanFreq, polangle)

            reffreq = f"{spwMeanFreq}GHz"
            logger.info("Ref freq %s", reffreq)
            setjy(vis=visname,
                field=id3C138,
                scalebychan=True,
                standard="manual",
                fluxdensity=[-8.26, 0.0, 0.0, 0.0],
                #spix=-0.57,  # between 1465MHz and 1565MHz
                reffreq=reffreq,
                polindex=[polindex],
                polangle=[polangle],
                rotmeas=0)

    msmd.done()


def main(args,taskvals):

    if os.path.exists(os.path.join(os.getcwd(), "caltables")):
        shutil.rmtree(os.path.join(os.getcwd(), "caltables"))

    taskvals,config = config_parser.parse_config(filename=CONFIG_PATH)
    visname = config['data']['vis'].strip("'")

    calfiles, caldir = bookkeeping.bookkeeping(visname)
    fields = bookkeeping.get_field_ids(config['fields'])

    spw = config['crosscal']['spw'].split(" ")[0]
    standard = config['crosscal']['standard'].split(" ")[0]
    dopol = config["run"]["dopol"].split(" ")[0]

    do_setjy(visname, spw, fields, standard, dopol)

if __name__ == '__main__':

    bookkeeping.run_script(main,logfile)
