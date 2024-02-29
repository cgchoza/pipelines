#Copyright (C) 2022 Inter-University Institute for Data Intensive Astronomy
#See processMeerKAT.py for license details.

#!/usr/bin/env python3
import sys
import os
import numpy as np

import processMeerKAT
import config_parser

from casatasks import *
from casatools import msmetadata,table,measures,quanta

logger = processMeerKAT.logger

def get_fields(msmd=None, config=None):

    """Extract field numbers from config file, including calibrators for bandpass, flux, phase & amplitude, and the target. Only the
    target allows for multiple field IDs, while all others extract the field with the most scans and put all other IDs as target fields.

    Arguments:
    ----------
    MS : str
        Input MeasurementSet (relative or absolute path).

    Returns:
    --------
    fieldIDs : dict
        fluxfield : int
            Field for total flux calibration.
        bpassfield : int
            Field for bandpass calibration.
        phasecalfield : int
            Field for phase calibration.
        targetfields : int
            Target field."""

    fieldIDs = {}

    #Set default for any missing intent as field for intent CALIBRATE_FLUX

    fluxcal = config['fields']['fluxfield']
    if fluxcal == '':
        logger.error('A flux calibrator must be set in the config file.')
        return fieldIDs
    else:
        default = fluxcal

    fieldIDs['fluxfield'] = fluxcal.strip("'")
    fieldIDs['bpassfield'] = config["fields"]["bpassfield"].strip("'")
    fieldIDs['phasecalfield'] = config["fields"]["phasecalfield"].strip("'")
    fieldIDs['targetfields'] = config["fields"]["targetfields"].strip("'")
    fieldIDs['polfield'] = config["fields"]["polcalfield"].strip("'")
        
    # fieldIDs['fluxfield'] = get_field(MS,'CALIBRATE_FLUX','fluxfield',extra_fields)
    # fieldIDs['bpassfield'] = get_field(MS,'CALIBRATE_BANDPASS','bpassfield',extra_fields,default=default)
    # fieldIDs['phasecalfield'] = get_field(MS,phasecal_intent,'phasecalfield',extra_fields,default=default)
    # fieldIDs['targetfields'] = get_field(MS,'TARGET','targetfields',extra_fields,default=default,multiple=True)

    return fieldIDs

def check_refant(MS, refant,config,msmd=None, warn=True):

    """Check if reference antenna exists, otherwise throw an error or display a warning.

    Arguments:
    ----------
    MS : str
        Input MeasurementSet (relative or absolute path).
    refant: str
        Input reference antenna.
    config : str
        Path to config file.
    warn : bool, optional
        Warn the user? If False, raise ValueError."""

    try:
        refant = int(refant)
    except ValueError: # It's not an int, but a str
        pass

    msmd.open(MS)
    if type(refant) is str:
        ants = msmd.antennanames()
    else:
        ants = msmd.antennaids()

    if refant not in ants:
        err = "Reference antenna '{0}' isn't present in input dataset '{1}'. Antennas present are: {2}. Try 'm052' or 'm005' if present, or ensure 'calcrefant=True' and 'calc_refant.py' script present in '{3}'.".format(refant,MS,ants,config)
        if warn:
            logger.warning(err)
        else:
            raise ValueError(err)
    else:
        logger.info("Using reference antenna '{0}'.".format(refant))

def parang_coverage(vis, calfield, msmd=None):

    """Check whether the parallactic angle coverage of the phase calibrator field is > 30 degrees, necessary to do polarisation calibration.

    Arguments:
    ----------
    vis : str
        Input MeasurementSet (relative or absolute path).
    calfield : int
        Phase calibrator field ID.

    Returns:
    --------
    delta_parang : float
        The parallactic angle coverage of the phase calibrator field."""

    tb.open(vis+'::ANTENNA')
    pos = tb.getcol('POSITION')
    meanpos = np.mean(pos, axis=1)
    frame = tb.getcolkeyword('POSITION','MEASINFO')['Ref']
    units = tb.getcolkeyword('POSITION','QuantumUnits')
    mpos  = me.position(frame,
                    str(meanpos[0])+units[0],
                    str(meanpos[1])+units[1],
                    str(meanpos[2])+units[2])
    me.doframe(mpos)
    tb.close()

    # _geodetic_ latitude
    latr=me.measure(mpos,'WGS84')['m1']['value']
    tb.open(vis+'::FIELD')
    srcid = tb.getcol('SOURCE_ID')
    dirs=tb.getcol('DELAY_DIR')[:,0,:]
    tb.close()
    tb.open(vis,nomodify=True)
    st = tb.query('FIELD_ID=='+str(calfield))

    # get time stamps of first and last row
    nrows = st.nrows()
    tbeg = st.getcol('TIME', startrow=0, nrow=1)[0]
    tend = st.getcol('TIME', startrow=nrows-1, nrow=1)[0]

    # calculate parallactic angles for first and last time
    parang = np.zeros(2)

    # calculate parallactic angle
    rah = dirs[0,calfield]*12.0/np.pi
    decr = dirs[1,calfield]

    for itim, ts in enumerate([tbeg, tend]):
        tm = me.epoch('UTC',str(ts)+'s')
        last = me.measure(tm,'LAST')['m0']['value']
        last -= np.floor(last)  # days
        last *= 24.0  # hours
        ha = last-rah  # hours
        har = ha*2.0*np.pi/24.0
        parang[itim] = np.arctan2((np.cos(latr)*np.sin(har)), (np.sin(latr)*np.cos(decr) - np.cos(latr)*np.sin(decr)*np.cos(har)))

    delta_parang = np.rad2deg(parang[1] - parang[0])
    logger.debug("Delta parang: {0}".format(delta_parang))
    tb.close()

    return np.abs(delta_parang)


def get_xy_field(visname, fields, msmd=None):
    """
    From the input MS determine which field should
    be used for XY-phase calibration (if required).

    In the following order :
    3C286
    3C138
    secondaryfield (nominally dpolfield)
    """

    msmd.open(visname)
    fieldnames = msmd.fieldnames()
    msmd.done()

    # Use 3C286 or 3C138 if present in the data
    calibrator_3C286 = set(["3C286", "1328+307", "1331+305", "J1331+3030"]).intersection(set(fieldnames))
    calibrator_3C138 = set(["3C138", "0518+165", "0521+166", "J0521+1638"]).intersection(set(fieldnames))

    if calibrator_3C286:
        xyfield = list(calibrator_3C286)[0]
    elif calibrator_3C138:
        xyfield = list(calibrator_3C138)[0]
    else:
        xyfield = fields.dpolfield

    return xyfield