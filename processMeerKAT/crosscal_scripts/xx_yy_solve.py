#Copyright (C) 2022 Inter-University Institute for Data Intensive Astronomy
#See processMeerKAT.py for license details.

import sys
import os
import shutil

import config_parser
import bookkeeping

from casatasks import *
logfile=casalog.logfile()
from crosscal_scripts.config import CONFIG_PATH

import logging
from time import gmtime
logging.Formatter.converter = gmtime
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s: %(message)s", level=logging.INFO)

def do_parallel_cal(visname, fields, calfiles, referenceant, caldir,
        minbaselines, standard):

    if not os.path.isdir(caldir):
        os.makedirs(caldir)
    elif not os.path.isdir(caldir+'_round1'):
        os.rename(caldir,caldir+'_round1')
        os.makedirs(caldir)

    logger.info(" starting antenna-based delay (kcorr)\n -> %s" % calfiles.kcorrfile)
    gaincal(vis=visname, caltable = calfiles.kcorrfile, field
            = fields.kcorrfield, refant = referenceant,
            minblperant = minbaselines, solnorm = False,  gaintype = 'K',
            solint = 'inf', combine = '', parang = False, append = False)
    bookkeeping.check_file(calfiles.kcorrfile)

    logger.info(" starting bandpass -> %s" % calfiles.bpassfile)
    bandpass(vis=visname, caltable = calfiles.bpassfile,
            field = fields.bpassfield, refant = referenceant,
            minblperant = minbaselines, solnorm = False,  solint = 'inf',
            combine = 'scan', bandtype = 'B', fillgaps = 8,
            gaintable = calfiles.kcorrfile, gainfield = fields.kcorrfield,
            parang = False, append = False)
    bookkeeping.check_file(calfiles.bpassfile)

    logger.info(" starting gain calibration\n -> %s" % calfiles.gainfile)
    gaincal(vis=visname, caltable = calfiles.gainfile,
            field = fields.gainfields, refant = referenceant,
            minblperant = minbaselines, solnorm = False,  gaintype = 'G',
            solint = 'inf', combine = '', calmode='ap',
            gaintable=[calfiles.kcorrfile, calfiles.bpassfile],
            gainfield=[fields.kcorrfield, fields.bpassfield],
            parang = False, append = False)
    bookkeeping.check_file(calfiles.gainfile)

    # Only run fluxscale if bootstrapping
    if len(fields.gainfields.split(',')) > 1:
        fluxscale(vis=visname, caltable=calfiles.gainfile,
                reference=[fields.fluxfield], transfer='',
                fluxtable=calfiles.fluxfile, append=False, display=False,
                listfile = os.path.join(caldir,'fluxscale_xx_yy.txt'))
        bookkeeping.check_file(calfiles.fluxfile)



taskvals,config = config_parser.parse_config(filename=CONFIG_PATH)
visname = config['data']['vis'].strip("'")

calfiles, caldir = bookkeeping.bookkeeping(visname)
fields = bookkeeping.get_field_ids(config['fields'])

minbaselines = taskvals['crosscal']['minbaselines']
standard = taskvals['crosscal']['standard']
refant = taskvals['crosscal']['refant']

do_parallel_cal(visname, fields, calfiles, f"'{refant}'", caldir, minbaselines, standard)
