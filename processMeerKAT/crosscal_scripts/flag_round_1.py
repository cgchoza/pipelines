#Copyright (C) 2022 Inter-University Institute for Data Intensive Astronomy
#See processMeerKAT.py for license details.

import sys
import os

import config_parser
import bookkeeping
# sys.path.append('/home/cchoza/pipelines/processMeerKAT/')
from crosscal_scripts.config import CONFIG_PATH

def do_pre_flag(visname, fields, badfreqranges, badants):

    clip = [0., 50.]

    if badfreqranges!='[]':
        badspw = '*:' + ',*:'.join(badfreqranges)
        flagdata(vis=visname, mode='manual', spw=badspw)

    if badants != '[]':
        badants = badants.split("#")[0].strip(" ").strip("[").strip("]")
        flagdata(vis=visname, mode='manual', antenna=badants)

    flagdata(vis=visname, mode='manual', autocorr=True, action='apply',
            flagbackup=True, savepars=False, writeflags=True)

    #Manually clip all fields in config
    allfields = ','.join(set([i for i in (','.join([fields.gainfields] + [fields.targetfield] + [fields.extrafields]).split(',')) if i])) #remove duplicate and empty fields

    flagdata(vis=visname, mode="clip", field=allfields,
            clipminmax=clip, datacolumn="DATA",clipoutside=True,
            clipzeros=True, extendpols=True, action="apply",flagbackup=True,
            savepars=False, overwrite=True, writeflags=True)
    
    #tfcrop calfields and targetfield with different cutoffs / fits
    calfields = ','.join(set([i for i in (','.join([fields.gainfields] + [fields.extrafields]).split(',')) if i])) #remove duplicate and empty fields

    flagdata(vis=visname, mode='tfcrop', field=calfields,
            ntime='scan', timecutoff=5.0, freqcutoff=5.0, timefit='line',
            freqfit='line', extendflags=False, timedevscale=5., freqdevscale=5.,
            extendpols=True, growaround=False, action='apply', flagbackup=True,
            overwrite=True, writeflags=True, datacolumn='DATA')

    flagdata(vis=visname, mode='tfcrop', field=fields.targetfield,
            ntime='scan', timecutoff=6.0, freqcutoff=6.0, timefit='poly',
            freqfit='poly', extendflags=False, timedevscale=5., freqdevscale=5.,
            extendpols=True, growaround=False, action='apply', flagbackup=True,
            overwrite=True, writeflags=True, datacolumn='DATA')

    # Conservatively extend flags for all fields in config
    flagdata(vis=visname, mode='extend', field=allfields,
            datacolumn='data', clipzeros=True, ntime='scan', extendflags=False,
            extendpols=True, growtime=80., growfreq=80., growaround=False,
            flagneartime=False, flagnearfreq=False, action='apply',
            flagbackup=True, overwrite=True, writeflags=True)

    flagdata(vis=visname, mode='summary', datacolumn='DATA',
            name=visname+'.flag.summary')



########### RUN IT DOWN HERE ##############

taskvals,config = config_parser.parse_config(filename=CONFIG_PATH)
visname = config['data']['vis'].strip("'")

calfiles, caldir = bookkeeping.bookkeeping(visname)
fields = bookkeeping.get_field_ids(config['fields'])

badfreqranges = config['crosscal']['badfreqranges']
badants = config['crosscal']['badants']

do_pre_flag(visname, fields, badfreqranges, badants)
