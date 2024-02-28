#!/usr/bin/env python3

# Removing all SLURM arguments from this script! Will explore other necessary scripts and remove SLURM as needed
__version__ = 'ATA0.0'

license = """
    Process MeerKAT data via CASA MeasurementSet.
    Copyright (C) 2022 Inter-University Institute for Data Intensive Astronomy.
    support@ilifu.ac.za

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import os
import sys
import re
import config_parser  # Circular import??
import bookkeeping
from shutil import copyfile
from copy import deepcopy
import logging
from time import gmtime
from datetime import datetime
logging.Formatter.converter = gmtime
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s: %(message)s")

# #Set global limits for current ilifu cluster configuration
# TOTAL_NODES_LIMIT = 79
# CPUS_PER_NODE_LIMIT = 32
# NTASKS_PER_NODE_LIMIT = CPUS_PER_NODE_LIMIT
# MEM_PER_NODE_GB_LIMIT = 232 #237568 MB
# MEM_PER_NODE_GB_LIMIT_HIGHMEM = 480 #491520 MB

#Set global values for paths and file names
THIS_PROG = __file__
SCRIPT_DIR = os.path.dirname(THIS_PROG)
LOG_DIR = 'logs'
CALIB_SCRIPTS_DIR = 'crosscal_scripts'
AUX_SCRIPTS_DIR = 'aux_scripts'
SELFCAL_SCRIPTS_DIR = 'selfcal_scripts'
CONFIG = 'default_config.txt'
TMP_CONFIG = '.config.tmp'
# MASTER_SCRIPT = 'submit_pipeline.sh'
SPW_PREFIX = '*:'

#Set global values for field, crosscal and SLURM arguments copied to config file, and some of their default values
FIELDS_CONFIG_KEYS = ['fluxfield','bpassfield','phasecalfield','targetfields','extrafields']
CROSSCAL_CONFIG_KEYS = ['minbaselines','chanbin','width','timeavg','createmms','keepmms','spw','nspw','calcrefant','refant','standard','badants','badfreqranges']
SELFCAL_CONFIG_KEYS = ['nloops','loop','cell','robust','imsize','wprojplanes','niter','threshold','uvrange','nterms','gridder','deconvolver','solint','calmode','discard_nloops','gaintype','outlier_threshold','flag','outlier_radius']
IMAGING_CONFIG_KEYS = ['cell', 'robust', 'imsize', 'wprojplanes', 'niter', 'threshold', 'multiscale', 'nterms', 'gridder', 'deconvolver', 'restoringbeam', 'stokes', 'mask', 'rmsmap','outlierfile', 'pbthreshold', 'pbband']
# SLURM_CONFIG_STR_KEYS = ['container','mpi_wrapper','partition','time','name','dependencies','exclude','account','reservation']
# SLURM_CONFIG_KEYS = ['nodes','ntasks_per_node','mem','plane','submit','precal_scripts','postcal_scripts','scripts','verbose','modules'] + SLURM_CONFIG_STR_KEYS
# CONTAINER = '/idia/software/containers/casa-6.5.0-modular.sif'
# MPI_WRAPPER = 'mpirun'
RUN_CONFIG_KEYS = ['verbose', 'scripts']
PRECAL_SCRIPTS = [('calc_refant.py',False,''),('partition.py',True,'')] #Scripts run before calibration at top level directory when nspw > 1
POSTCAL_SCRIPTS = [('concat.py',False,''),('plotcal_spw.py', False, ''),('selfcal_part1.py',True,''),('selfcal_part2.py',False,''),('science_image.py', True, '')] #Scripts run after calibration at top level directory when nspw > 1
SCRIPTS = [ ('validate_input.py',False,''),
            ('flag_round_1.py',True,''),
            ('calc_refant.py',False,''),
            ('setjy.py',True,''),
            ('xx_yy_solve.py',False,''),
            ('xx_yy_apply.py',True,''),
            ('flag_round_2.py',True,''),
            ('xx_yy_solve.py',False,''),
            ('xx_yy_apply.py',True,''),
            ('split.py',True,''),
            ('quick_tclean.py',True,'')]


def check_path(path,update=False):

    """Check in specific location for a script or container, including in bash path, and in this pipeline's calibration
    scripts directory (SCRIPT_DIR/{CALIB_SCRIPTS_DIR,AUX_SCRIPTS_DIR}/). If path isn't found, raise IOError, otherwise return the path.

    Arguments:
    ----------
    path : str
        Check for script or container at this path.
    update : bool, optional
        Update the path according to where the file is found.

    Returns:
    --------
    path : str
        Path to script or container (if path found and update=True)."""

    newpath = path

    #Attempt to find path firstly in CWD, then directory up, then pipeline directories, then bash path.
    if os.path.exists(path) and path[0] != '/':
        newpath = '{0}/{1}'.format(os.getcwd(),path)
    if not os.path.exists(path) and path != '':
        if os.path.exists('../{0}'.format(path)):
            newpath = '../{0}'.format(path)
        elif os.path.exists('{0}/{1}'.format(SCRIPT_DIR,path)):
            newpath = '{0}/{1}'.format(SCRIPT_DIR,path)
        elif os.path.exists('{0}/{1}/{2}'.format(SCRIPT_DIR,CALIB_SCRIPTS_DIR,path)):
            newpath = '{0}/{1}/{2}'.format(SCRIPT_DIR,CALIB_SCRIPTS_DIR,path)
        elif os.path.exists('{0}/{1}/{2}'.format(SCRIPT_DIR,AUX_SCRIPTS_DIR,path)):
            newpath = '{0}/{1}/{2}'.format(SCRIPT_DIR,AUX_SCRIPTS_DIR,path)
        elif os.path.exists('{0}/{1}/{2}'.format(SCRIPT_DIR,SELFCAL_SCRIPTS_DIR,path)):
            newpath = '{0}/{1}/{2}'.format(SCRIPT_DIR,SELFCAL_SCRIPTS_DIR,path)
        elif os.path.exists(check_bash_path(path)):
            newpath = check_bash_path(path)
        else:
            #If it still doesn't exist, throw error
            raise IOError('File "{0}" not found.'.format(path))

    if update:
        return newpath
    else:
        return path

def check_bash_path(fname):

    """Check if file is in your bash path and executable (i.e. executable from command line), and prepend path to it if so.

    Arguments:
    ----------
    fname : str
        Filename to check.

    Returns:
    --------
    fname : str
        Potentially updated filename with absolute path prepended."""

    PATH = os.environ['PATH'].split(':')
    for path in PATH:
        if os.path.exists('{0}/{1}'.format(path,fname)):
            if not os.access('{0}/{1}'.format(path,fname), os.X_OK):
                raise IOError('"{0}" found in "{1}" but file is not executable.'.format(fname,path))
            else:
                fname = '{0}/{1}'.format(path,fname)
            break

    return fname


def get_config_kwargs(config,section,expected_keys):

    """Return kwargs from config section. Check section exists, and that all expected keys are present, otherwise raise KeyError.

    Arguments:
    ----------
    config : str
        Path to config file.
    section : str
        Config section from which to extract kwargs.
    expected_keys : list
        List of expected keys.

    Returns:
    --------
    kwargs : dict
        Keyword arguments from this config section."""

    config_dict = config_parser.parse_config(config)[0]

    #Ensure section exists, otherwise raise KeyError
    if section not in config_dict.keys():
        raise KeyError("Config file '{0}' has no section [{1}]. Please insert section or build new config with [-B --build].".format(config,section))

    kwargs = config_dict[section]

    #Check for any unknown keys and display warning
    unknown_keys = list(set(kwargs) - set(expected_keys))
    if len(unknown_keys) > 0:
        logger.warning("Unknown keys {0} present in section [{1}] in '{2}'.".format(unknown_keys,section,config))

    #Check that expected keys are present, otherwise raise KeyError
    missing_keys = list(set(expected_keys) - set(kwargs))
    if len(missing_keys) > 0:
        raise KeyError("Keys {0} missing from section [{1}] in '{2}'. Please add these keywords to '{2}', or else run [-B --build] step again.".format(missing_keys,section,config))

    return kwargs

def setup_logger(config,verbose=False):

    """Setup logger at debug or info level according to whether verbose option selected (via command line or config file).

    Arguments:
    ----------
    config : str
        Path to config file.
    verbose : bool
        Verbose output? This will display all logger debug output."""

    #Overwrite with verbose mode if set to True in config file
    if not verbose:
        config_dict = config_parser.parse_config(config)[0]
        if 'run' in config_dict.keys() and 'run' in config_dict['slurm']:
            verbose = config_dict['run']['verbose']

    loglevel = logging.DEBUG if verbose else logging.INFO
    logger.setLevel(loglevel)

def main():

    # We're turning this into a run-engine instead

    #Parse command-line arguments, and setup logger
    args = parse_args()

    #Mutually exclusive arguments - display version, build config file or run pipeline
    if args.version:
        logger.info('This is version {0}'.format(__version__))
    if args.license:
        logger.info(license)
    


if __name__ == "__main__":
    main()
