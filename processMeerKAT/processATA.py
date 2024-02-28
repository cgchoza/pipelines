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
import config_parser
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

def parse_args():

    """Parse arguments into this script.

    Returns:
    --------
    args : class ``argparse.ArgumentParser``
        Known and validated arguments."""

    def parse_scripts(val):

        """Format individual arguments passed into a list for [ -S --scripts] argument, including paths and boolean values.

        Arguments/Returns:
        ------------------
        val : bool or str
            Path to script or container, or boolean representing whether that script is threadsafe (for MPI)."""

        if val.lower() in ('true','false'):
            return (val.lower() == 'true')
        else:
            return check_path(val)

    parser = argparse.ArgumentParser(prog=THIS_PROG,description='Process MeerKAT data via CASA MeasurementSet. Version: {0}'.format(__version__))

    parser.add_argument("-M","--MS",metavar="path", required=False, type=str, help="Path to MeasurementSet.")
    parser.add_argument("-C","--config",metavar="path", default=CONFIG, required=False, type=str, help="Relative (not absolute) path to config file.")
    parser.add_argument("-D","--plane", metavar="num", required=False, type=int, default=1,
                            help="Distribute tasks of this block size before moving onto next node [default: 1; max: ntasks-per-node].")
    parser.add_argument("-p","--partition", metavar="name", required=False, type=str, default="Main", help="SLURM partition to use [default: 'Main'].")
    parser.add_argument("-T","--time", metavar="time", required=False, type=str, default="12:00:00", help="Time limit to use for all jobs, in the form d-hh:mm:ss [default: '12:00:00'].")
    parser.add_argument("-S","--scripts", action='append', metavar='script', required=False, type=parse_scripts, default=SCRIPTS, help="Run pipeline with these scripts, in this order")
    parser.add_argument("-b","--precal_scripts", action='append', nargs=3, metavar=('script','threadsafe','container'), required=False, type=parse_scripts, default=PRECAL_SCRIPTS, help="Same as [-S --scripts], but run before calibration.")
    parser.add_argument("-a","--postcal_scripts", action='append', nargs=3, metavar=('script','threadsafe','container'), required=False, type=parse_scripts, default=POSTCAL_SCRIPTS, help="Same as [-S --scripts], but run after calibration.")
    parser.add_argument("-n","--name", metavar="unique", required=False, type=str, default='', help="Unique name to give this pipeline run (e.g. 'run1_'), appended to the start of all job names. [default: ''].")
    parser.add_argument("-d","--dependencies", metavar="list", required=False, type=str, default='', help="Comma-separated list (without spaces) of SLURM job dependencies (only used when nspw=1). [default: ''].")

    parser.add_argument("-l","--local", action="store_true", required=False, default=False, help="Build config file locally (i.e. without calling srun) [default: False].")
    parser.add_argument("-v","--verbose", action="store_true", required=False, default=False, help="Verbose output? [default: False].")
    parser.add_argument("-q","--quiet", action="store_true", required=False, default=False, help="Activate quiet mode, with suppressed output [default: False].")
    parser.add_argument("-P","--dopol", action="store_true", required=False, default=False, help="Perform polarization calibration in the pipeline [default: False].")
    parser.add_argument("-2","--do2GC", action="store_true", required=False, default=False, help="Perform (2GC) self-calibration in the pipeline [default: False].")
    parser.add_argument("-I","--science_image", action="store_true", required=False, default=False, help="Create a science image [default: False].")
    parser.add_argument("-x","--nofields", action="store_true", required=False, default=False, help="Do not read the input MS to extract field IDs [default: False].")
    parser.add_argument("-j","--justrun", action="store_true", required=False, default=False, help="Just run the pipeline, don't rebuild each job script if it exists [default: False].")

    #add mutually exclusive group - don't want to build config, run pipeline, or display version at same time
    run_args = parser.add_mutually_exclusive_group(required=True)
    run_args.add_argument("-B","--build", action="store_true", required=False, default=False, help="Build config file using input MS.")
    run_args.add_argument("-R","--run", action="store_true", required=False, default=False, help="Run pipeline with input config file.")
    run_args.add_argument("-V","--version", action="store_true", required=False, default=False, help="Display the version of this pipeline and quit.")
    run_args.add_argument("-L","--license", action="store_true", required=False, default=False, help="Display this program's license and quit.")

    args, unknown = parser.parse_known_args()

    if len(unknown) > 0:
        parser.error('Unknown input argument(s) present - {0}'.format(unknown))

    if args.run:
        if args.config is None:
            parser.error("You must input a config file [--config] to run the pipeline.")
        if not os.path.exists(args.config):
            parser.error("Input config file '{0}' not found. Please set [-C --config] or write a new one with [-B --build].".format(args.config))

    #if user inputs a list of scripts, remove the default list
    if len(args.scripts) > len(SCRIPTS):
        [args.scripts.pop(0) for i in range(len(SCRIPTS))]
    if len(args.precal_scripts) > len(PRECAL_SCRIPTS):
        [args.precal_scripts.pop(0) for i in range(len(PRECAL_SCRIPTS))]
    if len(args.postcal_scripts) > len(POSTCAL_SCRIPTS):
        [args.postcal_scripts.pop(0) for i in range(len(POSTCAL_SCRIPTS))]

    
    return args


# def write_master(filename,config,scripts=[],submit=False,dir='jobScripts',pad_length=5,verbose=False, echo=True, dependencies='',slurm_kwargs={}):

#     """Write master pipeline submission script, calling various sbatch files, and writing ancillary job scripts.

#     Arguments:
#     ----------
#     filename : str
#         Name of master pipeline submission script.
#     config : str
#         Path to config file.
#     scripts : list, optional
#         List of sbatch scripts to call in order.
#     submit : bool, optional
#         Submit jobs to SLURM queue immediately?
#     dir : str, optional
#         Name of directory to output ancillary job scripts.
#     pad_length : int, optional
#         Length to pad the SLURM sacct output columns.
#     verbose : bool, optional
#         Verbose output (inserted into master script)?
#     echo : bool, optional
#         Echo the pupose of each job script for the user?
#     dependencies : str, optional
#         Comma-separated list of SLURM job dependencies.
#     slurm_kwargs : list, optional
#         Parameters parsed from [slurm] section of config."""

#     master = open(filename,'w')
#     master.write('#!/bin/bash\n')
#     timestamp = config_parser.get_key(config,'run','timestamp')
#     if timestamp == '':
#         timestamp = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
#         config_parser.overwrite_config(config, conf_dict={'timestamp' : "'{0}'".format(timestamp)}, conf_sec='run', sec_comment='# Internal variables for pipeline execution')

#     #Copy config file to TMP_CONFIG and inform user
#     if verbose:
#         master.write("\necho Copying \'{0}\' to \'{1}\', and using this to run pipeline.\n".format(config,TMP_CONFIG))
#     master.write('cp {0} {1}\n'.format(config, TMP_CONFIG))

#     #Hack to perform correct number of selfcal loops
#     if config_parser.has_section(config,'selfcal') and 'selfcal_part1.sbatch' in scripts and 'selfcal_part2.sbatch' in scripts:
#         start_loop = config_parser.get_key(config, 'selfcal', 'loop')
#         selfcal_loops = config_parser.get_key(config, 'selfcal', 'nloops') - start_loop
#         idx = scripts.index('selfcal_part2.sbatch')

#         #check that we're doing nloops in order, otherwise don't duplicate scripts
#         if idx == scripts.index('selfcal_part1.sbatch') + 1:
#             init_scripts = scripts[:idx+1]
#             final_scripts = scripts[idx+1:]
#             init_scripts.extend(['selfcal_part1.sbatch','selfcal_part2.sbatch']*(selfcal_loops-1))
#             init_scripts.append('selfcal_part1.sbatch')
#             scripts = init_scripts + final_scripts

#     command = 'sbatch'

#     if dependencies != '':
#         master.write('\n#Run after these dependencies\nDep={0}\n'.format(dependencies))
#         command += " -d afterok:${Dep//,/:} --kill-on-invalid-dep=yes"
#     master.write('\n#{0}\n'.format(scripts[0]))
#     if verbose:
#         master.write('echo Submitting {0} to SLURM queue with following command:\necho {1} {0}.\n'.format(scripts[0],command))
#     master.write("IDs=$({0} {1} | cut -d ' ' -f4)\n".format(command,scripts[0]))
#     scripts.pop(0)


#     #Submit each script with dependency on all previous scripts, and extract job IDs
#     for script in scripts:
#         command = "sbatch -d afterok:${IDs//,/:} --kill-on-invalid-dep=yes"
#         master.write('\n#{0}\n'.format(script))
#         if verbose:
#             master.write('echo Submitting {0} to SLURM queue with following command\necho {1} {0}.\n'.format(script,command))
#         master.write("IDs+=,$({0} {1} | cut -d ' ' -f4)\n".format(command,script))

#     master.write('\n#Output message and create {0} directory\n'.format(dir))
#     master.write('echo Submitted sbatch jobs with following IDs: $IDs\n') #DON'T CHANGE as this output is relied on by bash sed expression in write_spw_master()
#     master.write('mkdir -p {0}\n'.format(dir))

#     #Add time as extn to this pipeline run, to give unique filenames
#     master.write('\n#Add time as extn to this pipeline run, to give unique filenames')
#     master.write("\nDATE={0}".format(timestamp))
#     extn = '_$DATE.sh'

#     #Copy contents of config file to jobScripts directory
#     master.write('\n#Copy contents of config file to {0} directory\n'.format(dir))
#     master.write('cp {0} {1}/{2}_$DATE.txt\n'.format(config,dir,os.path.splitext(config)[0]))

#     #Write each job script - kill script, summary script, error script, and timing script
#     write_all_bash_jobs_scripts(master,extn,IDs='IDs',dir=dir,echo=echo,pad_length=pad_length,slurm_kwargs=slurm_kwargs)

#     #Close master submission script and make executable
#     master.close()
#     os.chmod(filename, 509)

#     #Submit script or output that it will not run
#     if submit:
#         if echo:
#             logger.info('Running master script "{0}"'.format(filename))
#         os.system('./{0}'.format(filename))
#     else:
#         logger.info('Master script "{0}" written in "{1}", but will not run.'.format(filename,os.path.split(os.getcwd())[-1]))



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
