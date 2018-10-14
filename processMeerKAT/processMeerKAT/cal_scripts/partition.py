
"""
Runs partition on the input MS
"""
import sys

sys.path.append('/data/users/krishna/pipeline/processMeerKAT/processMeerKAT')

#from processMeerKAT import config_parser
import config_parser

# Get the name of the config file
args = config_parser.parse_args()
print(args)

# Parse config file
taskvals, config = config_parser.parse_config(args['config'])

# Partition
visname = taskvals['data']['vis']
mvis = visname.replace('.ms', '.mms')
partition(vis=visname, outputvis=mvis, createmms=True, datacolumn='DATA')
