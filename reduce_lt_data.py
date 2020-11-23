from astropy.io import fits
from aspired import image_reduction
from aspired import spectral_reduction
import json
import numpy as np
from scipy import interpolate as itp


# Load the configuration file
with open('reduction_config.json') as f:
    config = json.load(f)
    base_folder = config['base_folder']
    log_file = config['log_file']
    log_level = config['log_level']
    db_name = config['db_name']
    data_table = config['data_table']
    reduction_table = config['reduction_table']
    raw_data_folder = config['raw_data_folder']
    reduced_data_folder = config['reduced_data_folder']
    slack_config = config['slack_config']

with open('slack_cfg.json') as slack_f:
    slack_config = json.load(slack_f)
    webhook_url = decode(key, slack_config['webhook_url'])
    slack_token = decode(key, slack_config['token'])


# Set-up logger
logger = logging.getLogger()
if log_level == "CRITICAL":
    logging.basicConfig(level=logging.CRITICAL)
if log_level == "ERROR":
    logging.basicConfig(level=logging.ERROR)
if log_level == "WARNING":
    logging.basicConfig(level=logging.WARNING)
if log_level == "INFO":
    logging.basicConfig(level=logging.INFO)
if log_level == "DEBUG":
    logging.basicConfig(level=logging.DEBUG)
formatter = logging.Formatter(
    '[%(asctime)s] %(levelname)s [%(filename)s:%(lineno)d] %(message)s',
    datefmt='%a, %d %b %Y %H:%M:%S')

# Prepare the log format for ingestion
data_column_key = ('filename', 'remote', 'local', 'downloaded')
reduction_column_key = ('filename', 'arc', 'sensitivity', 'stacked', 'output')

if base_folder[0] == '~':
    base_folder = os.path.join(expanduser("~"), base_folder[2:])

# get the base paths
if not os.path.isabs(base_folder):
    # if the script is called
    try:
        base_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   base_folder)
    # if run interactively
    except:
        base_folder = os.path.join(os.path.abspath(os.getcwd()), base_folder)

db_path = os.path.join(base_folder, db_name)
raw_data_folder_path = os.path.join(base_folder, raw_data_folder)
reduced_data_folder_path = os.path.join(base_folder, reduced_data_folder)
log_file_path = os.path.join(base_folder, log_file)

# Configure to log into file
fh = logging.FileHandler(log_file_path, 'a+')
fh.setFormatter(formatter)
logger.addHandler(fh)

# Connect to the database
db_connector = sqlite3.connect(db_path)
logging.info('Connected to the {}.'.format(db_name))

# Line list
atlas = [
    4193.5, 4385.77, 4500.98, 4524.68, 4582.75, 4624.28, 4671.23, 4697.02,
    4734.15, 4807.02, 4921.48, 5028.28, 5618.88, 5823.89, 5893.29, 5934.17,
    6182.42, 6318.06, 6472.841, 6595.56, 6668.92, 6728.01, 6827.32, 6976.18,
    7119.60, 7257.9, 7393.8, 7584.68, 7642.02, 7740.31, 7802.65, 7887.40,
    7967.34, 8057.258
]
element = ['Xe'] * len(atlas)

while 1:

