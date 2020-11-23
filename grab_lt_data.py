import io
import json
import getpass
import logging
import numpy as np
import os
import requests
import sqlite3
import sys
import time
import urllib
from datetime import datetime, timedelta, tzinfo
from os.path import expanduser

from encode import decode

key = getpass.getpass(prompt='key: ')


class GMTm8(tzinfo):
    def utcoffset(self, dt):
        return timedelta(hours=-8)

    def dst(self, dt):
        # Ignore daylight saving
        return timedelta(0)

    def tzname(self, dt):
        return "GMT-8"

def mkdir_recursive(path):
    sub_path = os.path.dirname(path)
    if not os.path.exists(sub_path):
        mkdir_recursive(sub_path)
    if not os.path.exists(path):
        os.mkdir(path)

def url_is_alive(url, username, password):
    """
    Checks that a given URL is reachable.

    Parameter
    ---------
    url:
        A URL

    Return
    ------
    bool

    """

    r = requests.get(url, auth=(username, password))

    return r


# Load the configuration file
with open('config.json') as f:
    config = json.load(f)
    base_folder = config['base_folder']
    log_file = config['log_file']
    log_level = config['log_level']
    db_name = config['db_name']
    obslog_table = config['obslog_table']
    data_table = config['data_table']
    raw_data_folder = config['raw_data_folder']
    time_start = config['time_start']
    time_end = config['time_end']
    login_config = config['login_config']

with open('login_cfg.json') as login_f:
    login_config = json.load(login_f)
    username = decode(key, login_config['username'])
    password = decode(key, login_config['password'])

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
obslog_column_key = ('utc', 'name', 'propid', 'ra', 'dec', 'airmass',
                     'instrument', 'filter', 'binning', 'grating', 'exptime',
                     'seeing', 'sky', 'filename', 'groupid', 'err', 'qa')
data_column_key = ('filename', 'remote', 'local', 'downloaded')

# operating hours
time_start = 18
time_end = 7

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

# Create the db table if not exist
check_obslog_exist_query = 'CREATE TABLE IF NOT EXISTS "{}" ({}'.format(
    obslog_table, ''.join(['"' + i + '" TEXT, '
                           for i in obslog_column_key]))[:-2] + ');'
logging.debug(check_obslog_exist_query)
db_connector.execute(check_obslog_exist_query)

check_data_exist_query = 'CREATE TABLE IF NOT EXISTS "{}" ({}'.format(
    data_table, ''.join(['"' + i + '" TEXT, '
                         for i in data_column_key]))[:-2] + ');'
logging.debug(check_data_exist_query)
db_connector.execute(check_data_exist_query)

# Create the unique column if not exist
db_connector.execute(
    'CREATE UNIQUE INDEX IF NOT EXISTS filename_index on {}("filename");'.
    format(obslog_table))
db_connector.execute(
    'CREATE UNIQUE INDEX IF NOT EXISTS filename_index on {}("filename");'.
    format(data_table))

# run indefinitely until being terminated externally
while 1:

    # Get the current UTC and UTC-8. The latter is because the filename is
    # named after the date of the beginning of the night, GMT-8 guarantee the
    # name to be correct up to morning twilight.
    time_now = datetime.utcnow()
    time_now_m8 = time_now.astimezone(GMTm8())
    logging.info('UTC now is {}.'.format(time_now))
    logging.info('UTC-8 now is {}.'.format(time_now_m8))

    year = time_now.year
    month = time_now.month
    day = time_now.day
    hour = time_now.hour
    minute = time_now.minute

    year_m8 = time_now_m8.year
    month_m8 = time_now_m8.month
    day_m8 = time_now_m8.day

    something_added = False

    # If the current time is between 18 hours and 7 hours, check if new data
    # is collected
    if (hour >= time_start) or (hour <= time_end):

        uri = 'https://telescope.livjm.ac.uk/DataProd/quicklook'
        propid = 'SPRATBias'
        ltobslog_filename = str(year_m8) + str(month_m8) + str(day_m8)
        obs_date = str(year) + '-' + str(month) + '-' + str(day)

        url = uri + '/{}/{}/{}.log'.format(propid, ltobslog_filename,
                                           ltobslog_filename)
        logging.info('Checking {}.'.format(url))
        u = url_is_alive(url, username, password)
        logging.info('HTTP status: {}.'.format(u.status_code))

        # If observing log exists, get it
        if u.status_code == 200:

            # Load text into array
            obslog = np.genfromtxt(io.StringIO(u.text),
                                   skip_header=4,
                                   dtype=None,
                                   names=obslog_column_key,
                                   encoding=None)

            for row in obslog:

                check_if_exist_query = 'SELECT * FROM {} WHERE "filename" = "{}";'.format(
                    obslog_table, row['filename'])

                # Add to the DB if not exist
                if db_connector.execute(
                        check_if_exist_query).fetchone() is None:

                    new_row = row.flatten()[0]
                    logging.info('Ingesting to {}.'.format(db_name))

                    # log the UTC in the obslog as info
                    # log each value in the obslog as debug
                    for i, value in enumerate(new_row):

                        # Change the time to include date
                        if i == 0:
                            new_row[i] = obs_date + 'T' + new_row[i]

                        if i == 13:
                            logging.info('    {}: {}'.format(
                                obslog_column_key[i], value))
                        else:
                            logging.debug('    {}: {}'.format(
                                obslog_column_key[i], value))

                    # log into the obslog table
                    insertion_query = 'INSERT INTO "{}" {} VALUES {}'.format(
                        obslog_table, obslog_column_key, new_row)
                    logging.debug(insertion_query)
                    db_connector.execute(insertion_query)

                    # create download path
                    download_path = uri + '/{}/{}/{}{}'.format(
                        propid, ltobslog_filename, row['filename'], '.fits.gz')
                    if not os.path.exists(
                            os.path.join(raw_data_folder_path,
                                         ltobslog_filename)):
                        mkdir_recursive(
                            os.path.join(raw_data_folder_path,
                                         ltobslog_filename))

                    storage_path = os.path.join(raw_data_folder_path,
                                                ltobslog_filename,
                                                download_path.split('/')[-1])
                    logging.info(download_path)
                    logging.info(storage_path)
                    # Save the downloaded file
                    try:

                        with urllib.request.urlopen(
                                download_path) as response, open(
                                    storage_path, 'wb') as out_file:

                            data = response.read()  # a `bytes` object
                            out_file.write(data)

                        logging.info(row['filename'] + '.fits.gz downloaded.')
                        downloaded = "True"

                    except Exception as e:

                        logging.error(e)
                        downloaded = "False"

                    # log into data table
                    insertion_query = 'INSERT INTO "{}" {} VALUES {}'.format(
                        data_table, data_column_key,
                        (row['filename'], download_path, storage_path,
                         downloaded))
                    logging.debug(insertion_query)
                    db_connector.execute(insertion_query)

                    something_added = True

                # If the row already exists in the DB, do nothing
                else:

                    pass

            if not something_added:
                logging.info('Nothing new.')

        else:

            # do nothing
            pass

        # sleep for 15 minutes
        logging.info('Sleeping for 15 minutes.')
        time.sleep(15 * 60)

    else:

        time_diff = time_start - hour - minute / 60.
        logging.info('Sleeping for {} hours.'.format(time_diff))

        while 1:

            utc_hour_now = datetime.utcnow().hour
            if (utc_hour_now > time_end) and (utc_hour_now < time_start):

                time.sleep(5 * 60)

            else:

                break

'''
if response.status_code != 200:
    raise ValueError(
        'Request to slack returned an error %s, the response is:\n%s'
        % (response.status_code, response.text)
)
'''
