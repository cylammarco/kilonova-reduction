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
    data_log_file = config['data_log_file']
    data_log_level = config['data_log_level']
    reduction_log_file = config['reduction_log_file']
    reduction_log_level = config['reduction_log_level']
    propid_list = config['propid']
    db_name = config['db_name']
    obslog_table = config['obslog_table']
    data_table = config['data_table']
    reduction_table = config['reduction_table']
    raw_data_folder = config['raw_data_folder']
    reduced_data_folder = config['reduced_data_folder']
    time_start = config['time_start']
    time_end = config['time_end']
    login_config = config['login_config']

if base_folder[0] == '~':
    base_folder = os.path.join(expanduser("~"), base_folder[2:])

# get the base paths
if not os.path.isabs(base_folder):
    # if the script is called
    try:
        base_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   base_folder)
    # if run interactively
    except Exception as e:
        print('This script should not be run interactively for deployment. '
              'log level is set to DEBUG.')
        log_level = 'DEBUG'
        logging.basicConfig(level=logging.DEBUG)
        base_folder = os.path.join(os.path.abspath(os.getcwd()), base_folder)

if not os.path.exists(base_folder):
    mkdir_recursive(base_folder)

db_path = os.path.join(base_folder, db_name)
raw_data_folder_path = os.path.join(base_folder, raw_data_folder)
reduced_data_folder_path = os.path.join(base_folder, reduced_data_folder)
log_file_path = os.path.join(base_folder, data_log_file)

username = []
password = []
for login_config_path in login_config:

    with open(login_config_path) as login_f:

        login_config = json.load(login_f)
        username.append(decode(key, login_config['username']))
        password.append(decode(key, login_config['password']))

# Prepare the log format for ingestion
obslog_column_key = ('utc', 'name', 'propid', 'ra', 'dec', 'airmass',
                     'instrument', 'filter', 'binning', 'grating', 'exptime',
                     'seeing', 'sky', 'filename', 'groupid', 'err', 'qa')
data_column_key = ('filename', 'url', 'localpath', 'frametype')
reduction_column_key = ('filename', 'arc', 'sensitivity', 'stacked', 'output')

# Set-up logger
logger = logging.getLogger()
if data_log_level == "CRITICAL":
    logging.basicConfig(level=logging.CRITICAL)
if data_log_level == "ERROR":
    logging.basicConfig(level=logging.ERROR)
if data_log_level == "WARNING":
    logging.basicConfig(level=logging.WARNING)
if data_log_level == "INFO":
    logging.basicConfig(level=logging.INFO)
if data_log_level == "DEBUG":
    logging.basicConfig(level=logging.DEBUG)
formatter = logging.Formatter(
    '[%(asctime)s] %(levelname)s [%(filename)s:%(lineno)d] %(message)s',
    datefmt='%a, %d %b %Y %H:%M:%S')

# Configure to log into file
fh = logging.FileHandler(log_file_path, 'a+')
fh.setFormatter(formatter)
logger.addHandler(fh)

logging.debug("Obslog table columns: {}".format(obslog_column_key))
logging.debug("Raw data table columns: {}".format(data_column_key))
logging.debug("Reduced data table columns: {}".format(reduction_column_key))

logging.debug("The data grabbing starts at {} hours.".format(time_start))
logging.debug("The data grabbing stops at {} hours.".format(time_end))

# Connect to the database
db_connector = sqlite3.connect(db_path)
logging.info('Connected to the {}.'.format(db_name))

# Create the DB tables if not exist
check_obslog_table_exist_query = 'CREATE TABLE IF NOT EXISTS "{}" ({}'.format(
    obslog_table, ''.join(['"' + i + '" TEXT, '
                           for i in obslog_column_key]))[:-2] + ');'
logging.debug(check_obslog_table_exist_query)
db_connector.execute(check_obslog_table_exist_query)

check_data_table_exist_query = 'CREATE TABLE IF NOT EXISTS "{}" ({}'.format(
    data_table, ''.join(['"' + i + '" TEXT, '
                         for i in data_column_key]))[:-2] + ');'
logging.debug(check_data_table_exist_query)
db_connector.execute(check_data_table_exist_query)

check_reduction_table_exist_query = 'CREATE TABLE IF NOT EXISTS "{}" ({}'.format(
    reduction_table, ''.join(
        ['"' + i + '" TEXT, ' for i in reduction_column_key]))[:-2] + ');'
logging.debug(check_reduction_table_exist_query)
db_connector.execute(check_reduction_table_exist_query)

# Create the unique column if not exist
db_connector.execute(
    'CREATE UNIQUE INDEX IF NOT EXISTS filename_index on {}("filename");'.
    format(obslog_table))
db_connector.execute(
    'CREATE UNIQUE INDEX IF NOT EXISTS filename_index on {}("filename");'.
    format(data_table))
db_connector.execute(
    'CREATE UNIQUE INDEX IF NOT EXISTS filename_index on {}("filename");'.
    format(reduction_table))

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

    day = '24'
    day_m8 = '24'

    something_added = False

    # If the current time is between 18 hours and 7 hours, check if new data
    # is collected
    if (hour >= time_start) or (hour <= time_end):

        uri = 'https://telescope.livjm.ac.uk/DataProd/quicklook'

        for nth_prop, propid in enumerate(propid_list):

            ltobslog_filename = str(year_m8) + str(month_m8) + str(day_m8)
            obs_date = str(year) + '-' + str(month) + '-' + str(day)

            url = uri + '/{}/{}/{}.log'.format(propid, ltobslog_filename,
                                               ltobslog_filename)
            logging.info('Checking {}.'.format(url))
            u = url_is_alive(url, username[nth_prop], password[nth_prop])
            logging.info('HTTP status: {}.'.format(u.status_code))

            # If observing log exists, get it
            if u.status_code == 200:

                # Load text into array
                obslog = np.genfromtxt(io.StringIO(u.text),
                                       skip_header=4,
                                       dtype=None,
                                       names=obslog_column_key,
                                       encoding=None)
                logging.debug(obslog)

                # Create folders to hold raw and reduced data if not exist
                if not os.path.exists(
                        os.path.join(raw_data_folder_path,
                                        ltobslog_filename)):
                    mkdir_recursive(
                        os.path.join(raw_data_folder_path,
                                        ltobslog_filename))

                if not os.path.exists(
                        os.path.join(reduced_data_folder_path,
                                    ltobslog_filename)):
                    mkdir_recursive(
                        os.path.join(reduced_data_folder_path,
                                    ltobslog_filename))

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
                        insertion_query = 'INSERT INTO "{}" {} VALUES {};'.format(
                            obslog_table, obslog_column_key, new_row)
                        logging.debug(insertion_query)
                        db_connector.execute(insertion_query)
                        db_connector.commit()

                        # create download path
                        download_path = uri + '/{}/{}/{}{}'.format(
                            propid, ltobslog_filename, row['filename'],
                            '.fits.gz')

                        storage_path = os.path.join(
                            raw_data_folder_path, ltobslog_filename,
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

                            logging.info(row['filename'] +
                                         '.fits.gz downloaded.')
                            downloaded = "True"

                        except Exception as e:

                            logging.error(e)
                            downloaded = "False"

                        file_type_letter = row['filename'].split('_')[1]
                        if file_type_letter == 'a':
                            frame_type = 'arc'
                        elif file_type_letter == 'b':
                            frame_type = 'bias'
                        elif file_type_letter == 'd':
                            frame_type = 'dark'
                        elif file_type_letter == 'e':
                            frame_type = 'science'
                        elif file_type_letter == 'q':
                            frame_type = 'acquisition'
                        elif file_type_letter == 's':
                            frame_type = 'standard'
                        elif file_type_letter == 'w':
                            frame_type = 'flat'
                        else:
                            logging.error('Unknown frame type: {}.'.format(
                                file_type_letter))

                        # log into data table
                        data_insertion_query = 'INSERT INTO "{}" {} VALUES {};'.format(
                            data_table, data_column_key,
                            (row['filename'], download_path, storage_path,
                             frame_type))
                        logging.debug(data_insertion_query)
                        db_connector.execute(data_insertion_query)
                        db_connector.commit()

                        something_added = True

                        # When an arc is loaded, look for all the immediate preceeeding light frames and
                        # ingest them all into the reduction table
                        if frame_type == 'arc':

                            obs_night, frame_number = row['filename'].split(
                                '_')[2:4]

                            # Query for the preceeding frames
                            preceeding_light_frames_query = 'SELECT filename from {} WHERE filename LIKE "v\_%\_{}\_{}%" ESCAPE "\\";'.format(
                                data_table, obs_night,
                                str(int(frame_number) - 1))
                            logging.debug(preceeding_light_frames_query)
                            light_frame_list = list(
                                np.array(
                                    db_connector.execute(
                                        preceeding_light_frames_query).
                                    fetchall()).flatten())

                            for light_frame_filename in light_frame_list:
                                # log into reduction table
                                logging.debug(
                                    'light frame: {} is being ingested to table {}.'
                                    .format(light_frame_filename,
                                            reduction_table))
                                reduction_insertion_query = 'INSERT INTO "{}" {} VALUES {};'.format(
                                    reduction_table, reduction_column_key,
                                    (light_frame_filename, row['filename'],
                                     'UNKNOWN', 'False', 'UNKNOWN'))
                                logging.debug(reduction_insertion_query)
                                db_connector.execute(reduction_insertion_query)
                                db_connector.commit()

                    # If the row already exists in the DB, do nothing
                    else:

                        pass

                if not something_added:
                    logging.info('Nothing new.')

            else:

                # do nothing
                pass

            # sleep for 5 minutes
            logging.info('Sleeping for 5 minutes.')
            time.sleep(5 * 60)

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
