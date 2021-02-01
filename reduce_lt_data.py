from astropy.io import fits
from aspired import spectral_reduction
import getpass
import logging
import json
import numpy as np
import os
import sqlite3
import time
from scipy import interpolate as itp
from datetime import datetime
from os.path import expanduser

from encode import decode
from slackcomm import post_message_to_slack, post_file_to_slack

key = getpass.getpass(prompt='key: ')

# Load the configuration file
with open('config.json') as f:
    config = json.load(f)
    base_folder = config['base_folder']
    data_log_file = config['data_log_file']
    data_log_level = config['data_log_level']
    reduction_log_file = config['reduction_log_file']
    reduction_log_level = config['reduction_log_level']
    db_name = config['db_name']
    data_table = config['data_table']
    reduction_table = config['reduction_table']
    raw_data_folder = config['raw_data_folder']
    reduced_data_folder = config['reduced_data_folder']
    time_start = config['time_start']
    time_end = config['time_end']
    slack_config = config['slack_config']

# Set-up logger
logger = logging.getLogger()
if reduction_log_level == "CRITICAL":
    logging.basicConfig(level=logging.CRITICAL)
if reduction_log_level == "ERROR":
    logging.basicConfig(level=logging.ERROR)
if reduction_log_level == "WARNING":
    logging.basicConfig(level=logging.WARNING)
if reduction_log_level == "INFO":
    logging.basicConfig(level=logging.INFO)
if reduction_log_level == "DEBUG":
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
    except Exception as e:
        print('This script should not be run interactively for deployment. '
              'log level is set to DEBUG.')
        print(e)
        log_level = 'DEBUG'
        logging.basicConfig(level=logging.DEBUG)
        base_folder = os.path.join(os.path.abspath(os.getcwd()), base_folder)

db_path = os.path.join(base_folder, db_name)
raw_data_folder_path = os.path.join(base_folder, raw_data_folder)
reduced_data_folder_path = os.path.join(base_folder, reduced_data_folder)
log_file_path = os.path.join(base_folder, reduction_log_file)

with open(slack_config) as slack_f:
    slack_cfg = json.load(slack_f)
    webhook_url = decode(key, slack_cfg['webhook_url'])
    slack_token = decode(key, slack_cfg['token'])

# Configure to log into file
fh = logging.FileHandler(log_file_path, 'a+')
fh.setFormatter(formatter)
logger.addHandler(fh)

# Connect to the database
db_connector = sqlite3.connect(db_path)
logging.info('Connected to the {}.'.format(db_name))

# Create the db table if not exist
check_reduction_exist_query = 'CREATE TABLE IF NOT EXISTS "{}" ({}'.format(
    reduction_table, ''.join(
        ['"' + i + '" TEXT, ' for i in reduction_column_key]))[:-2] + ');'
logging.debug(check_reduction_exist_query)
db_connector.execute(check_reduction_exist_query)

# Create the unique column if not exist
db_connector.execute(
    'CREATE UNIQUE INDEX IF NOT EXISTS filename_index on {} ("filename");'.
    format(reduction_table))

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

    # Get the current UTC and UTC-8. The latter is because the filename is
    # named after the date of the beginning of the night, GMT-8 guarantee the
    # name to be correct up to morning twilight.
    time_now = datetime.utcnow()
    logging.info('UTC now is {}.'.format(time_now))

    year = time_now.year
    month = time_now.month
    day = time_now.day
    hour = time_now.hour
    minute = time_now.minute

    logging.debug(
        db_connector.execute(
            "SELECT * FROM {} LIMIT 1;".format(reduction_table)).fetchall())

    # Check for data that is not reduced
    check_if_reduced_query =\
        'SELECT * FROM {} WHERE "output" = "UNKNOWN";'.format(reduction_table)
    logging.info(check_if_reduced_query)
    list_to_reduce = db_connector.execute(check_if_reduced_query).fetchall()
    logging.info(list_to_reduce)

    # Add to the DB if not exist
    if list_to_reduce == []:

        logging.info('Sleeping for 5 minutes.')
        time.sleep(5 * 60)

    else:

        for row in list_to_reduce:

            light_name = row[0].split('.')[0]
            arc_name = row[1].split('.')[0]
            obsnight = light_name.split('_')[2]
            extension = '.fits.gz'

            light_file_path = os.path.join(raw_data_folder_path, obsnight,
                                           light_name + extension)
            arc_file_path = os.path.join(raw_data_folder_path, obsnight,
                                         arc_name + extension)

            # Load the light frame
            twodspec = spectral_reduction.TwoDSpec(light_file_path,
                                                   cosmicray=True,
                                                   readnoise=5.7)

            twodspec.ap_trace(nspec=1,
                              nwindow=10,
                              display=False,
                              filename=os.path.join(reduced_data_folder_path,
                                                    obsnight,
                                                    light_name + '_aptrace'),
                              save_iframe=True)

            twodspec.ap_extract(apwidth=7,
                                skywidth=5,
                                skysep=3,
                                skydeg=1,
                                optimal=True,
                                display=False,
                                filename=os.path.join(
                                    reduced_data_folder_path, obsnight,
                                    light_name + '_apextract'),
                                save_iframe=True)

            twodspec.add_arc(arc_file_path)

            twodspec.extract_arc_spec(display=False,
                                      filename=os.path.join(
                                          reduced_data_folder_path, obsnight,
                                          light_name + '_arc_spec'),
                                      save_iframe=True)

            onedspec = spectral_reduction.OneDSpec()
            onedspec.from_twodspec(twodspec, stype='science')

            # Find the peaks of the arc
            onedspec.find_arc_lines(display=False,
                                    stype='science',
                                    filename=os.path.join(
                                        reduced_data_folder_path, obsnight,
                                        light_name + '_arc_lines'),
                                    save_iframe=True)

            # Configure the wavelength calibrator
            onedspec.initialise_calibrator(stype='science+')
            onedspec.set_hough_properties(num_slopes=500,
                                          xbins=100,
                                          ybins=100,
                                          min_wavelength=3500,
                                          max_wavelength=8000,
                                          stype='science')

            onedspec.set_ransac_properties(filter_close=True, stype='science')

            onedspec.load_user_atlas(elements=element,
                                     wavelengths=atlas,
                                     stype='science')
            onedspec.do_hough_transform()

            # Solve for the pixel-to-wavelength solution
            onedspec.fit(max_tries=500, stype='science', display=False)

            # Apply the wavelength calibration and display it
            onedspec.apply_wavelength_calibration(stype='science')

            # Get the sensitivity curve from the pre-calibrated data
            # Dividing HDU #5 by #4
            calibrated_data = fits.open(light_file_path)
            sensitivity = calibrated_data[4].data[0] / calibrated_data[3].data[
                0]

            wave_bin = float(calibrated_data[4].header['CDELT1'])
            wave_start = float(
                calibrated_data[4].header['CRVAL1']) + wave_bin / 2.
            wave_end = wave_start + (int(calibrated_data[3].header['NAXIS1']) -
                                     1) * wave_bin

            wave = np.linspace(wave_start, wave_end,
                               int(calibrated_data[3].header['NAXIS1']))

            # interpolate the senstivity curve with wavelength
            sensitivity_itp = itp.interp1d(wave,
                                           np.log10(sensitivity),
                                           fill_value='extrapolate')

            onedspec.add_sensitivity_func(sensitivity_itp)

            onedspec.apply_flux_calibration(stype='science')

            # Inspect reduced spectrum
            onedspec.inspect_reduced_spectrum(
                stype='science',
                display=False,
                filename=os.path.join(reduced_data_folder_path, obsnight,
                                      light_name + '_reduced'),
                save_png=True,
                save_pdf=True,
                save_iframe=True)

            # Save as a FITS file
            onedspec.save_fits(output='flux_resampled',
                               filename=os.path.join(reduced_data_folder_path,
                                                     obsnight,
                                                     light_name + '_reduced'),
                               stype='science',
                               overwrite=True)

            # [1] Update the DB
            update_query = 'UPDATE "{}" SET "sensitivity" = "{}", "output" = \
                "{}" WHERE "filename" = "{}";'.format(
                reduction_table, row[0],
                os.path.join(reduced_data_folder_path, obsnight,
                             light_name + '_reduced'), row[0])
            logging.info(update_query)

            # [2] Post message to Slack
            post_msg_res = post_message_to_slack(
                webhook_url,
                "{} is reduced and is now available at {}.".format(
                    row[0], os.path.join(reduced_data_folder_path, obsnight)))
            logging.info(post_msg_res)

            # [3 + 4] Post images to Slack
            post_file_res = post_file_to_slack(
                slack_token,
                os.path.join(reduced_data_folder_path, obsnight,
                             light_name + '_reduced_0.png'))
            post_file_res = post_file_to_slack(
                slack_token,
                os.path.join(reduced_data_folder_path, obsnight,
                             light_name + '_reduced_0.html'))
            logging.info(post_file_res)

            # Commit changes if all 4 actions are successful
            db_connector.execute(update_query)
            db_connector.commit()
            logging.info('Commited to the DB.')
