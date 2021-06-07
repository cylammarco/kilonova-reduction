# kilonova-reduction

This is an automated quicklook data grabber and reducer. Because SPRAT data is already flattened and ready for reduction, the Primary HDU can be used directly for extraction. An arc is always taking at the end of an observing group (observing block in a more common operation term). It checks the quicklook every 5 minutes between two set hours (default: 1900 to 0700) and a reduction script is run every 5 minutes and if it finds a new arc being downloaded, it will reduce all the data from that group and attempt to stack after extraction (weighted by the uncertainties).

The reduced plots are sent to slack alongside with the paths to the reduced data.

This is designed for immediate data reduction, so if it was stuck/not running overnight, when it is restarted, it does not attempt to check if there is missing data.

## How to use it

1. Modify the config.json accordingly.

2. Run the script grab_lt_data.py headless, it will prompt for a password.

3. Run the script reduce_lt_data.py headless, it will prompt for a password.

4. Wait for data to be collected

log files will be generated and stored in where the `config.json` sets them to be, with the keywords `data_log_file` and `reduction_log_file`. See the explanation below:

### What do the keywords in config.json mean?

* `base_folder` - The path to where all the logs, raw and reduced data are going. If a relative path is given, it will go to `~/base_folder`.
* `data_log_file` - The name of the log file for data grabbing from the LT quicklook.
* `data_log_level` - The log level - 'INFO', 'DEBUG' (There are 5 levels, but only these two are used)
* `reduction_log_file` - The name of the log file for data reduction.
* `reduction_log_level` - The log level - 'INFO', 'DEBUG' (There are 5 levels, but only these two are used)
* `propid` - List of Proposal IDs.
* `db_name` - The name of the SQLite database for holding the `obslog_table`, `data_table` and `reduction_table`.
* `obslog_table` - Basically the obslog from LT quicklook, except the time UTC included the date (the LT one only includes the time).
* `data_table` - The table listing the downloaded files.
* `reduction_table` - The table listing the reduction status of the files.
* `raw_data_folder` - Folder in the `base_folder` to put the raw data.
* `reduced_data_folder` - Folder in the `base_folder` to put the reduced data.
* `time_start` - The hours to start the data grabbing.
* `time_end` - The hours to end the data grabbing.
* `login_config` - List of files containing encrypted login details for downloading the raw data.
* `slack_config` - The slack webhool url and token.
