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
from os.path import splitext
import sys

def extract_floyds(light_fits,
                   flat_fits,
                   arc_fits,
                   coeff_red=None,
                   coeff_blue=None):
    light_data = light_fits.data
    light_header = fits.Header(light_fits.header)

    red = spectral_reduction.TwoDSpec(light_data,
                                      header=light_header,
                                      spatial_mask=red_spatial_mask,
                                      spec_mask=red_spec_mask,
                                      cosmicray=True,
                                      sigclip=2.,
                                      readnoise=3.5,
                                      gain=2.3,
                                      log_level='INFO',
                                      log_file_name=None)

    blue = spectral_reduction.TwoDSpec(light_data,
                                       header=light_header,
                                       spatial_mask=blue_spatial_mask,
                                       spec_mask=blue_spec_mask,
                                       cosmicray=True,
                                       sigclip=2.,
                                       readnoise=3.5,
                                       gain=2.3,
                                       log_level='INFO',
                                       log_file_name=None)

    # Add the arcs before rectifying the image, which will apply the
    # rectification to the arc frames too
    blue.add_arc(arc_fits.data, fits.Header(arc_fits.header))
    blue.apply_twodspec_mask_to_arc()
    red.add_arc(arc_fits.data, fits.Header(arc_fits.header))
    red.apply_twodspec_mask_to_arc()

    # Get the trace to rectify the image
    red.ap_trace(nspec=1,
                 ap_faint=20,
                 trace_width=20,
                 shift_tol=50,
                 fit_deg=5,
                 display=False)
    red.compute_rectification(upsample_factor=10, coeff=coeff_red)
    red.apply_rectification()
    # Need to store the traces for fringe correction before overwriting them
    # with the new traces
    trace_red = red.spectrum_list[0].trace
    trace_sigma_red = red.spectrum_list[0].trace_sigma

    # Get the trace again for the rectified image and then extract
    red.ap_trace(nspec=1, trace_width=20, fit_deg=3, display=False)
    red.ap_extract(apwidth=10, spec_id=0, display=False)

    # Do the same with the blue
    blue.ap_trace(nspec=1,
                  ap_faint=20,
                  trace_width=20,
                  shift_tol=50,
                  fit_deg=5,
                  display=False)
    blue.compute_rectification(upsample_factor=10,
                               bin_size=7,
                               n_bin=[3, 1],
                               coeff=coeff_blue)
    blue.apply_rectification()

    blue.ap_trace(nspec=1,
                  percentile=30,
                  trace_width=20,
                  fit_deg=3,
                  display=False)
    blue.ap_extract(apwidth=10, display=False)

    red.extract_arc_spec(spec_width=15, display=False)
    blue.extract_arc_spec(spec_width=15, display=False)

    # Get the blue and red traces for fringe removal
    trace_red_rectified = red.spectrum_list[0].trace
    trace_sigma_red_rectified = red.spectrum_list[0].trace_sigma

    # Extract the red flat
    flat_data = flat_fits.data
    flat_header = fits.Header(flat_fits.header)

    flat = spectral_reduction.TwoDSpec(flat_data,
                                       header=flat_header,
                                       spatial_mask=red_spatial_mask,
                                       spec_mask=red_spec_mask,
                                       cosmicray=True,
                                       readnoise=3.5,
                                       gain=2.3,
                                       log_level='INFO',
                                       log_file_name=None)

    # Force extraction from the flat for fringe correction
    flat.add_trace(trace_red, trace_sigma_red)
    flat.compute_rectification(coeff=coeff_red)
    flat.apply_rectification()
    flat.add_trace(trace_red_rectified, trace_sigma_red_rectified)
    flat.ap_extract(apwidth=10, skywidth=0, display=False)

    return red, blue, flat


def calibrate_red(science, standard, standard_name):
    #
    # Start handling 1D spectra here
    #
    # Need to add fringe subtraction here
    red = spectral_reduction.OneDSpec(log_level='INFO', log_file_name=None)

    # Red spectrum first
    red.from_twodspec(standard, stype='standard')
    red.from_twodspec(science, stype='science')

    # Find the peaks of the arc
    red.find_arc_lines(display=False, prominence=100, stype='science')
    red.find_arc_lines(display=False, prominence=20, stype='standard')

    # Configure the wavelength calibrator
    red.initialise_calibrator(stype='science+standard')

    red.add_user_atlas(elements=element_Hg_red,
                       wavelengths=atlas_Hg_red,
                       stype='science+standard')
    red.add_user_atlas(elements=element_Ar_red,
                       wavelengths=atlas_Ar_red,
                       stype='science+standard')

    red.set_hough_properties(num_slopes=2000,
                             xbins=100,
                             ybins=100,
                             min_wavelength=4500,
                             max_wavelength=10500,
                             stype='science+standard')
    red.set_ransac_properties(stype='science+standard')
    red.do_hough_transform(stype='science+standard')

    # Solve for the pixel-to-wavelength solution
    red.fit(max_tries=2000, stype='science+standard', display=False, savefig=True)

    # Apply the wavelength calibration and display it
    red.apply_wavelength_calibration(wave_start=5000,
                                     wave_end=11000,
                                     wave_bin=1,
                                     stype='science+standard')

    red.load_standard(standard_name)
    red.compute_sensitivity()
    red.apply_flux_calibration()
    '''

    L745_fringe_count = L745_twodspec_flat.spectrum_list[0].count

    L745_fringe_continuum = lowess(L745_fringe_count,
                                np.arange(len(L745_fringe_count)),
                                frac=0.04,
                                return_sorted=False)
    L745_fringe_normalised = L745_fringe_count / L745_fringe_continuum

    L745_red_count = L745_twodspec_red.spectrum_list[0].count
    L745_red_continuum = lowess(L745_red_count,
                                np.arange(len(L745_red_count)),
                                frac=0.04,
                                return_sorted=False)
    L745_red_normalised = L745_red_count / L745_red_continuum
    L745_sed_correction = L745_fringe_continuum / L745_red_continuum
    L745_sed_correction /= np.nanmean(L745_sed_correction)

    L745_factor = (np.nanpercentile(L745_fringe_normalised[1000:1800], 95)) / (
        np.nanpercentile(L745_red_normalised[1000:1800], 5))
    L745_factor_mean = np.nanmean(L745_factor)

    L745_fringe_correction =\
        L745_fringe_normalised / L745_factor_mean *\
            L745_red_continuum * L745_sed_correction

    # Apply the flat correction
    L745_twodspec_red.spectrum_list[0].count -= L745_fringe_correction

    '''

    return red


def calibrate_blue(science, standard, standard_name):
    # Blue spectrum here
    blue = spectral_reduction.OneDSpec(log_level='INFO', log_file_name=None)

    blue.from_twodspec(standard, stype='standard')
    blue.from_twodspec(science, stype='science')

    blue.find_arc_lines(prominence=10, display=False, stype='science')
    blue.find_arc_lines(prominence=5, display=False, stype='standard')

    blue.initialise_calibrator(stype='science+standard')

    blue.add_user_atlas(elements=element_Hg_blue,
                        wavelengths=atlas_Hg_blue,
                        stype='science+standard')
    blue.add_user_atlas(elements=element_Zn_blue,
                        wavelengths=atlas_Zn_blue,
                        stype='science+standard')

    blue.set_hough_properties(num_slopes=2000,
                              xbins=200,
                              ybins=200,
                              min_wavelength=3000,
                              max_wavelength=6000,
                              stype='science+standard')
    blue.set_ransac_properties(filter_close=True, stype='science+standard')
    blue.do_hough_transform(stype='science+standard')

    # Solve for the pixel-to-wavelength solution
    blue.fit(max_tries=1000, stype='science+standard', display=False, savefig=True)

    # Apply the wavelength calibration and display it
    blue.apply_wavelength_calibration(wave_start=3000,
                                      wave_end=5800,
                                      wave_bin=1,
                                      stype='science+standard')

    blue.load_standard(standard_name)
    blue.compute_sensitivity()
    blue.apply_flux_calibration()

    return blue


# Line list
atlas_Hg_red = [5460.7348, 5769.5982, 5790.6630]
atlas_Ar_red = [
    6965.4307, 7067.2175, 7147.0416, 7272.9359, 7383.9805, 7503.8691,
    7635.1056, 7723.7599, 7948.1764, 8014.7857, 8115.3108, 8264.5225,
    8424.6475, 8521.4422, 9122.9674, 9224.4992, 9657.7863
]
element_Hg_red = ['Hg'] * len(atlas_Hg_red)
element_Ar_red = ['Ar'] * len(atlas_Ar_red)

atlas_Hg_blue = [
    3650.153, 4046.563, 4077.8314, 4358.328, 4916.068, 5460.7348, 5769.5982
]
atlas_Zn_blue = [4078.14, 4298.3249, 4722.15, 4810.53, 5181.9819]
element_Hg_blue = ['Hg'] * len(atlas_Hg_blue)
element_Zn_blue = ['Zn'] * len(atlas_Zn_blue)

# Set the frame
red_spatial_mask = np.arange(0, 330)
blue_spatial_mask = np.arange(335, 512)
red_spec_mask = np.arange(0, 1900)
blue_spec_mask = np.arange(500, 2060)

filelist = np.genfromtxt(sys.argv[1], delimiter=',', dtype='U', autostrip=True)

obstype = np.array(filelist[:, 0]).astype(str)
frametype = np.array(filelist[:, 1]).astype(str)
filename = np.array(filelist[:, 2]).astype(str)

#
# Standard frame here
#
standard_light_fits = fits.open(filename[(obstype == 'standard')
                                         & (frametype == 'light')][0])[1]
standard_flat_fits = fits.open(filename[(obstype == 'standard')
                                        & (frametype == 'flat')][0])[1]
standard_arc_fits = fits.open(filename[(obstype == 'standard')
                                       & (frametype == 'arc')][0])[1]
standard_name = standard_light_fits.header['OBJECT']

standard_red, standard_blue, standard_flat =\
    extract_floyds(standard_light_fits, standard_flat_fits, standard_arc_fits)

#
# Science frame here
#
science_light_fits = fits.open(filename[(obstype == 'science')
                                        & (frametype == 'light')][0])[1]
science_flat_fits = fits.open(filename[(obstype == 'science')
                                       & (frametype == 'flat')][0])[1]
science_arc_fits = fits.open(filename[(obstype == 'science')
                                      & (frametype == 'arc')][0])[1]

science_red, science_blue, science_flat =\
    extract_floyds(science_light_fits, science_flat_fits, science_arc_fits,
        coeff_red=standard_red.rec_coeff,
        coeff_blue=standard_blue.rec_coeff)

onedspec_blue = calibrate_blue(science_blue, standard_blue, standard_name)
onedspec_red = calibrate_red(science_red, standard_red, standard_name)

# Inspect
onedspec_blue.inspect_reduced_spectrum(
    wave_min=3000,
    wave_max=6000,
    stype='science',
    display=False,
    save_iframe=True,
    filename=splitext(filename[(obstype == 'science')
                      & (frametype == 'light')][0])[0] + '_blue')

onedspec_red.inspect_reduced_spectrum(
    wave_min=5000,
    wave_max=11000,
    stype='science',
    display=False,
    save_iframe=True,
    filename=splitext(filename[(obstype == 'science')
                      & (frametype == 'light')][0])[0] + '_red')

onedspec_blue.inspect_reduced_spectrum(
    wave_min=3000,
    wave_max=6000,
    stype='standard',
    display=False,
    save_iframe=True,
    filename=splitext(filename[(obstype == 'standard')
                      & (frametype == 'light')][0])[0] + '_blue')

onedspec_red.inspect_reduced_spectrum(
    wave_min=5000,
    wave_max=11000,
    stype='standard',
    display=False,
    save_iframe=True,
    filename=splitext(filename[(obstype == 'standard')
                      & (frametype == 'light')][0])[0] + '_red')
