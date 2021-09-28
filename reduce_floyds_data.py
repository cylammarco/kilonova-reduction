import copy
import sys

from astropy.io import fits
from astroscrappy import detect_cosmics
from aspired import spectral_reduction
from matplotlib import pyplot as plt
from spectres import spectres
from statsmodels.nonparametric.smoothers_lowess import lowess
import numpy as np


def flux_diff(ratio, a, b):
    diff = a * ratio - b
    mask = diff < np.nanpercentile(diff, 95)
    return np.nansum(diff[mask]**2.)


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
    blue.apply_mask_to_arc()
    red.add_arc(arc_fits.data, fits.Header(arc_fits.header))
    red.apply_mask_to_arc()

    # Get the trace to rectify the image
    red.ap_trace(nspec=1,
                 ap_faint=20,
                 percentile=10,
                 trace_width=50,
                 shift_tol=35,
                 fit_deg=5,
                 display=False)
    red.get_rectification(upsample_factor=5, coeff=coeff_red, display=True)
    red.apply_rectification()
    # Need to store the traces for fringe correction before overwriting them
    # with the new traces
    trace_red = copy.deepcopy(red.spectrum_list[0].trace)
    trace_sigma_red = copy.deepcopy(red.spectrum_list[0].trace_sigma)

    # Get the trace again for the rectified image and then extract
    red.ap_trace(nspec=1,
                 ap_faint=20,
                 percentile=10,
                 trace_width=20,
                 fit_deg=3,
                 display=False)
    red.ap_extract(apwidth=10, spec_id=0, display=False)

    # Do the same with the blue
    blue.ap_trace(nspec=1,
                  ap_faint=10,
                  percentile=10,
                  trace_width=20,
                  shift_tol=50,
                  fit_deg=5,
                  display=False)
    blue.get_rectification(upsample_factor=5, coeff=coeff_blue, display=True)
    blue.apply_rectification()

    blue.ap_trace(nspec=1,
                  percentile=10,
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
    flat.get_rectification(coeff=red.rec_coeff)
    flat.apply_rectification()
    flat.add_trace(trace_red_rectified, trace_sigma_red_rectified)
    flat.ap_extract(apwidth=10, skywidth=0, display=False)

    return red, blue, flat


# only red specs are used
def fringe_correction(target, flat):
    flat_count = flat.spectrum_list[0].count
    flat_continuum = lowess(flat_count,
                            np.arange(len(flat_count)),
                            frac=0.04,
                            return_sorted=False)
    flat_continuum_divided = flat_count / flat_continuum
    flat_continuum_divided[:500] = 1.0
    # Apply the flat correction
    target.spectrum_list[0].count /= flat_continuum_divided


def calibrate_red(science, standard, standard_name):
    #
    # Start handling 1D spectra here
    #
    # Need to add fringe subtraction here
    red_onedspec = spectral_reduction.OneDSpec(log_level='INFO',
                                               log_file_name=None)

    # Red spectrum first
    red_onedspec.from_twodspec(standard, stype='standard')
    red_onedspec.from_twodspec(science, stype='science')

    # Find the peaks of the arc
    red_onedspec.find_arc_lines(display=True,
                                prominence=0.25,
                                top_n_peaks=25,
                                stype='science')
    red_onedspec.find_arc_lines(display=False,
                                prominence=0.25,
                                top_n_peaks=25,
                                stype='standard')

    # Configure the wavelength calibrator
    red_onedspec.initialise_calibrator(stype='science+standard')

    red_onedspec.add_user_atlas(elements=element_Hg_red,
                                wavelengths=atlas_Hg_red,
                                stype='science+standard')
    red_onedspec.add_user_atlas(elements=element_Ar_red,
                                wavelengths=atlas_Ar_red,
                                stype='science+standard')

    red_onedspec.set_hough_properties(num_slopes=5000,
                                      xbins=200,
                                      ybins=200,
                                      min_wavelength=4750,
                                      max_wavelength=10750,
                                      stype='science+standard')
    red_onedspec.set_ransac_properties(minimum_matches=13,
                                       stype='science+standard')
    red_onedspec.do_hough_transform(stype='science+standard')

    # Solve for the pixel-to-wavelength solution
    red_onedspec.fit(max_tries=5000, stype='science+standard', display=True)

    # Apply the wavelength calibration and display it
    red_onedspec.apply_wavelength_calibration(wave_start=5000,
                                              wave_end=11000,
                                              wave_bin=1,
                                              stype='science+standard')

    red_onedspec.load_standard(standard_name)
    red_onedspec.get_sensitivity()
    red_onedspec.apply_flux_calibration()
    red_onedspec.get_telluric_profile()
    red_onedspec.apply_telluric_correction()

    return red_onedspec


def calibrate_blue(science, standard, standard_name):
    # Blue spectrum here
    blue = spectral_reduction.OneDSpec(log_level='INFO', log_file_name=None)

    blue.from_twodspec(standard, stype='standard')
    blue.from_twodspec(science, stype='science')

    blue.find_arc_lines(prominence=0.25,
                        top_n_peaks=10,
                        display=True,
                        stype='science')
    blue.find_arc_lines(prominence=0.25,
                        top_n_peaks=10,
                        display=False,
                        stype='standard')

    blue.initialise_calibrator(stype='science+standard')

    blue.add_user_atlas(elements=element_Hg_blue,
                        wavelengths=atlas_Hg_blue,
                        stype='science+standard')
    blue.add_user_atlas(elements=element_Ar_blue,
                        wavelengths=atlas_Ar_blue,
                        stype='science+standard')
    blue.add_user_atlas(elements=element_Zn_blue,
                        wavelengths=atlas_Zn_blue,
                        stype='science+standard')

    blue.set_hough_properties(num_slopes=2000,
                              xbins=200,
                              ybins=200,
                              min_wavelength=3250,
                              max_wavelength=5850,
                              stype='science+standard')
    blue.set_ransac_properties(filter_close=True, stype='science+standard')
    blue.do_hough_transform(stype='science+standard')

    # Solve for the pixel-to-wavelength solution
    blue.fit(max_tries=5000, stype='science+standard', display=True)

    # Apply the wavelength calibration and display it
    blue.apply_wavelength_calibration(wave_start=3200,
                                      wave_end=5600,
                                      wave_bin=1,
                                      stype='science+standard')

    blue.load_standard(standard_name)
    blue.get_sensitivity()
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

atlas_Ar_blue = [4158.590, 4200.674, 4272.169, 4300.101, 4510.733]
atlas_Hg_blue = [
    3650.153, 4046.563, 4077.8314, 4358.328, 5460.7348, 5769.598, 5790.663
]
atlas_Zn_blue = [4722.1569]
element_Hg_blue = ['Hg'] * len(atlas_Hg_blue)
element_Ar_blue = ['Ar'] * len(atlas_Ar_blue)
element_Zn_blue = ['Zn'] * len(atlas_Zn_blue)

# Set the frame
red_spatial_mask = np.arange(0, 330)
blue_spatial_mask = np.arange(335, 512)
red_spec_mask = np.arange(0, 1800)
blue_spec_mask = np.arange(500, 2060)

filelist = np.genfromtxt(sys.argv[1], delimiter=',', dtype='U', autostrip=True)

obstype = np.array(filelist[:, 0]).astype(str)
frametype = np.array(filelist[:, 1]).astype(str)
filename = np.array(filelist[:, 2]).astype(str)

#
# Standard fits here (NOT just an array, need the )
#
standard_light_fits = fits.open(filename[(obstype == 'standard')
                                         & (frametype == 'light')][0])[1]
standard_flat_fits = fits.open(filename[(obstype == 'standard')
                                        & (frametype == 'flat')][0])[1]
standard_arc_fits = fits.open(filename[(obstype == 'standard')
                                       & (frametype == 'arc')][0])[1]
standard_name = standard_light_fits.header['OBJECT']

#
# Science frame here
#
science_light_fits = fits.open(filename[(obstype == 'science')
                                        & (frametype == 'light')][0])[1]
science_flat_fits = fits.open(filename[(obstype == 'science')
                                       & (frametype == 'flat')][0])[1]
science_arc_fits = fits.open(filename[(obstype == 'science')
                                      & (frametype == 'arc')][0])[1]
science_name = science_light_fits.header['OBJECT']

science_light_fits.data = detect_cosmics(
    np.array(science_light_fits.data).astype('float') / 2.3,
    gain=2.3,
    readnoise=3.5,
    sigclip=3.5,
    fsmode='convolve',
    psfmodel='gaussy',
    psfsize=11)[1]

standard_light_fits.data = detect_cosmics(
    np.array(standard_light_fits.data).astype('float') / 2.3,
    gain=2.3,
    readnoise=3.5,
    sigclip=3.5,
    fsmode='convolve',
    psfmodel='gaussy',
    psfsize=11)[1]

rec_coeff_blue = None
rec_coeff_red = None
# rec_coeff_blue = (-1.28091530e+02, 1.11427370e-01, 5.19032671e-07)
# rec_coeff_red = (-1.52689504e+02, 1.31426850e-01, -1.25352165e-06)

standard_red, standard_blue, standard_flat = extract_floyds(
    standard_light_fits,
    standard_flat_fits,
    standard_arc_fits,
    coeff_red=rec_coeff_red,
    coeff_blue=rec_coeff_blue)
science_red, science_blue, science_flat = extract_floyds(
    science_light_fits,
    science_flat_fits,
    science_arc_fits,
    coeff_red=rec_coeff_red,
    coeff_blue=rec_coeff_blue)

fringe_correction(science_red, science_flat)
fringe_correction(standard_red, standard_flat)

onedspec_blue = calibrate_blue(science_blue, standard_blue, standard_name)
onedspec_red = calibrate_red(science_red, standard_red, standard_name)

onedspec_red.get_telluric_profile()
onedspec_red.apply_telluric_correction()

# Inspect
onedspec_blue.inspect_reduced_spectrum(
    wave_min=3000,
    wave_max=6000,
    stype='science',
    display=True,
    save_fig=True,
    filename=filename[(obstype == 'science')
                      & (frametype == 'light')][0].split('.')[0] + '_blue')

onedspec_red.inspect_reduced_spectrum(
    wave_min=5000,
    wave_max=10500,
    stype='science',
    display=True,
    save_fig=True,
    filename=filename[(obstype == 'science')
                      & (frametype == 'light')][0].split('.')[0] + '_red')

onedspec_blue.inspect_reduced_spectrum(
    wave_min=3000,
    wave_max=6000,
    stype='standard',
    display=True,
    save_fig=True,
    filename=filename[(obstype == 'standard')
                      & (frametype == 'light')][0].split('.')[0] + '_blue')

onedspec_red.inspect_reduced_spectrum(
    wave_min=5000,
    wave_max=10500,
    stype='standard',
    display=True,
    save_fig=True,
    filename=filename[(obstype == 'standard')
                      & (frametype == 'light')][0].split('.')[0] + '_red')

wave_red = onedspec_red.science_spectrum_list[0].wave_resampled
wave_blue = onedspec_blue.science_spectrum_list[0].wave_resampled

flux_red = onedspec_red.science_spectrum_list[0].flux_resampled
flux_blue = onedspec_blue.science_spectrum_list[0].flux_resampled

flux_red_err = onedspec_red.science_spectrum_list[0].flux_err_resampled
flux_blue_err = onedspec_blue.science_spectrum_list[0].flux_err_resampled

# trim the last ~100A from the blue and the first ~100A from the red
# in the combined spectrum
red_limit = 5000
blue_limit = 5800

blue_mask = (wave_blue >= red_limit) & (wave_blue <= blue_limit)
red_mask = (wave_red >= red_limit) & (wave_red <= blue_limit)

# resample the red to match blue resolution
flux_red_resampled, flux_red_resampled_err = spectres(wave_blue[blue_mask],
                                                      wave_red[red_mask],
                                                      flux_red[red_mask],
                                                      flux_red_err[red_mask])

flux_weighted_combine = (flux_red_resampled / flux_red_resampled_err +
                         flux_blue[blue_mask] / flux_blue_err[blue_mask]) / (
                             1 / flux_red_resampled_err +
                             1 / flux_blue_err[blue_mask])

plt.figure(1, figsize=(16, 8))
plt.clf()
plt.plot(wave_blue, flux_blue, color='blue', label='Blue arm (ASPIRED)')
plt.plot(wave_red, flux_red, color='red', label='Red arm (ASPIRED)')
plt.xlim(3300., 10500.)
plt.ylim(
    0,
    max(np.nanpercentile(flux_red_resampled, 99.5),
        np.nanpercentile(flux_blue, 99.5)))
plt.xlabel('Wavelength / A')
plt.ylabel('Flux / (erg / s / cm / cm / A)')
plt.legend()
plt.grid()
plt.tight_layout()
plt.title(science_name)
plt.show()
