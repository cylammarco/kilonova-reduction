from astropy.io import fits
from aspired import spectral_reduction
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from spectres import spectres
from statsmodels.nonparametric.smoothers_lowess import lowess
import numpy as np
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
                                      sigclip=3.,
                                      readnoise=3.5,
                                      gain=2.3,
                                      log_level='WARNING',
                                      log_file_name=None)

    blue = spectral_reduction.TwoDSpec(light_data,
                                       header=light_header,
                                       spatial_mask=blue_spatial_mask,
                                       spec_mask=blue_spec_mask,
                                       cosmicray=True,
                                       sigclip=3.,
                                       readnoise=3.5,
                                       gain=2.3,
                                       log_level='WARNING',
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
    red.compute_rectification(upsample_factor=10,
                              coeff=coeff_red,
                              display=True,
                              save_iframe=True)
    red.apply_rectification()
    # Need to store the traces for fringe correction before overwriting them
    # with the new traces
    trace_red = red.spectrum_list[0].trace
    trace_sigma_red = red.spectrum_list[0].trace_sigma

    # Get the trace again for the rectified image and then extract
    red.ap_trace(nspec=1, trace_width=20, fit_deg=3, display=False)
    red.ap_extract(apwidth=15, spec_id=0, display=True)

    # Do the same with the blue
    blue.ap_trace(nspec=1,
                  ap_faint=20,
                  trace_width=20,
                  shift_tol=50,
                  fit_deg=5,
                  display=False)
    blue.compute_rectification(upsample_factor=10,
                               n_bin=[5, 2],
                               coeff=coeff_blue,
                               display=True,
                               save_iframe=True)
    blue.apply_rectification()

    blue.ap_trace(nspec=1,
                  percentile=30,
                  trace_width=20,
                  fit_deg=3,
                  display=False)
    blue.ap_extract(apwidth=15, display=True)

    red.extract_arc_spec(spec_width=30, display=False)
    blue.extract_arc_spec(spec_width=30, display=False)

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
                                       log_level='WARNING',
                                       log_file_name=None)

    # Force extraction from the flat for fringe correction
    flat.add_trace(trace_red, trace_sigma_red)
    # This line will throw a warning about "Arc frame not available"
    flat.compute_rectification(coeff=coeff_red)
    flat.apply_rectification()
    flat.add_trace(trace_red_rectified, trace_sigma_red_rectified)
    flat.ap_extract(apwidth=10, skywidth=0, display=False)

    return red, blue, flat


def match_fringe_amplitude(coeff, red_continuum_subtracted, fringe_normalised):
    factor = np.polyval(coeff, np.arange(850, 1600))
    diff = red_continuum_subtracted[
        850:1600] - fringe_normalised[850:1600] * factor
    return np.nansum(diff**2.)


def calibrate_red(science,
                  standard,
                  standard_name,
                  science_flat=None,
                  standard_flat=None):
    #
    # Start handling 1D spectra here
    #
    red_onedspec = spectral_reduction.OneDSpec(log_level='WARNING',
                                               log_file_name=None)

    # Finge subtraction
    if science_flat is not None:

        science_fringe_count = science_flat.spectrum_list[0].count

        science_fringe_continuum = lowess(science_fringe_count,
                                          np.arange(len(science_fringe_count)),
                                          frac=0.02,
                                          return_sorted=False)
        science_fringe_normalised = science_fringe_count - science_fringe_continuum

        science_red_count = science_red.spectrum_list[0].count
        science_red_continuum = lowess(science_red_count,
                                       np.arange(len(science_red_count)),
                                       frac=0.02,
                                       return_sorted=False)
        science_red_continuum_subtracted = science_red_count - science_red_continuum

        science_coeff = minimize(match_fringe_amplitude, (1, 1),
                                 args=(science_red_continuum_subtracted,
                                       science_fringe_normalised),
                                 method='Nelder-Mead',
                                 options={
                                     'maxiter': 10000
                                 }).x
        science_factor = np.polyval(science_coeff,
                                    np.arange(len(science_fringe_normalised)))
        science_fringe_correction = science_fringe_normalised * science_factor
        science_fringe_correction[:450] = 0.

        # Apply the flat correction
        science.spectrum_list[0].count -= science_fringe_correction

    if (science_flat is not None) & (standard_flat is None):

        standard_flat = science_flat

    if standard_flat is not None:

        standard_fringe_count = standard_flat.spectrum_list[0].count

        standard_fringe_continuum = lowess(standard_fringe_count,
                                           np.arange(
                                               len(standard_fringe_count)),
                                           frac=0.02,
                                           return_sorted=False)
        standard_fringe_normalised = standard_fringe_count - standard_fringe_continuum

        standard_red_count = standard_red.spectrum_list[0].count
        standard_red_continuum = lowess(standard_red_count,
                                        np.arange(len(standard_red_count)),
                                        frac=0.02,
                                        return_sorted=False)
        standard_red_continuum_subtracted = standard_red_count - standard_red_continuum

        standard_coeff = minimize(match_fringe_amplitude, (1, 1),
                                  args=(standard_red_continuum_subtracted,
                                        standard_fringe_normalised),
                                  method='Nelder-Mead',
                                  options={
                                      'maxiter': 10000
                                  }).x
        standard_factor = np.polyval(
            standard_coeff, np.arange(len(standard_fringe_normalised)))
        standard_fringe_correction = standard_fringe_normalised * standard_factor
        standard_fringe_correction[:450] = 0.

        # Apply the flat correction
        science.spectrum_list[0].count -= standard_fringe_correction

    # Red spectrum first
    red_onedspec.from_twodspec(standard, stype='standard')
    red_onedspec.from_twodspec(science, stype='science')

    # Find the peaks of the arc
    red_onedspec.find_arc_lines(display=False,
                                prominence=0.25,
                                top_n_peaks=25,
                                stype='science')
    red_onedspec.find_arc_lines(display=False,
                                prominence=0.25,
                                top_n_peaks=25,
                                stype='standard')

    # Configure the wavelength calibrator
    red_onedspec.initialise_calibrator(stype='science+standard')

    red_onedspec.add_user_atlas(
        elements=element_Hg_red,
        wavelengths=atlas_Hg_red,
        pressure=float(science.header['REFPRES']) * 100,
        temperature=float(science.header['REFTEMP'] + 273),
        relative_humidity=float(science.header['REFHUMID']),
        stype='science+standard')
    red_onedspec.add_user_atlas(
        elements=element_Ar_red,
        wavelengths=atlas_Ar_red,
        pressure=float(standard.header['REFPRES']) * 100,
        temperature=float(standard.header['REFTEMP'] + 273),
        relative_humidity=float(standard.header['REFHUMID']),
        stype='science+standard')

    red_onedspec.set_hough_properties(num_slopes=2000,
                                      xbins=100,
                                      ybins=100,
                                      min_wavelength=4500,
                                      max_wavelength=10500,
                                      stype='science+standard')
    red_onedspec.set_ransac_properties(stype='science+standard')
    red_onedspec.do_hough_transform(stype='science+standard')

    # Solve for the pixel-to-wavelength solution
    red_onedspec.fit(max_tries=1000,
                     stype='science+standard',
                     display=True,
                     savefig=True)

    # Apply the wavelength calibration and display it
    red_onedspec.apply_wavelength_calibration(wave_start=4750,
                                              wave_end=11000,
                                              wave_bin=2,
                                              stype='science+standard')

    red_onedspec.load_standard(standard_name)
    red_onedspec.compute_sensitivity(frac=0.25)
    red_onedspec.inspect_sensitivity()
    red_onedspec.apply_flux_calibration()

    return red_onedspec


def calibrate_blue(science, standard, standard_name):
    # Blue spectrum here
    blue_onedspec = spectral_reduction.OneDSpec(log_level='WARNING',
                                                log_file_name=None)

    blue_onedspec.from_twodspec(standard, stype='standard')
    blue_onedspec.from_twodspec(science, stype='science')

    blue_onedspec.find_arc_lines(prominence=0.25,
                                 top_n_peaks=10,
                                 display=True,
                                 stype='science')
    blue_onedspec.find_arc_lines(prominence=0.25,
                                 top_n_peaks=10,
                                 display=True,
                                 stype='standard')

    blue_onedspec.initialise_calibrator(stype='science+standard')

    blue_onedspec.add_user_atlas(
        elements=element_Hg_blue,
        wavelengths=atlas_Hg_blue,
        pressure=float(science.header['REFPRES']) * 100,
        temperature=float(science.header['REFTEMP'] + 273),
        relative_humidity=float(science.header['REFHUMID']),
        stype='science+standard')
    blue_onedspec.add_user_atlas(
        elements=element_Ar_blue,
        wavelengths=atlas_Ar_blue,
        pressure=float(standard.header['REFPRES']) * 100,
        temperature=float(standard.header['REFTEMP'] + 273),
        relative_humidity=float(standard.header['REFHUMID']),
        stype='science+standard')

    blue_onedspec.set_hough_properties(num_slopes=2000,
                                       xbins=100,
                                       ybins=100,
                                       min_wavelength=3000,
                                       max_wavelength=6000,
                                       stype='science+standard')
    blue_onedspec.set_ransac_properties(filter_close=True,
                                        sample_size=3,
                                        top_n_candidate=4,
                                        stype='science+standard')
    blue_onedspec.do_hough_transform(stype='science+standard')

    # Solve for the pixel-to-wavelength solution
    blue_onedspec.fit(max_tries=1000,
                      stype='science+standard',
                      display=True,
                      savefig=True)

    # Apply the wavelength calibration and display it
    blue_onedspec.apply_wavelength_calibration(wave_start=3350,
                                               wave_end=6000,
                                               wave_bin=2,
                                               stype='science+standard')

    blue_onedspec.load_standard(standard_name)
    blue_onedspec.compute_sensitivity(smooth=True,
                                      slength=15,
                                      use_lowess=False)
    blue_onedspec.inspect_sensitivity()
    blue_onedspec.apply_flux_calibration()

    return blue_onedspec


# Line list
atlas_Hg_red = [5460.7348]
atlas_Ar_red = [
    6965.4307, 7067.2175, 7147.0416, 7272.9359, 7383.9805, 7503.8691,
    7635.1056, 7723.7599, 7948.1764, 8014.7857, 8115.3108, 8264.5225,
    8424.6475, 8521.4422, 8667.945, 9122.9674, 9224.4992, 9354.220, 9657.7863,
    9784.503, 10470.054
]
element_Hg_red = ['Hg'] * len(atlas_Hg_red)
element_Ar_red = ['Ar'] * len(atlas_Ar_red)

atlas_Ar_blue = [4158.590, 4200.674, 4300.101]
atlas_Hg_blue = [
    3650.153, 4046.563, 4077.8314, 4358.328, 5460.7348, 5769.598, 5790.663
]
#atlas_Zn_blue = [4680.14, 4722.15, 4810.53]
element_Ar_blue = ['Ar'] * len(atlas_Ar_blue)
element_Hg_blue = ['Hg'] * len(atlas_Hg_blue)
#element_Zn_blue = ['Zn'] * len(atlas_Zn_blue)

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
# Standard fits here (NOT just an array, need the )
#
standard_light_fits = fits.open(filename[(obstype == 'standard')
                                         & (frametype == 'light')][0])[1]
standard_flat_fits = fits.open(filename[(obstype == 'standard')
                                        & (frametype == 'flat')][0])[1]
standard_arc_fits = fits.open(filename[(obstype == 'standard')
                                       & (frametype == 'arc')][0])[1]
standard_name = standard_light_fits.header['OBJECT']

standard_red, standard_blue, standard_flat = extract_floyds(
    standard_light_fits, standard_flat_fits, standard_arc_fits)

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

science_red, science_blue, science_flat = extract_floyds(
    science_light_fits, science_flat_fits, science_arc_fits)

onedspec_blue = calibrate_blue(science_blue, standard_blue, standard_name)
onedspec_red = calibrate_red(science_red, standard_red, standard_name,
                             science_flat, standard_flat)

# Inspect
onedspec_blue.inspect_reduced_spectrum(
    wave_min=3000,
    wave_max=6000,
    stype='science',
    display=True,
    save_iframe=True,
    filename=filename[(obstype == 'science')
                      & (frametype == 'light')][0].split('.')[0] + '_blue')

onedspec_red.inspect_reduced_spectrum(
    wave_min=5000,
    wave_max=11000,
    stype='science',
    display=True,
    save_iframe=True,
    filename=filename[(obstype == 'science')
                      & (frametype == 'light')][0].split('.')[0] + '_red')

onedspec_blue.inspect_reduced_spectrum(
    wave_min=3000,
    wave_max=6000,
    stype='standard',
    display=True,
    save_iframe=True,
    filename=filename[(obstype == 'standard')
                      & (frametype == 'light')][0].split('.')[0] + '_blue')

onedspec_red.inspect_reduced_spectrum(
    wave_min=5000,
    wave_max=11000,
    stype='standard',
    display=True,
    save_iframe=True,
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
red_limit = 5200
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

flux_combined = np.concatenate(
    (flux_blue[wave_blue < red_limit], flux_weighted_combine,
     flux_red[wave_red > blue_limit]))
wave_combined = np.concatenate(
    (wave_blue[wave_blue < blue_limit], wave_red[wave_red >= blue_limit]))

plt.figure(1, figsize=(16, 8))
plt.clf()
plt.plot(wave_blue, flux_blue, color='blue', label='Blue arm (ASPIRED)')
plt.plot(wave_red, flux_red, color='red', label='Red arm (ASPIRED)')
plt.plot(wave_combined, flux_combined, color='black', label='Combined')
plt.xlim(min(wave_blue), max(wave_red))
plt.ylim(0, max(flux_combined))
plt.xlabel('Wavelength / A')
plt.ylabel('Flux / (erg / s / cm / cm / A)')
plt.legend()
plt.grid()
plt.tight_layout()
plt.title(science_name)
plt.show()
