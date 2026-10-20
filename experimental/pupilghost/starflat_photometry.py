"""
Star-flat test of the pupil ghost: measure its impact on aperture photometry and verify its mitigation.

Input: a dithered set of BANZAI-reduced (e91) frames of one star field (see submit_sophia_observation.py --dither), and
the ghost model written by model_pupilghost.py --fullres (pupilghost_<filter>.fits.fz, ghost_profile_<filter>.dat).

Every frame is measured in two versions that differ only in the flat field:
    A  as delivered (flat fielded with the ghost-contaminated flat)
    B  A * (1 + G), i.e. flat fielded with flat / (1 + G)
The sky in B again contains the ghost; it is removed by a fine background mesh plus a local annulus.

Analysis per filter and version:
  1. self calibration   m_ij = m_i + ZP_j + S(x, y),  S = radial bins about the ghost centre + non-radial quadratic
  2. template fit       S(r) = k * 1.0857 * G(r) + c;  expect k_A ~ 1 (impact), k_B ~ 0 (mitigation)
  3. paired-star test   per star: <m - ZP> where |G| is extremal  minus  <m - ZP> in the ghost-free zone
  4. ghost centre       fit the centre of the template to the A residuals; must match the flat-derived centre
  5. repeatability      per-star rms around m_i + ZP_j for A and B
  6. optional: second aperture (PSF change across the field) and colour split (--colour-catalogue)

Frames without an astrometric solution from BANZAI (WCSERR != 0, or no ra / dec in the catalogue) are registered
against the stars of the solved frames: pointing offset from a 2D histogram of star offsets, then a full linear WCS
fitted to the matched stars.

Naming: i indexes stars, j indexes frames. "Deviation" is a measured magnitude minus the model m_i + ZP_j.

Usage:
    python starflat_photometry.py --ghost-model-dir plots/pupilghost --outdir plots/starflat frames/*e91.fits.fz
"""
import argparse
import json
import logging
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sep
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points

from lcocommissioning.common.astrometric_registration import frame_has_valid_astrometry, register_frame
from scipy.optimize import minimize
from scipy.spatial import cKDTree

_log = logging.getLogger(__name__)

# Magnitude change per fractional flux change: d(mag) = 2.5 / ln(10) * d(flux) / flux = 1.0857 * d(flux) / flux
MAG_PER_FRACTIONAL_FLUX = 2.5 / np.log(10)

# Photometry versions: A uses the delivered (ghosted) flat, B the ghost-corrected flat
FLAT_VERSIONS = ('A', 'B')
VERSION_COLOURS = {'A': 'C3', 'B': 'C0'}

# Radial range in which to search for the strongest ghost feature for the paired-star test. The central core
# (r < 2.5') is excluded because too few stars land there.
PAIRED_TEST_SEARCH_RANGE_ARCMIN = (2.5, 12.5)


def photometry_key(flat_version, aperture_index):
    """ Column suffix of one photometry variant, e.g. 'A0' = flat version A, first aperture."""
    return f'{flat_version}{aperture_index}'


# ----------------------------------------------------------------------------------------------------------------------
# Ghost model
# ----------------------------------------------------------------------------------------------------------------------
class GhostModel:
    """ Pupil ghost model of one filter: full-resolution image of (1 + G) and the radial profile G(r)."""

    def __init__(self, ghost_model_dir, filter_name):
        self.filter_name = filter_name

        # Full-resolution correction image 1 + G and the ghost centre in unbinned pixel coordinates
        image_filename = os.path.join(ghost_model_dir, f'pupilghost_{filter_name}.fits.fz')
        with fits.open(image_filename) as hdulist:
            image_hdu = hdulist[1] if len(hdulist) > 1 else hdulist[0]
            self.one_plus_ghost_image = image_hdu.data.astype(np.float32)
            self.centre_x = image_hdu.header['GHOSTCX']
            self.centre_y = image_hdu.header['GHOSTCY']

        # Radial profile G(r), r in arcmin from the ghost centre
        profile = np.loadtxt(os.path.join(ghost_model_dir, f'ghost_profile_{filter_name}.dat'))
        self.profile_radius_arcmin = profile[:, 0]
        self.profile_ghost_fraction = profile[:, 1]

    def ghost_fraction_at_pixel(self, x_pixel, y_pixel):
        """ G at the given pixel positions, read from the full-resolution model image."""
        column = np.clip(np.round(x_pixel).astype(int), 0, self.one_plus_ghost_image.shape[1] - 1)
        row = np.clip(np.round(y_pixel).astype(int), 0, self.one_plus_ghost_image.shape[0] - 1)
        return self.one_plus_ghost_image[row, column] - 1

    def ghost_fraction_at_radius(self, radius_arcmin):
        """ G(r) from the radial profile; zero beyond the end of the profile."""
        return np.interp(radius_arcmin, self.profile_radius_arcmin, self.profile_ghost_fraction, right=0.)


# ----------------------------------------------------------------------------------------------------------------------
# Photometry of individual frames
# ----------------------------------------------------------------------------------------------------------------------
def read_frame(filename):
    """ Read science image, bad pixel mask and header of a BANZAI frame."""
    with fits.open(filename) as hdulist:
        science_hdu = hdulist['SCI']
        image = np.ascontiguousarray(science_hdu.data, dtype=np.float32)  # sep needs native byte order
        header = science_hdu.header
        if 'BPM' in hdulist:
            bad_pixel_mask = hdulist['BPM'].data > 0
        else:
            bad_pixel_mask = np.zeros(image.shape, bool)
    return image, bad_pixel_mask, header


def measure_frame(filename, ghost_model, args):
    """ Detect stars on version A, then do aperture photometry on versions A and B at identical positions.

    :param ghost_model: model used to create version B; None means B is identical to A (control run).
    :return: dict of per-star arrays (positions, fluxes per photometry variant), dict of frame metadata
    """
    image, bad_pixel_mask, header = read_frame(filename)
    pixel_scale_arcsec = header['PIXSCALE']
    exposure_time = header['EXPTIME']
    saturation_level = header.get('SATURATE', np.inf)

    # Treat non-finite pixels as bad
    bad_pixel_mask = bad_pixel_mask | ~np.isfinite(image)
    image[~np.isfinite(image)] = 0

    # Detect sources on the background-subtracted version A image
    background = sep.Background(image, mask=bad_pixel_mask, bw=args.background_mesh_size, bh=args.background_mesh_size,
                                fw=3, fh=3)
    background_subtracted = image - background
    sources = sep.extract(background_subtracted, args.detection_threshold, err=background.globalrms,
                          mask=bad_pixel_mask, minarea=9)
    is_clean_source = (sources['flag'] == 0) & (sources['b'] > 0.5)
    peak_above_background = sources['peak'] + background.globalback

    # Seeing: FWHM from the half-light radius of bright, unsaturated stars (FWHM = 2 * r50 for a Gaussian)
    is_psf_star = is_clean_source & (sources['peak'] > 20 * background.globalrms) & (peak_above_background < 0.5 * saturation_level)
    half_light_radius, _ = sep.flux_radius(background_subtracted, sources['x'][is_psf_star], sources['y'][is_psf_star],
                                           6 * sources['a'][is_psf_star], 0.5, subpix=5)
    fwhm_pixel = float(2 * np.median(half_light_radius))

    # Aperture and local background annulus radii scale with the seeing
    aperture_radii_pixel = [factor * fwhm_pixel for factor in args.aperture_radii_fwhm]
    annulus_inner_pixel = args.annulus_fwhm[0] * fwhm_pixel
    annulus_outer_pixel = args.annulus_fwhm[1] * fwhm_pixel

    # Select isolated, unsaturated stars away from the image edge:
    # no other detection within the outer annulus radius, peak below the linearity limit
    source_positions = np.c_[sources['x'], sources['y']]
    nearest_neighbour_distance = cKDTree(source_positions).query(source_positions, k=2)[0][:, 1]
    is_selected = is_clean_source
    is_selected &= nearest_neighbour_distance > annulus_outer_pixel
    is_selected &= peak_above_background < args.max_peak_fraction * saturation_level
    edge_margin = annulus_outer_pixel + 2
    is_selected &= (sources['x'] > edge_margin) & (sources['x'] < image.shape[1] - edge_margin)
    is_selected &= (sources['y'] > edge_margin) & (sources['y'] < image.shape[0] - edge_margin)
    x_pixel = sources['x'][is_selected]
    y_pixel = sources['y'][is_selected]

    star_measurements = {'x': x_pixel, 'y': y_pixel}

    # The two flat-field versions of this frame
    flat_version_images = {'A': image}
    if ghost_model is not None:
        flat_version_images['B'] = image * ghost_model.one_plus_ghost_image
    else:
        flat_version_images['B'] = image

    # Aperture photometry per flat version and aperture; the background mesh is recomputed for each version because
    # the sky of version B contains the ghost pattern
    for flat_version, version_image in flat_version_images.items():
        version_background = sep.Background(version_image, mask=bad_pixel_mask, bw=args.background_mesh_size,
                                            bh=args.background_mesh_size, fw=3, fh=3)
        version_background_subtracted = version_image - version_background
        pixel_noise = version_background.rms()
        for aperture_index, aperture_radius in enumerate(aperture_radii_pixel):
            flux, flux_error, photometry_flag = sep.sum_circle(
                version_background_subtracted, x_pixel, y_pixel, aperture_radius, err=pixel_noise, mask=bad_pixel_mask,
                bkgann=(annulus_inner_pixel, annulus_outer_pixel), gain=1.0)
            is_good_photometry = (photometry_flag == 0) & (flux > 0)
            key = photometry_key(flat_version, aperture_index)
            star_measurements[f'flux_{key}'] = np.where(is_good_photometry, flux / exposure_time, np.nan)
            star_measurements[f'flux_error_{key}'] = flux_error / exposure_time
        del version_background, version_background_subtracted, pixel_noise  # free memory before the next version

    # Sky coordinates for cross-matching stars between frames
    wcs = WCS(header)
    star_sky_positions = wcs.pixel_to_world(x_pixel, y_pixel)
    star_measurements['ra'] = star_sky_positions.ra.deg
    star_measurements['dec'] = star_sky_positions.dec.deg

    frame_centre = wcs.pixel_to_world(image.shape[1] / 2, image.shape[0] / 2)
    frame_metadata = dict(file=os.path.basename(filename), filter=header['FILTER'], exptime=exposure_time,
                          instrument=header.get('INSTRUME', 'camera'),
                          wcs=wcs, has_valid_astrometry=frame_has_valid_astrometry(filename, header),
                          airmass=header.get('AIRMASS', np.nan), mjd=header.get('MJD-OBS', np.nan),
                          fwhm_px=fwhm_pixel, fwhm_arcsec=fwhm_pixel * pixel_scale_arcsec, nstars=int(len(x_pixel)),
                          pixscale=pixel_scale_arcsec, shape=image.shape,
                          ra_c=frame_centre.ra.deg, dec_c=frame_centre.dec.deg)
    _log.info(f"{frame_metadata['file']}: {frame_metadata['filter']} FWHM {frame_metadata['fwhm_arcsec']:.2f}\" "
              f"{frame_metadata['nstars']} stars")
    return star_measurements, frame_metadata


# ----------------------------------------------------------------------------------------------------------------------
# Astrometric registration of frames without a BANZAI solution
# ----------------------------------------------------------------------------------------------------------------------
def register_frames(per_frame_measurements, frame_metadata, match_radius_arcsec):
    """ Give frames without a BANZAI astrometric solution a WCS fitted against the stars of the solved frames.

    If no frame is solved, all frames are registered relative to the first one; a consistent relative registration is
    all the cross-matching needs. Frames that cannot be registered are dropped.

    :return: per-frame measurements and metadata of the usable frames
    """
    is_solved = [metadata['has_valid_astrometry'] for metadata in frame_metadata]
    if all(is_solved):
        return per_frame_measurements, frame_metadata
    reference_frames = [index for index, solved in enumerate(is_solved) if solved]
    if not reference_frames:
        _log.warning("no frame has an astrometric solution; registering all frames relative to the first frame")
        reference_frames = [0]
        is_solved = [index == 0 for index in range(len(frame_metadata))]
    _, reference_positions = match_stars_between_frames([per_frame_measurements[index] for index in reference_frames],
                                                        match_radius_arcsec)
    _log.info(f"registering {is_solved.count(False)} frames without astrometric solution against "
              f"{len(reference_positions)} stars of {len(reference_frames)} solved frames")

    usable_measurements, usable_metadata = [], []
    for frame_measurements, metadata, solved in zip(per_frame_measurements, frame_metadata, is_solved):
        if not solved:
            registration = register_frame(frame_measurements['x'], frame_measurements['y'], metadata['wcs'], reference_positions)
            if registration is None:
                _log.warning(f"{metadata['file']}: registration failed, frame dropped")
                continue
            sky_positions = registration['wcs'].pixel_to_world(frame_measurements['x'], frame_measurements['y'])
            frame_measurements['ra'] = sky_positions.ra.deg
            frame_measurements['dec'] = sky_positions.dec.deg
            metadata['wcs'] = registration['wcs']
            metadata['registration'] = {key: value for key, value in registration.items() if key != 'wcs'}
            _log.info(f"{metadata['file']}: registered with {registration['number_of_matches']} stars, "
                      f"rms {registration['rms_arcsec']:.2f}\", pointing offset E {registration['pointing_offset_arcsec'][0]:+.1f}\" "
                      f"N {registration['pointing_offset_arcsec'][1]:+.1f}\"")
        usable_measurements.append(frame_measurements)
        usable_metadata.append(metadata)
    return usable_measurements, usable_metadata


def match_stars_between_frames(per_frame_measurements, match_radius_arcsec):
    """ Assign a common star id to detections across frames by sky position.

    Stars of the first frame start the master list; detections of later frames that do not match a known star
    within match_radius_arcsec are appended to it as new stars.

    :return: list (one entry per frame) of star id arrays, SkyCoord of the master list
    """
    master_list = None
    star_ids_per_frame = []
    for frame_measurements in per_frame_measurements:
        frame_positions = SkyCoord(frame_measurements['ra'], frame_measurements['dec'], unit='deg')
        if master_list is None:
            master_list = frame_positions
            star_ids_per_frame.append(np.arange(len(frame_positions)))
            continue

        nearest_master_index, separation, _ = frame_positions.match_to_catalog_sky(master_list)
        is_known_star = separation.arcsec < match_radius_arcsec
        number_of_new_stars = int((~is_known_star).sum())
        new_star_ids = np.arange(len(master_list), len(master_list) + number_of_new_stars)

        star_ids = np.where(is_known_star, nearest_master_index, -1)
        star_ids[~is_known_star] = new_star_ids
        star_ids_per_frame.append(star_ids)

        master_list = SkyCoord(np.r_[master_list.ra.deg, frame_positions.ra.deg[~is_known_star]],
                               np.r_[master_list.dec.deg, frame_positions.dec.deg[~is_known_star]], unit='deg')
    return star_ids_per_frame, master_list


def build_measurement_table(per_frame_measurements, frame_metadata, ghost_model, args):
    """ Combine all frames into one flat table with one row per (star, frame) measurement.

    Columns: star, frame, x, y, ra, dec, r (arcmin from the ghost centre), G (ghost fraction at the star position),
    and per photometry variant: flux_*, flux_error_*, mag_*, mag_error_*.
    """
    star_ids_per_frame, master_list = match_stars_between_frames(per_frame_measurements, args.match_radius)

    # Concatenate the per-frame arrays
    columns = {name: [] for name in ['star', 'frame', 'x', 'y', 'ra', 'dec']}
    for frame_index, (frame_measurements, star_ids) in enumerate(zip(per_frame_measurements, star_ids_per_frame)):
        columns['star'].append(star_ids)
        columns['frame'].append(np.full(len(star_ids), frame_index))
        for name in ('x', 'y', 'ra', 'dec'):
            columns[name].append(frame_measurements[name])
        for name in frame_measurements:
            if name.startswith('flux_'):
                columns.setdefault(name, []).append(frame_measurements[name])
    measurements = {name: np.concatenate(values) for name, values in columns.items()}

    # Position relative to the pupil ghost
    pixel_scale_arcsec = frame_metadata[0]['pixscale']
    measurements['r'] = np.hypot(measurements['x'] - ghost_model.centre_x,
                                 measurements['y'] - ghost_model.centre_y) * pixel_scale_arcsec / 60
    measurements['G'] = ghost_model.ghost_fraction_at_pixel(measurements['x'], measurements['y'])

    # Instrumental magnitudes; an error floor accounts for systematics not captured by photon noise
    flux_columns = [name for name in measurements if name.startswith('flux_') and not name.startswith('flux_error_')]
    for flux_column in flux_columns:
        key = flux_column[len('flux_'):]
        flux = measurements[flux_column]
        flux_error = measurements[f'flux_error_{key}']
        measurements[f'mag_{key}'] = -2.5 * np.log10(flux)
        measurements[f'mag_error_{key}'] = np.sqrt((MAG_PER_FRACTIONAL_FLUX * flux_error / flux) ** 2 + args.error_floor ** 2)

    # Keep good, high signal-to-noise measurements of stars that were observed often enough
    signal_to_noise = measurements['flux_A0'] / measurements['flux_error_A0']
    is_kept = np.isfinite(measurements['mag_A0']) & np.isfinite(measurements['mag_B0'])
    is_kept &= signal_to_noise > args.min_snr
    observations_per_star = np.bincount(measurements['star'][is_kept], minlength=measurements['star'].max() + 1)
    is_kept &= observations_per_star[measurements['star']] >= args.min_observations
    measurements = {name: values[is_kept] for name, values in measurements.items()}

    # Renumber the remaining stars 0 .. N-1
    _, measurements['star'] = np.unique(measurements['star'], return_inverse=True)
    return measurements, master_list


# ----------------------------------------------------------------------------------------------------------------------
# Self calibration (star flat)
# ----------------------------------------------------------------------------------------------------------------------
def nonradial_quadratic_terms(measurements, image_shape):
    """ Non-radial low-order illumination terms x, y, x^2 - y^2, xy (normalised coordinates).

    x^2 + y^2 is omitted on purpose: it is radial and therefore degenerate with the radial star-flat bins.
    """
    x_normalised = (measurements['x'] - image_shape[1] / 2) / (image_shape[1] / 2)
    y_normalised = (measurements['y'] - image_shape[0] / 2) / (image_shape[0] / 2)
    return np.c_[x_normalised, y_normalised, x_normalised ** 2 - y_normalised ** 2, x_normalised * y_normalised]


def weighted_mean_per_group(group_index, values, weights, number_of_groups):
    """ Weighted mean of values within each group (star or frame)."""
    # Zero-weight entries may hold NaN (invalid measurements); keep them out of the sum, since 0 * NaN = NaN
    values_with_weight = np.where(weights > 0, values, 0.)
    weighted_sum = np.bincount(group_index, weights * values_with_weight, number_of_groups)
    sum_of_weights = np.bincount(group_index, weights, number_of_groups)
    return weighted_sum / np.maximum(sum_of_weights, 1e-30)


def self_calibrate(measurements, magnitude, magnitude_error, image_shape, radial_bin_edges, normalisation_zone,
                   fit_star_flat=True, max_iterations=30, clip_sigma=5.):
    """ Fit  magnitude_ij = m_i + ZP_j + S(x, y)  by alternating weighted least squares with sigma clipping.

    m_i   mean magnitude of star i
    ZP_j  zero point of frame j (transparency, exposure); ZP of the first frame is fixed to 0
    S     star flat: radial bins about the ghost centre + non-radial quadratic terms; S = 0 on average in the
          normalisation zone. With fit_star_flat=False, S = 0 and only m_i + ZP_j are fitted.
    """
    star_index = measurements['star']
    frame_index = measurements['frame']
    number_of_stars = star_index.max() + 1
    number_of_frames = frame_index.max() + 1
    # Measurements without a valid magnitude (e.g. a failed second aperture) take no part in the fit
    is_valid = np.isfinite(magnitude) & np.isfinite(magnitude_error) & (magnitude_error > 0)
    magnitude = np.where(is_valid, magnitude, 0.)
    magnitude_error = np.where(is_valid, magnitude_error, 1.)
    inverse_variance = np.where(is_valid, 1 / magnitude_error ** 2, 0.)

    # Design matrix of the star flat: one indicator column per radial bin, then the non-radial terms
    radial_bin_index = np.digitize(measurements['r'], radial_bin_edges) - 1
    number_of_radial_bins = len(radial_bin_edges) - 1
    is_inside_radial_range = (radial_bin_index >= 0) & (radial_bin_index < number_of_radial_bins)
    radial_bin_indicator = np.eye(number_of_radial_bins)[np.clip(radial_bin_index, 0, number_of_radial_bins - 1)]
    radial_bin_indicator *= is_inside_radial_range[:, None]
    star_flat_design = np.c_[radial_bin_indicator, nonradial_quadratic_terms(measurements, image_shape)]

    radial_bin_centres = 0.5 * (radial_bin_edges[1:] + radial_bin_edges[:-1])
    is_normalisation_bin = (radial_bin_centres > normalisation_zone[0]) & (radial_bin_centres < normalisation_zone[1])

    # Start values
    is_used = is_valid.copy()
    frame_zeropoint = np.zeros(number_of_frames)
    star_flat = np.zeros(len(magnitude))
    star_flat_coefficients = np.zeros(star_flat_design.shape[1])
    star_flat_covariance = np.zeros((star_flat_design.shape[1], star_flat_design.shape[1]))

    for iteration in range(max_iterations):
        weights = inverse_variance * is_used

        # 1. star magnitudes, given zero points and star flat
        star_magnitude = weighted_mean_per_group(star_index, magnitude - frame_zeropoint[frame_index] - star_flat,
                                                 weights, number_of_stars)
        # 2. frame zero points, given star magnitudes and star flat; fix the first frame to 0
        frame_zeropoint = weighted_mean_per_group(frame_index, magnitude - star_magnitude[star_index] - star_flat,
                                                  weights, number_of_frames)
        frame_zeropoint -= frame_zeropoint[0]

        # 3. star flat, given star magnitudes and zero points (weighted linear least squares)
        if fit_star_flat:
            residual_without_star_flat = magnitude - star_magnitude[star_index] - frame_zeropoint[frame_index]
            sqrt_weights = np.sqrt(weights[is_used])
            star_flat_coefficients, *_ = np.linalg.lstsq(star_flat_design[is_used] * sqrt_weights[:, None],
                                                         residual_without_star_flat[is_used] * sqrt_weights, rcond=None)
            # define S = 0 in the ghost-free normalisation zone; the offset moves into m_i in the next iteration
            star_flat_coefficients[:number_of_radial_bins] -= np.mean(
                star_flat_coefficients[:number_of_radial_bins][is_normalisation_bin])
            star_flat = star_flat_design @ star_flat_coefficients

        # 4. sigma clipping on normalised residuals, with a robust (MAD) scale of at least 1
        residual = magnitude - star_magnitude[star_index] - frame_zeropoint[frame_index] - star_flat
        normalised_residual = residual / magnitude_error
        robust_scale = 1.4826 * np.median(np.abs(normalised_residual[is_used]))
        is_used_next = is_valid & (np.abs(normalised_residual) < clip_sigma * max(robust_scale, 1))
        if iteration > 3 and np.array_equal(is_used_next, is_used):
            break
        is_used = is_used_next

    # Goodness of fit; uncertainties of the star-flat coefficients are scaled up if the reduced chi^2 exceeds 1
    number_of_free_parameters = number_of_stars + number_of_frames + (star_flat_design.shape[1] if fit_star_flat else 0)
    reduced_chi2 = np.sum((residual[is_used] / magnitude_error[is_used]) ** 2) / max(is_used.sum() - number_of_free_parameters, 1)
    if fit_star_flat:
        weighted_design = star_flat_design[is_used] * np.sqrt(inverse_variance[is_used])[:, None]
        star_flat_covariance = np.linalg.pinv(weighted_design.T @ weighted_design) * max(reduced_chi2, 1)

    measurements_per_radial_bin = np.bincount(np.clip(radial_bin_index[is_used], 0, number_of_radial_bins - 1),
                                              minlength=number_of_radial_bins)
    if fit_star_flat:
        radial_star_flat_variance = np.where(measurements_per_radial_bin > 0,
                                             np.diag(star_flat_covariance)[:number_of_radial_bins], np.nan)
    else:
        radial_star_flat_variance = np.full(number_of_radial_bins, np.nan)

    # Per-star repeatability: rms of the residuals of each star
    sum_of_squared_residuals = np.bincount(star_index[is_used], residual[is_used] ** 2, number_of_stars)
    used_measurements_per_star = np.bincount(star_index[is_used], minlength=number_of_stars)
    star_rms = np.sqrt(sum_of_squared_residuals / np.maximum(used_measurements_per_star, 1))

    return dict(star_magnitude=star_magnitude, frame_zeropoint=frame_zeropoint, star_flat=star_flat,
                residual=np.where(is_valid, residual, np.nan), is_used=is_used, radial_bin_index=radial_bin_index,
                radial_bin_centres=radial_bin_centres,
                radial_star_flat=star_flat_coefficients[:number_of_radial_bins],
                radial_star_flat_error=np.sqrt(np.abs(radial_star_flat_variance)),
                nonradial_coefficients=star_flat_coefficients[number_of_radial_bins:],
                measurements_per_radial_bin=measurements_per_radial_bin, reduced_chi2=reduced_chi2,
                star_rms=star_rms)


def fit_ghost_template(calibration, measurements, max_radius_arcmin, min_measurements_per_bin=30):
    """ Fit the radial star flat with the ghost prediction:  S(r_bin) = k * 1.0857 * <G>_bin + c.

    <G>_bin is the mean ghost fraction of the measurements in each radial bin. k = 1 means the photometric error
    is fully explained by the ghost model, k = 0 means no ghost signature.
    """
    number_of_radial_bins = len(calibration['radial_bin_centres'])
    radial_bin_index = np.clip(calibration['radial_bin_index'], 0, number_of_radial_bins - 1)
    is_used = calibration['is_used']

    # Mean ghost fraction per radial bin, from the actual star positions in that bin
    ghost_sum_per_bin = np.bincount(radial_bin_index[is_used], measurements['G'][is_used], number_of_radial_bins)
    measurements_per_bin = np.bincount(radial_bin_index[is_used], minlength=number_of_radial_bins)
    mean_ghost_per_bin = ghost_sum_per_bin / np.maximum(measurements_per_bin, 1)

    # Use well-populated bins with a valid uncertainty
    radial_star_flat = calibration['radial_star_flat']
    radial_star_flat_error = calibration['radial_star_flat_error']
    is_fitted_bin = (calibration['measurements_per_radial_bin'] >= min_measurements_per_bin)
    is_fitted_bin &= calibration['radial_bin_centres'] < max_radius_arcmin
    is_fitted_bin &= np.isfinite(radial_star_flat_error) & (radial_star_flat_error > 0)

    # Weighted linear least squares for (k, c)
    design = np.c_[MAG_PER_FRACTIONAL_FLUX * mean_ghost_per_bin[is_fitted_bin], np.ones(is_fitted_bin.sum())]
    bin_weight = 1 / radial_star_flat_error[is_fitted_bin]
    weighted_design = design * bin_weight[:, None]
    (slope_k, offset_c), *_ = np.linalg.lstsq(weighted_design, radial_star_flat[is_fitted_bin] * bin_weight, rcond=None)
    chi2 = np.sum(((radial_star_flat[is_fitted_bin] - design @ [slope_k, offset_c]) * bin_weight) ** 2)
    degrees_of_freedom = max(is_fitted_bin.sum() - 2, 1)
    parameter_covariance = np.linalg.inv(weighted_design.T @ weighted_design) * max(chi2 / degrees_of_freedom, 1)

    # Reference: chi^2 of the null hypothesis "S(r) = constant", i.e. no ghost signature at all
    constant_level = np.sum(radial_star_flat[is_fitted_bin] * bin_weight ** 2) / np.sum(bin_weight ** 2)
    chi2_no_ghost = np.sum(((radial_star_flat[is_fitted_bin] - constant_level) * bin_weight) ** 2)

    is_inside_ghost = is_fitted_bin & (calibration['radial_bin_centres'] < 12.5)
    return dict(k=float(slope_k), k_err=float(np.sqrt(parameter_covariance[0, 0])), c=float(offset_c),
                chi2red=float(chi2 / degrees_of_freedom),
                chi2red_null=float(chi2_no_ghost / max(is_fitted_bin.sum() - 1, 1)),
                mean_ghost_per_bin=mean_ghost_per_bin, is_fitted_bin=is_fitted_bin,
                rms_Sr_pct=float(np.std(radial_star_flat[is_inside_ghost]) / MAG_PER_FRACTIONAL_FLUX * 100))


def paired_star_test(measurements, calibration, ghost_model, magnitude, normalisation_zone, zone_half_width=0.75,
                     number_of_bootstrap_samples=500, random_seed=1):
    """ Model-independent test: for each star, mean (magnitude - ZP) where |G| is largest minus the same in the
    ghost-free normalisation zone. Expected: 1.0857 * (G_ghost_zone - G_ghost_free_zone) for version A, 0 for B.
    """
    # Ghost zone: radius of the strongest ghost feature (largest |G|) within the search range
    candidate_radii = np.arange(*PAIRED_TEST_SEARCH_RANGE_ARCMIN, 0.05)
    candidate_ghost = ghost_model.ghost_fraction_at_radius(candidate_radii)
    ghost_zone_radius = candidate_radii[np.argmax(np.abs(candidate_ghost))]
    is_in_ghost_zone = np.abs(measurements['r'] - ghost_zone_radius) < zone_half_width
    is_in_ghost_free_zone = (measurements['r'] > normalisation_zone[0]) & (measurements['r'] < normalisation_zone[1])

    # Mean zero-point corrected magnitude and mean ghost fraction of every star in each zone
    zeropoint_corrected_magnitude = magnitude - calibration['frame_zeropoint'][measurements['frame']]
    number_of_stars = measurements['star'].max() + 1
    zone_statistics = {}
    for zone_name, is_in_zone in (('ghost', is_in_ghost_zone), ('ghost_free', is_in_ghost_free_zone)):
        is_selected = is_in_zone & calibration['is_used']
        star_index = measurements['star'][is_selected]
        count = np.bincount(star_index, minlength=number_of_stars)
        zone_statistics[zone_name] = dict(
            count=count,
            mean_magnitude=np.bincount(star_index, zeropoint_corrected_magnitude[is_selected], number_of_stars) / np.maximum(count, 1),
            mean_ghost=np.bincount(star_index, measurements['G'][is_selected], number_of_stars) / np.maximum(count, 1))

    # Stars observed in both zones
    is_in_both_zones = (zone_statistics['ghost']['count'] > 0) & (zone_statistics['ghost_free']['count'] > 0)
    magnitude_difference = (zone_statistics['ghost']['mean_magnitude']
                            - zone_statistics['ghost_free']['mean_magnitude'])[is_in_both_zones]
    predicted_difference = MAG_PER_FRACTIONAL_FLUX * (zone_statistics['ghost']['mean_ghost']
                                                      - zone_statistics['ghost_free']['mean_ghost'])[is_in_both_zones]

    # Uncertainty of the median by bootstrap resampling of the stars
    random_generator = np.random.default_rng(random_seed)
    if len(magnitude_difference) > 2:
        bootstrap_medians = [np.median(random_generator.choice(magnitude_difference, len(magnitude_difference)))
                             for _ in range(number_of_bootstrap_samples)]
    else:
        bootstrap_medians = [np.nan]
    has_stars = len(magnitude_difference) > 0
    return dict(r_zone=float(ghost_zone_radius), nstars=int(is_in_both_zones.sum()),
                median_dmag=float(np.median(magnitude_difference)) if has_stars else np.nan,
                err=float(np.std(bootstrap_medians)),
                predicted=float(np.median(predicted_difference)) if has_stars else np.nan,
                magnitude_difference=magnitude_difference)


def fit_ghost_centre(measurements, calibration_without_star_flat, ghost_model, pixel_scale_arcsec, image_shape):
    """ Fit the centre of the ghost template to the residuals of the m_i + ZP_j fit (no star flat).

    For a trial centre, the residuals are fitted linearly with  k * 1.0857 * G(r) + constant + linear x/y terms.
    The centre with the lowest sum of squared residuals wins. Several starting points guard against local minima.
    """
    is_used = calibration_without_star_flat['is_used']
    residual = calibration_without_star_flat['residual'][is_used]
    x_pixel = measurements['x'][is_used]
    y_pixel = measurements['y'][is_used]

    def sum_of_squared_residuals(trial_centre):
        radius_arcmin = np.hypot(x_pixel - trial_centre[0], y_pixel - trial_centre[1]) * pixel_scale_arcsec / 60
        ghost_fraction = ghost_model.ghost_fraction_at_radius(radius_arcmin)
        design = np.c_[MAG_PER_FRACTIONAL_FLUX * ghost_fraction, np.ones_like(ghost_fraction),
                       x_pixel / image_shape[1], y_pixel / image_shape[0]]
        coefficients, *_ = np.linalg.lstsq(design, residual, rcond=None)
        return np.sum((residual - design @ coefficients) ** 2)

    starting_points = [(ghost_model.centre_x, ghost_model.centre_y),
                       (image_shape[1] / 2, image_shape[0] / 2),
                       (ghost_model.centre_x + 300, ghost_model.centre_y),
                       (ghost_model.centre_x, ghost_model.centre_y + 300)]
    best_result = None
    for starting_point in starting_points:
        result = minimize(sum_of_squared_residuals, starting_point, method='Nelder-Mead', options=dict(xatol=1, fatol=1e-10))
        if best_result is None or result.fun < best_result.fun:
            best_result = result
    return float(best_result.x[0]), float(best_result.x[1])


def colour_split_test(measurements, args, image_shape, radial_bin_edges):
    """ Repeat the template fit separately for blue and red stars, split at the median colour of an external
    catalogue (CSV with columns ra,dec,color). A colour-independent k shows the effect is achromatic for stars.
    """
    catalogue = np.genfromtxt(args.colour_catalogue, delimiter=',', names=True)

    # Mean position of every star, matched to the catalogue
    number_of_stars = measurements['star'].max() + 1
    measurements_per_star = np.bincount(measurements['star'], minlength=number_of_stars)
    star_ra = np.bincount(measurements['star'], measurements['ra'], number_of_stars) / measurements_per_star
    star_dec = np.bincount(measurements['star'], measurements['dec'], number_of_stars) / measurements_per_star
    catalogue_index, separation, _ = SkyCoord(star_ra, star_dec, unit='deg').match_to_catalog_sky(
        SkyCoord(catalogue['ra'], catalogue['dec'], unit='deg'))
    star_colour = np.where(separation.arcsec < args.match_radius, catalogue['color'][catalogue_index], np.nan)
    median_colour = np.nanmedian(star_colour)

    results = dict(split_colour=float(median_colour))
    for colour_group, is_star_in_group in (('blue', star_colour < median_colour), ('red', star_colour >= median_colour)):
        # Subset of the measurement table, with stars and frames renumbered
        is_in_group = is_star_in_group[measurements['star']]
        group_measurements = {name: values[is_in_group] for name, values in measurements.items()}
        _, group_measurements['star'] = np.unique(group_measurements['star'], return_inverse=True)
        _, group_measurements['frame'] = np.unique(group_measurements['frame'], return_inverse=True)
        for flat_version in FLAT_VERSIONS:
            key = photometry_key(flat_version, 0)
            calibration = self_calibrate(group_measurements, group_measurements[f'mag_{key}'],
                                         group_measurements[f'mag_error_{key}'], image_shape, radial_bin_edges,
                                         args.normalisation_zone)
            template = fit_ghost_template(calibration, group_measurements, args.max_template_radius)
            results[f'{colour_group}_{flat_version}'] = dict(k=template['k'], k_err=template['k_err'],
                                                            nstars=int(group_measurements['star'].max() + 1))
    return results


def evaluate_success_criteria(results, is_control_run=False):
    """ Pass/fail of the pre-defined success criteria."""
    version_a = results['A']
    version_b = results['B']
    if is_control_run:
        return {'control: no ghost signature, |k_A| < 0.15': bool(abs(version_a['k']) < 0.15)}

    predicted_paired_difference = abs(version_a['paired_predicted_A_mmag'])
    paired_a_deviation_from_prediction = abs(version_a['paired_dmag_mmag'] - version_a['paired_predicted_A_mmag'])
    criteria = {
        'impact detected: k_A = 1.00 +/- 0.15': bool(abs(version_a['k'] - 1) < 0.15),
        'impact significant: k_A / err > 10': bool(version_a['k_snr'] > 10),
        'mitigation: |k_B| < 0.15': bool(abs(version_b['k']) < 0.15),
        'ghost centre within 1 arcmin of flat-derived centre': bool(results['ghost_centre_offset_arcmin'] < 1),
        'repeatability improved: rms_B < rms_A': bool(version_b['star_rms_median_mmag'] < version_a['star_rms_median_mmag']),
        'paired test A agrees with prediction (2 sigma or 20%)': bool(
            paired_a_deviation_from_prediction < max(2 * version_a['paired_err_mmag'], 0.2 * predicted_paired_difference)),
        'paired test B consistent with zero (< 2 sigma or < 15% of prediction)': bool(
            abs(version_b['paired_dmag_mmag']) < max(2 * version_b['paired_err_mmag'], 0.15 * predicted_paired_difference)),
    }
    if version_a.get('k_aperture2') is not None:
        criteria['aperture independent: |k_A(ap1) - k_A(ap2)| < 0.1'] = bool(
            abs(version_a['k'] - version_a['k_aperture2']) < 0.1)
    return criteria


# ----------------------------------------------------------------------------------------------------------------------
# Analysis of one filter
# ----------------------------------------------------------------------------------------------------------------------
def analyse_filter(filter_name, files, args):
    """ Photometry, self calibration, tests and plots for all frames of one filter. Returns the results dict."""
    # The ghost template is the filter's own model, or another filter's model for a control run
    template_filter = args.template or filter_name
    is_control_run = template_filter != filter_name
    if not os.path.exists(os.path.join(args.ghost_model_dir, f'pupilghost_{template_filter}.fits.fz')):
        _log.error(f"no ghost model for template filter {template_filter} in {args.ghost_model_dir}")
        return None
    ghost_model = GhostModel(args.ghost_model_dir, template_filter)
    if is_control_run:
        _log.info(f"{filter_name}: control run, testing against the {template_filter} ghost template; "
                  f"version B is not corrected")

    # Photometry of every frame; in a control run, version B stays uncorrected
    per_frame_measurements, frame_metadata = [], []
    for filename in files:
        frame_measurements, metadata = measure_frame(filename, None if is_control_run else ghost_model, args)
        per_frame_measurements.append(frame_measurements)
        frame_metadata.append(metadata)
    per_frame_measurements, frame_metadata = register_frames(per_frame_measurements, frame_metadata, args.match_radius)
    files = [metadata['file'] for metadata in frame_metadata]
    measurements, _ = build_measurement_table(per_frame_measurements, frame_metadata, ghost_model, args)
    image_shape = frame_metadata[0]['shape']
    pixel_scale_arcsec = frame_metadata[0]['pixscale']
    number_of_stars = int(measurements['star'].max() + 1)
    _log.info(f"{filter_name}: {number_of_stars} stars, {len(measurements['star'])} measurements in {len(files)} frames")

    # Radial bins from the ghost centre out to the detector corners
    max_radius_arcmin = np.hypot(*image_shape) / 2 * pixel_scale_arcsec / 60
    radial_bin_edges = np.arange(0, max_radius_arcmin + args.radial_bin_width, args.radial_bin_width)

    results = dict(filter=filter_name, template=template_filter, nframes=len(files), nstars=number_of_stars,
                   nmeas=int(len(measurements['star'])), ghost_centre_model=[ghost_model.centre_x, ghost_model.centre_y],
                   fwhm_arcsec_median=float(np.median([metadata['fwhm_arcsec'] for metadata in frame_metadata])))

    # Self calibration with and without star flat, template fit and paired test, per flat version and aperture
    calibrations_with_star_flat = {}
    calibrations_without_star_flat = {}
    template_fits = {}
    paired_tests = {}
    for flat_version in FLAT_VERSIONS:
        for aperture_index in range(len(args.aperture_radii_fwhm)):
            key = photometry_key(flat_version, aperture_index)
            magnitude = measurements[f'mag_{key}']
            magnitude_error = measurements[f'mag_error_{key}']
            calibrations_with_star_flat[key] = self_calibrate(measurements, magnitude, magnitude_error, image_shape,
                                                              radial_bin_edges, args.normalisation_zone)
            template_fits[key] = fit_ghost_template(calibrations_with_star_flat[key], measurements, args.max_template_radius)
            calibrations_without_star_flat[key] = self_calibrate(measurements, magnitude, magnitude_error, image_shape,
                                                                 radial_bin_edges, args.normalisation_zone,
                                                                 fit_star_flat=False)

        main_key = photometry_key(flat_version, 0)
        second_key = photometry_key(flat_version, 1)
        has_second_aperture = len(args.aperture_radii_fwhm) > 1
        paired_tests[flat_version] = paired_star_test(measurements, calibrations_with_star_flat[main_key], ghost_model,
                                                      measurements[f'mag_{main_key}'], args.normalisation_zone)
        star_rms = calibrations_without_star_flat[main_key]['star_rms']
        template = template_fits[main_key]
        paired = paired_tests[flat_version]
        results[flat_version] = dict(
            k=template['k'], k_err=template['k_err'], k_snr=template['k'] / template['k_err'],
            chi2red_template=template['chi2red'], chi2red_noghost=template['chi2red_null'],
            radial_rms_pct=template['rms_Sr_pct'],
            k_aperture2=template_fits[second_key]['k'] if has_second_aperture else None,
            k_aperture2_err=template_fits[second_key]['k_err'] if has_second_aperture else None,
            paired_dmag_mmag=paired['median_dmag'] * 1000, paired_err_mmag=paired['err'] * 1000,
            paired_predicted_A_mmag=paired['predicted'] * 1000, paired_nstars=paired['nstars'],
            paired_r_zone=paired['r_zone'],
            star_rms_median_mmag=float(np.median(star_rms) * 1000),
            selfcal_chi2red=float(calibrations_with_star_flat[main_key]['reduced_chi2']),
            # radial star flat of the main aperture (bins with enough measurements), for comparisons between data sets
            radial_star_flat=dict(
                radius_arcmin=calibrations_with_star_flat[main_key]['radial_bin_centres'][template['is_fitted_bin']].tolist(),
                star_flat_mmag=(calibrations_with_star_flat[main_key]['radial_star_flat'][template['is_fitted_bin']] * 1000).tolist(),
                error_mmag=(calibrations_with_star_flat[main_key]['radial_star_flat_error'][template['is_fitted_bin']] * 1000).tolist()))

    # Does the photometric error pattern sit where the flat-field ghost sits?
    fitted_centre_x, fitted_centre_y = fit_ghost_centre(measurements, calibrations_without_star_flat['A0'], ghost_model,
                                                        pixel_scale_arcsec, image_shape)
    results['ghost_centre_fit_A'] = [fitted_centre_x, fitted_centre_y]
    results['ghost_centre_offset_arcmin'] = float(np.hypot(fitted_centre_x - ghost_model.centre_x,
                                                           fitted_centre_y - ghost_model.centre_y) * pixel_scale_arcsec / 60)

    if args.colour_catalogue:
        results['colour'] = colour_split_test(measurements, args, image_shape, radial_bin_edges)

    results['criteria'] = evaluate_success_criteria(results, is_control_run)
    camera = frame_metadata[0]['instrument']
    results['camera'] = camera
    plot_summary(filter_name, camera, measurements, calibrations_with_star_flat, calibrations_without_star_flat, template_fits,
                 paired_tests, ghost_model, results, args, image_shape)
    plot_deviations(filter_name, camera, measurements, calibrations_without_star_flat, ghost_model, args, image_shape,
                    radial_bin_edges)
    return results


# ----------------------------------------------------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------------------------------------------------
def binned_detector_map(x_pixel, y_pixel, values, image_shape, cells_per_axis, min_count):
    """ Mean of values in a cells_per_axis x cells_per_axis grid over the detector; NaN where count < min_count."""
    cell_column = np.clip((x_pixel / image_shape[1] * cells_per_axis).astype(int), 0, cells_per_axis - 1)
    cell_row = np.clip((y_pixel / image_shape[0] * cells_per_axis).astype(int), 0, cells_per_axis - 1)
    cell_index = cell_row * cells_per_axis + cell_column
    count = np.bincount(cell_index, minlength=cells_per_axis ** 2)
    mean_value = np.bincount(cell_index, values, cells_per_axis ** 2) / np.maximum(count, 1)
    mean_value[count < min_count] = np.nan
    return mean_value.reshape(cells_per_axis, cells_per_axis)


def plot_summary(filter_name, camera, measurements, calibrations_with_star_flat, calibrations_without_star_flat, template_fits,
                 paired_tests, ghost_model, results, args, image_shape):
    """ Six-panel summary: radial star flat, residual maps A and B, paired test, repeatability, radial coverage."""
    fig, axes = plt.subplots(2, 3, figsize=(19, 11))
    radius_grid = np.linspace(0, calibrations_with_star_flat['A0']['radial_bin_centres'].max(), 400)

    # Panel 1: radial star flat of A and B with the ghost prediction
    axis = axes[0, 0]
    for flat_version in FLAT_VERSIONS:
        key = photometry_key(flat_version, 0)
        calibration = calibrations_with_star_flat[key]
        template = template_fits[key]
        is_fitted = template['is_fitted_bin']
        axis.errorbar(calibration['radial_bin_centres'][is_fitted], calibration['radial_star_flat'][is_fitted] * 1000,
                      calibration['radial_star_flat_error'][is_fitted] * 1000, fmt='o', ms=4,
                      color=VERSION_COLOURS[flat_version],
                      label=f"{flat_version}: k = {template['k']:.3f} ± {template['k_err']:.3f}")
    predicted_star_flat = MAG_PER_FRACTIONAL_FLUX * ghost_model.ghost_fraction_at_radius(radius_grid) + template_fits['A0']['c']
    axis.plot(radius_grid, predicted_star_flat * 1000, 'k-', lw=1, label='prediction 1.086·G(r)')
    axis.axvspan(*args.normalisation_zone, color='grey', alpha=0.15)
    axis.axhline(0, color='k', lw=0.5)
    axis.set_xlabel('r from ghost centre [arcmin]')
    axis.set_ylabel('star flat S(r) [mmag]  (+ = stars too faint)')
    axis.set_title(f'{camera} {filter_name}: radial star flat')
    axis.legend(fontsize=8)
    axis.grid(alpha=0.3)

    # Panels 2 and 3: 2D maps of the residuals of the m_i + ZP_j fit (no star-flat model)
    colour_limit_mmag = max(5, 1000 * MAG_PER_FRACTIONAL_FLUX * np.ptp(ghost_model.profile_ghost_fraction) / 2)
    for panel_index, flat_version in enumerate(FLAT_VERSIONS):
        axis = axes[0, 1 + panel_index]
        calibration = calibrations_without_star_flat[photometry_key(flat_version, 0)]
        is_used = calibration['is_used']
        residual_map = binned_detector_map(measurements['x'][is_used], measurements['y'][is_used],
                                           calibration['residual'][is_used], image_shape, cells_per_axis=24, min_count=5)
        image = axis.imshow(residual_map * 1000, origin='lower', cmap='RdBu_r', vmin=-colour_limit_mmag,
                            vmax=colour_limit_mmag, extent=(0, image_shape[1], 0, image_shape[0]))
        axis.plot(ghost_model.centre_x, ghost_model.centre_y, 'k+', ms=12)
        if flat_version == 'A':
            axis.plot(*results['ghost_centre_fit_A'], 'gx', ms=10)
        axis.set_title(f'{flat_version}: photometric residual map [mmag] (+ model ghost centre, x fitted)')
        plt.colorbar(image, ax=axis, fraction=0.046)

    # Panel 4: paired-star test histograms
    axis = axes[1, 0]
    histogram_bins = np.linspace(-0.08, 0.08, 81)
    for flat_version in FLAT_VERSIONS:
        paired = paired_tests[flat_version]
        axis.hist(paired['magnitude_difference'], bins=histogram_bins, histtype='step',
                  color=VERSION_COLOURS[flat_version], lw=1.5,
                  label=f"{flat_version}: median {paired['median_dmag'] * 1000:+.1f} ± {paired['err'] * 1000:.1f} mmag "
                        f"(N={paired['nstars']})")
    axis.axvline(paired_tests['A']['predicted'], color='k', ls='--',
                 label=f"prediction A {paired_tests['A']['predicted'] * 1000:+.1f} mmag")
    axis.axvline(0, color='k', lw=0.5)
    axis.set_xlabel(f"Δm per star: r≈{paired_tests['A']['r_zone']:.1f}' minus ghost-free zone [mag]")
    axis.legend(fontsize=8)
    axis.set_title('paired-star test')

    # Panel 5: per-star repeatability
    axis = axes[1, 1]
    histogram_bins = np.linspace(0, 0.04, 81)
    for flat_version in FLAT_VERSIONS:
        star_rms = calibrations_without_star_flat[photometry_key(flat_version, 0)]['star_rms']
        axis.hist(star_rms, bins=histogram_bins, histtype='step', color=VERSION_COLOURS[flat_version], lw=1.5,
                  label=f"{flat_version}: median {np.median(star_rms) * 1000:.1f} mmag")
    axis.set_xlabel('per-star rms around m_i + ZP_j [mag]')
    axis.legend()
    axis.set_title('repeatability (no star-flat model)')

    # Panel 6: radial coverage of the measurements, with the ghost profile for orientation
    axis = axes[1, 2]
    axis.hist(measurements['r'], bins=np.arange(0, 19.5, 0.25), color='grey')
    axis.set_xlabel('r from ghost centre [arcmin]')
    axis.set_ylabel('measurements')
    ghost_axis = axis.twinx()
    ghost_axis.plot(radius_grid, ghost_model.ghost_fraction_at_radius(radius_grid) * 100, 'k-')
    ghost_axis.set_ylabel('G [%]')
    axis.set_title(f"coverage: {results['nstars']} stars, {results['nmeas']} measurements, {results['nframes']} frames")

    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, f'starflat_{camera}_{filter_name}.png'), dpi=90)
    plt.close()


def plot_deviations(filter_name, camera, measurements, calibrations_without_star_flat, ghost_model, args, image_shape,
                    radial_bin_edges):
    """ Photometric deviation of every measurement, dm = m_ij - m_i - ZP_j (no star-flat model), for A and B.

    The prediction 1.086 * G is passed through the identical m_i + ZP_j fit, since the per-star mean and per-frame
    zero point absorb part of the ghost signature; measured and predicted deviations are then directly comparable.
    """
    predicted_deviation = self_calibrate(measurements, MAG_PER_FRACTIONAL_FLUX * measurements['G'],
                                         measurements['mag_error_A0'], image_shape, radial_bin_edges,
                                         args.normalisation_zone, fit_star_flat=False, clip_sigma=1e9)['residual']
    colour_limit_mmag = 1000 * max(np.percentile(np.abs(predicted_deviation), 99), 0.005)
    cells_per_axis = 32

    # Ghost model sampled every 32 pixels, for contour overlays
    contour_step = 32
    contour_y, contour_x = np.mgrid[:ghost_model.one_plus_ghost_image.shape[0]:contour_step,
                                    :ghost_model.one_plus_ghost_image.shape[1]:contour_step]
    contour_ghost_percent = (ghost_model.one_plus_ghost_image[::contour_step, ::contour_step] - 1) * 100
    radial_bin_centres = 0.5 * (radial_bin_edges[1:] + radial_bin_edges[:-1])

    fig, axes = plt.subplots(2, 4, figsize=(24, 11.5))
    for row, flat_version in enumerate(FLAT_VERSIONS):
        calibration = calibrations_without_star_flat[photometry_key(flat_version, 0)]
        is_used = calibration['is_used']
        deviation_mmag = calibration['residual'][is_used] * 1000
        predicted_mmag = predicted_deviation[is_used] * 1000
        x_pixel = measurements['x'][is_used]
        y_pixel = measurements['y'][is_used]
        radius_arcmin = measurements['r'][is_used]
        version_colour = VERSION_COLOURS[flat_version]

        # Panel 1: individual measurements on the detector, ghost model contours
        axis = axes[row, 0]
        drawing_order = np.argsort(np.abs(deviation_mmag))  # largest deviations drawn on top
        image = axis.scatter(x_pixel[drawing_order], y_pixel[drawing_order], c=deviation_mmag[drawing_order], s=6,
                             cmap='RdBu_r', vmin=-colour_limit_mmag, vmax=colour_limit_mmag, lw=0)
        contours = axis.contour(contour_x, contour_y, contour_ghost_percent, levels=np.linspace(-3, 3, 13), colors='k',
                                linewidths=0.5, alpha=0.5)
        axis.clabel(contours, fontsize=6, fmt='%.1f%%')
        axis.plot(ghost_model.centre_x, ghost_model.centre_y, 'k+', ms=12)
        axis.set_xlim(0, image_shape[1])
        axis.set_ylim(0, image_shape[0])
        axis.set_aspect('equal')
        axis.set_title(f'{flat_version}: Δm of {is_used.sum()} measurements [mmag]; contours: ghost G [%]')
        plt.colorbar(image, ax=axis, fraction=0.046)

        # Panel 2: mean deviation in detector cells
        axis = axes[row, 1]
        deviation_map = binned_detector_map(x_pixel, y_pixel, deviation_mmag, image_shape, cells_per_axis, min_count=3)
        image = axis.imshow(deviation_map, origin='lower', cmap='RdBu_r', vmin=-colour_limit_mmag, vmax=colour_limit_mmag,
                            extent=(0, image_shape[1], 0, image_shape[0]))
        axis.contour(contour_x, contour_y, contour_ghost_percent, levels=[0], colors='k', linewidths=0.8, linestyles='--')
        axis.set_title(f'{flat_version}: mean Δm in {image_shape[1] // cells_per_axis} px cells [mmag]; dashed: G = 0')
        plt.colorbar(image, ax=axis, fraction=0.046)

        # Panel 3: deviation vs radius from the ghost centre, binned medians and prediction
        axis = axes[row, 2]
        axis.scatter(radius_arcmin, deviation_mmag, s=2, color='grey', alpha=0.4, lw=0, label='measurements')
        radial_bin_index = np.digitize(radius_arcmin, radial_bin_edges) - 1
        binned_median = np.full(len(radial_bin_centres), np.nan)
        binned_median_error = np.full(len(radial_bin_centres), np.nan)
        binned_prediction = np.full(len(radial_bin_centres), np.nan)
        for bin_index in range(len(radial_bin_centres)):
            is_in_bin = radial_bin_index == bin_index
            if is_in_bin.sum() >= 10:
                binned_median[bin_index] = np.median(deviation_mmag[is_in_bin])
                # standard error of the median for normally distributed data: 1.2533 * sigma / sqrt(N)
                binned_median_error[bin_index] = 1.2533 * np.std(deviation_mmag[is_in_bin]) / np.sqrt(is_in_bin.sum())
                binned_prediction[bin_index] = np.median(predicted_mmag[is_in_bin])
        axis.errorbar(radial_bin_centres, binned_median, binned_median_error, fmt='o', ms=4, color=version_colour,
                      label='binned median')
        axis.plot(radial_bin_centres, binned_prediction, 'k-', label='predicted from ghost model (A)')
        axis.axhline(0, color='k', lw=0.5)
        axis.axvspan(*args.normalisation_zone, color='grey', alpha=0.12)
        axis.set_ylim(-2.5 * colour_limit_mmag, 2.5 * colour_limit_mmag)
        axis.set_xlabel('r from ghost centre [arcmin]')
        axis.set_ylabel('Δm [mmag]  (+ = too faint)')
        axis.set_title(f'{flat_version}: deviation vs radius')
        axis.legend(fontsize=8, loc='lower right')
        axis.grid(alpha=0.3)

        # Panel 4: measured vs predicted deviation, with the weighted slope through the origin
        axis = axes[row, 3]
        weights = 1 / measurements[f'mag_error_{photometry_key(flat_version, 0)}'][is_used] ** 2
        slope = np.sum(weights * predicted_mmag * deviation_mmag) / np.sum(weights * predicted_mmag ** 2)
        prediction_bin_edges = np.linspace(-colour_limit_mmag, colour_limit_mmag, 17)
        prediction_bin_index = np.digitize(predicted_mmag, prediction_bin_edges) - 1
        binned_measured = [np.median(deviation_mmag[prediction_bin_index == bin_index])
                           if (prediction_bin_index == bin_index).sum() >= 10 else np.nan
                           for bin_index in range(len(prediction_bin_edges) - 1)]
        axis.scatter(predicted_mmag, deviation_mmag, s=2, color='grey', alpha=0.4, lw=0)
        axis.plot(0.5 * (prediction_bin_edges[1:] + prediction_bin_edges[:-1]), binned_measured, 'o', color=version_colour,
                  label='binned median')
        axis.plot([-colour_limit_mmag, colour_limit_mmag], [-colour_limit_mmag, colour_limit_mmag], 'k--', lw=1,
                  label='1:1 (ghost uncorrected)')
        axis.plot([-colour_limit_mmag, colour_limit_mmag], [-colour_limit_mmag * slope, colour_limit_mmag * slope], '-',
                  color=version_colour, label=f'fit slope {slope:.3f}')
        axis.axhline(0, color='k', lw=0.5)
        axis.set_xlim(-colour_limit_mmag, colour_limit_mmag)
        axis.set_ylim(-2.5 * colour_limit_mmag, 2.5 * colour_limit_mmag)
        axis.set_xlabel('predicted Δm from ghost model [mmag]')
        axis.set_ylabel('measured Δm [mmag]')
        axis.set_title(f'{flat_version}: measured vs predicted')
        axis.legend(fontsize=8, loc='upper left')
        axis.grid(alpha=0.3)

    fig.suptitle(f'{camera} {filter_name}: photometric deviations  Δm = m_ij − m_i − ZP_j   '
                 f'(A: ghosted flat, B: ghost-corrected flat)')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, f'starflat_{camera}_{filter_name}_deviations.png'), dpi=80)
    plt.close()


# ----------------------------------------------------------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------------------------------------------------------
def parseCommandLine():
    parser = argparse.ArgumentParser(description='Star-flat photometry test of pupil-ghost impact and mitigation',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('files', nargs='+', help='BANZAI e91 frames of one dithered field')
    parser.add_argument('--ghost-model-dir', required=True, help='output directory of model_pupilghost.py --fullres')
    parser.add_argument('--outdir', default='plots/starflat')
    parser.add_argument('--template', default=None,
                        help='ghost template filter; set to e.g. Y for a control filter (rp) without own ghost model')
    parser.add_argument('--aperture-radii-fwhm', type=float, nargs='+', default=[1.5, 2.5],
                        help='aperture radii in units of the frame FWHM; the first one is the main aperture')
    parser.add_argument('--annulus-fwhm', type=float, nargs=2, default=[4, 6],
                        help='inner and outer radius of the local background annulus [FWHM]')
    parser.add_argument('--background-mesh-size', type=int, default=64,
                        help='background mesh size [px]; must be much smaller than the ghost structure')
    parser.add_argument('--detection-threshold', type=float, default=10, help='source detection threshold [sigma]')
    parser.add_argument('--max-peak-fraction', type=float, default=0.5,
                        help='reject stars whose peak exceeds this fraction of SATURATE (non-linearity)')
    parser.add_argument('--min-snr', type=float, default=100, help='minimum signal-to-noise of a measurement')
    parser.add_argument('--min-observations', type=int, default=5, help='minimum number of frames a star must be in')
    parser.add_argument('--error-floor', type=float, default=0.002, help='magnitude error floor added in quadrature')
    parser.add_argument('--match-radius', type=float, default=1.5, help='cross-match radius between frames [arcsec]')
    parser.add_argument('--radial-bin-width', type=float, default=0.5, help='radial bin width of the star flat [arcmin]')
    parser.add_argument('--normalisation-zone', type=float, nargs=2, default=[12.5, 16.5],
                        help='radial range [arcmin] assumed ghost-free, where the star flat is defined as 0')
    parser.add_argument('--max-template-radius', type=float, default=18.5,
                        help='outer radius used in the ghost template fit [arcmin]')
    parser.add_argument('--colour-catalogue', default=None,
                        help='optional CSV with columns ra,dec,color for the colour split test')
    parser.add_argument('--loglevel', default='INFO', choices=['DEBUG', 'INFO', 'WARN'])
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.loglevel), format='%(asctime)s %(levelname)s %(message)s')
    return args


def print_results(filter_results):
    """ Human-readable summary of one filter's results on stdout."""
    print(f"\n=== {filter_results['camera']} {filter_results['filter']} (template {filter_results['template']}): {filter_results['nstars']} stars, "
          f"{filter_results['nmeas']} measurements, {filter_results['nframes']} frames, "
          f"FWHM {filter_results['fwhm_arcsec_median']:.2f}\"")
    for flat_version in FLAT_VERSIONS:
        version_results = filter_results[flat_version]
        print(f"  {flat_version}: k = {version_results['k']:+.3f} ± {version_results['k_err']:.3f} "
              f"(S/N {version_results['k_snr']:.1f})  radial rms {version_results['radial_rms_pct']:.3f}%  "
              f"paired {version_results['paired_dmag_mmag']:+.1f} ± {version_results['paired_err_mmag']:.1f} mmag "
              f"(pred. A {version_results['paired_predicted_A_mmag']:+.1f})  "
              f"star rms {version_results['star_rms_median_mmag']:.1f} mmag")
    fitted_centre = filter_results['ghost_centre_fit_A']
    model_centre = filter_results['ghost_centre_model']
    print(f"  ghost centre fit (A): {fitted_centre[0]:.0f},{fitted_centre[1]:.0f}  model "
          f"{model_centre[0]:.0f},{model_centre[1]:.0f}  offset {filter_results['ghost_centre_offset_arcmin']:.2f}'")
    for criterion, has_passed in filter_results['criteria'].items():
        print(f"  [{'PASS' if has_passed else 'FAIL'}] {criterion}")


def main():
    args = parseCommandLine()
    os.makedirs(args.outdir, exist_ok=True)

    # Full 4k frames contain many sources; raise sep's internal limits
    sep.set_extract_pixstack(10_000_000)
    sep.set_sub_object_limit(4096)

    # Analyse each filter separately
    files_by_filter = {}
    for filename in args.files:
        filter_name = fits.getval(filename, 'FILTER', extname='SCI')
        files_by_filter.setdefault(filter_name, []).append(filename)

    all_results = {}
    for filter_name, files in files_by_filter.items():
        filter_results = analyse_filter(filter_name, sorted(files), args)
        if filter_results is None:
            continue
        all_results[filter_name] = filter_results
        print_results(filter_results)

    cameras = '_'.join(sorted({filter_results['camera'] for filter_results in all_results.values()})) or 'camera'
    filters = '_'.join(all_results) or 'none'
    with open(os.path.join(args.outdir, f'starflat_{cameras}_{filters}_summary.json'), 'w') as summary_file:
        json.dump(all_results, summary_file, indent=1, default=float)


if __name__ == '__main__':
    main()
