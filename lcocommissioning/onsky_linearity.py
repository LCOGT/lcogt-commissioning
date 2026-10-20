"""
On-sky linearity test of a camera from a series of BANZAI-processed images of a star field with increasing exposure times.

Principle
---------
After dividing by the exposure time, every star's count rate must be the same in every image. A non-linear detector
breaks this in a way that depends on the signal level of the star in that image: bright stars in long exposures
reach high pixel levels, faint stars in short exposures stay low.

Model:  count rate_ij = R_i * Z_j * (1 + N(peak_ij))
    R_i   true count rate of star i
    Z_j   throughput of image j: transparency, airmass, and a possible exposure time (shutter) error
    N     non-linearity as a function of the peak pixel level of star i in image j

R_i and Z_j are determined only from measurements whose peak pixel lies in a reference range (default 1000 - 20000 e-),
which defines N = 0 there. Every measurement's deviation  rate_ij / (R_i * Z_j) - 1  then traces N.

Photometry
----------
BANZAI's CAT extension drops saturated stars and its fluxes can be biased, so stars are detected independently on the
deepest image with DAOFIND (photutils DAOStarFinder), and measured with fixed circular apertures (in pixels) at the same
sky positions in every image (forced photometry with photutils, re-centred per image). The local sky is the sigma-clipped
median of an annulus around each star. DAOFIND does not return saturated stars with flat-topped cores; these are added
from the clusters of saturated pixels in the BANZAI bad pixel mask.

Usage:
    onskylinearity --outdir plots/linearity /path/to/*-e91.fits.fz
"""
import argparse
import json
import logging
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.stats import SigmaClip
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from photutils.aperture import ApertureStats, CircularAnnulus, CircularAperture, aperture_photometry
from photutils.background import Background2D, MedianBackground
from photutils.centroids import centroid_com, centroid_sources
from photutils.detection import DAOStarFinder
from scipy import ndimage
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree

from lcocommissioning.common.astrometric_registration import frame_has_valid_astrometry, register_frame

_log = logging.getLogger(__name__)

# BANZAI bad pixel mask bits. Bit 8 marks cosmic rays: in the ep60 test data it flags small (~6 pixel) clusters spread
# evenly over the detector, whose number grows with exposure time.
BPM_BAD_PIXEL = 1
BPM_SATURATED = 2
BPM_COSMIC_RAY = 8

# Half width of the box [pixels] in which the peak pixel value of a star is searched
PEAK_SEARCH_HALF_WIDTH = 3

# Relative error floor added in quadrature to every count rate (systematics beyond photon and read noise)
RELATIVE_ERROR_FLOOR = 0.002

# Signal levels [e-] at which the fitted non-linearity is reported
REPORT_LEVELS = [100, 1000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 80000, 100000, 150000, 200000, 250000]


# ----------------------------------------------------------------------------------------------------------------------
# Reading images
# ----------------------------------------------------------------------------------------------------------------------
def read_image(filename):
    """ Read science image, bad pixel mask, error image, header and BANZAI catalogue of one BANZAI frame."""
    with fits.open(filename) as hdulist:
        header = hdulist['SCI'].header
        image = np.ascontiguousarray(hdulist['SCI'].data, dtype=np.float32)  # native byte order, float32
        if 'BPM' in hdulist:
            bad_pixel_mask = np.ascontiguousarray(hdulist['BPM'].data, dtype=np.uint8)
        else:
            bad_pixel_mask = np.zeros(image.shape, np.uint8)
        if 'ERR' in hdulist:
            error_image = np.ascontiguousarray(hdulist['ERR'].data, dtype=np.float32)
        else:
            error_image = None
        banzai_catalogue = Table(hdulist['CAT'].data) if 'CAT' in hdulist else None
    return image, bad_pixel_mask, error_image, header, banzai_catalogue


def image_metadata(filename, header):
    """ Exposure time, levels and observing conditions of one image."""
    saturation_level = header.get('SATURATE', np.inf)
    return dict(file=os.path.basename(filename), exptime=float(header['EXPTIME']), filter=header.get('FILTER'),
                object=header.get('OBJECT'), dateobs=header.get('DATE-OBS'), mjd=header.get('MJD-OBS', np.nan),
                airmass=header.get('AIRMASS', np.nan), fwhm_arcsec=header.get('L1FWHM', np.nan),
                pixscale=header.get('PIXSCALE', np.nan),
                instrument=header.get('INSTRUME'), saturate=float(saturation_level),
                maxlin=float(header.get('MAXLIN', saturation_level)))


# ----------------------------------------------------------------------------------------------------------------------
# Star list
# ----------------------------------------------------------------------------------------------------------------------
def estimate_background(image, mask, box_size=64):
    """ Smooth 2D background: sigma-clipped median in box_size x box_size pixel boxes, median-filtered over 3 x 3 boxes."""
    return Background2D(image, (box_size, box_size), filter_size=(3, 3), mask=mask, sigma_clip=SigmaClip(sigma=3.0),
                        bkg_estimator=MedianBackground())


def saturated_star_positions(image, bad_pixel_mask, min_pixels=3):
    """ Centres and peaks of clusters of saturated pixels (BANZAI BPM bit 2), i.e. of saturated stars.

    DAOFIND does not return stars with a flat-topped, saturated core; their position is taken as the centre of the
    saturated pixel cluster instead. Clusters smaller than min_pixels are ignored.
    """
    is_saturated = (bad_pixel_mask & BPM_SATURATED) > 0
    cluster_labels, number_of_clusters = ndimage.label(is_saturated)
    if number_of_clusters == 0:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    cluster_index = np.arange(1, number_of_clusters + 1)
    cluster_size = ndimage.sum(is_saturated, cluster_labels, cluster_index)
    centres = np.array(ndimage.center_of_mass(is_saturated, cluster_labels, cluster_index))
    peaks = np.array(ndimage.maximum(image, cluster_labels, cluster_index))
    is_kept = cluster_size >= min_pixels
    return centres[is_kept, 1], centres[is_kept, 0], peaks[is_kept]


def detect_sources(filename, args):
    """ Detect all sources on one (the deepest) image with DAOFIND, plus saturated stars from the bad pixel mask.

    :return: dict with pixel positions, flux, peak (above background) and detection significance (peak / background
             noise) of every source, their sky positions, the image shape and the background noise
    """
    image, bad_pixel_mask, _, header, _ = read_image(filename)
    is_bad_pixel = (bad_pixel_mask & BPM_BAD_PIXEL) > 0
    background = estimate_background(image, is_bad_pixel)
    background_subtracted = image - background.background
    background_noise = float(background.background_rms_median)
    fwhm_pixel = header['L1FWHM'] / header['PIXSCALE'] if 'L1FWHM' in header and np.isfinite(header['L1FWHM']) else 4.

    # DAOFIND: stars above detection_threshold x background noise, with the default sharpness and roundness limits
    star_finder = DAOStarFinder(threshold=args.detection_threshold * background_noise, fwhm=fwhm_pixel, exclude_border=True)
    stars = star_finder(background_subtracted, mask=is_bad_pixel)
    x_pixel = np.asarray(stars['x_centroid']) if stars is not None else np.zeros(0)
    y_pixel = np.asarray(stars['y_centroid']) if stars is not None else np.zeros(0)
    flux = np.asarray(stars['flux']) if stars is not None else np.zeros(0)
    peak = np.asarray(stars['peak']) if stars is not None else np.zeros(0)

    # Saturated stars missed by DAOFIND; their flux is estimated in a 3 FWHM radius aperture
    saturated_x, saturated_y, saturated_peak = saturated_star_positions(image, bad_pixel_mask)
    if len(saturated_x) > 0:
        if len(x_pixel) > 0:
            distance_to_daofind_star = cKDTree(np.c_[x_pixel, y_pixel]).query(np.c_[saturated_x, saturated_y])[0]
            is_new = distance_to_daofind_star > 2 * fwhm_pixel
        else:
            is_new = np.ones(len(saturated_x), bool)
        saturated_flux = aperture_photometry(background_subtracted,
                                             CircularAperture(np.c_[saturated_x, saturated_y][is_new], r=3 * fwhm_pixel))['aperture_sum']
        x_pixel = np.r_[x_pixel, saturated_x[is_new]]
        y_pixel = np.r_[y_pixel, saturated_y[is_new]]
        flux = np.r_[flux, np.asarray(saturated_flux)]
        peak = np.r_[peak, saturated_peak[is_new] - np.median(background.background)]
        _log.info(f"{os.path.basename(filename)}: {is_new.sum()} saturated stars added from the bad pixel mask")
    _log.info(f"{os.path.basename(filename)}: {len(x_pixel)} sources detected")
    return dict(x=x_pixel, y=y_pixel, flux=flux, peak=peak, detection_significance=peak / background_noise,
                sky_positions=WCS(header).pixel_to_world(x_pixel, y_pixel), image_shape=image.shape,
                background_noise=background_noise)


def nearest_neighbour_distance(sources):
    """ Distance [pixels] from every source to its nearest other source."""
    source_positions = np.c_[sources['x'], sources['y']]
    return cKDTree(source_positions).query(source_positions, k=2)[0][:, 1]


def select_isolated_stars(sources, isolation_radius, max_neighbour_flux_ratio):
    """ Sources without a relevant neighbour within isolation_radius, and far enough from the image edge.

    A neighbour is relevant if its flux exceeds max_neighbour_flux_ratio times the flux of the star. Fainter
    detections next to bright stars (DAOFIND finds features in the wings of bright stars) change the aperture flux by
    less than that ratio and leave the median sky of the annulus unchanged, so they do not disqualify a star.

    :return: sky positions and detection significance of the selected sources
    """
    image_shape = sources['image_shape']
    source_positions = np.c_[sources['x'], sources['y']]
    neighbours_within_radius = cKDTree(source_positions).query_ball_point(source_positions, isolation_radius)
    is_isolated = np.array([all(sources['flux'][neighbour] <= max_neighbour_flux_ratio * sources['flux'][star]
                                for neighbour in neighbours if neighbour != star)
                            for star, neighbours in enumerate(neighbours_within_radius)])
    edge_margin = isolation_radius + 2
    is_away_from_edge = (sources['x'] > edge_margin) & (sources['x'] < image_shape[1] - edge_margin)
    is_away_from_edge &= (sources['y'] > edge_margin) & (sources['y'] < image_shape[0] - edge_margin)
    is_selected = is_isolated & is_away_from_edge
    _log.info(f"{is_selected.sum()} of {len(is_selected)} sources isolated within {isolation_radius:g} px")
    return sources['sky_positions'][is_selected], sources['detection_significance'][is_selected]


# ----------------------------------------------------------------------------------------------------------------------
# Photometry
# ----------------------------------------------------------------------------------------------------------------------
def count_flagged_pixels_in_aperture(bad_pixel_mask, flag_bit, x_pixel, y_pixel, radius):
    """ Number of pixels with the given BPM bit set inside each circular aperture."""
    flag_image = ((bad_pixel_mask & flag_bit) > 0).astype(np.float32)
    aperture = CircularAperture(np.c_[x_pixel, y_pixel], r=radius)
    return np.asarray(aperture_photometry(flag_image, aperture, method='center')['aperture_sum'])


def aperture_photometry_with_local_sky(image, x_pixel, y_pixel, radius, annulus_radii, error_image=None, mask=None):
    """ Sky-subtracted flux in a circular aperture, with the sky per pixel from the sigma-clipped median of an annulus.

    Masked pixels are excluded from both the aperture sum and the sky; the sky is subtracted for the unmasked aperture
    area. The flux error combines the pixel errors in the aperture with the uncertainty of the sky level.

    :return: flux, flux error (NaN if no error image is given)
    """
    positions = np.c_[x_pixel, y_pixel]
    aperture = CircularAperture(positions, r=radius)
    annulus = CircularAnnulus(positions, r_in=annulus_radii[0], r_out=annulus_radii[1])

    sky_statistics = ApertureStats(image, annulus, mask=mask, sigma_clip=SigmaClip(sigma=3.0, maxiters=5))
    sky_per_pixel = np.asarray(sky_statistics.median)
    sky_pixel_noise = np.asarray(sky_statistics.std)
    sky_pixel_count = np.asarray(getattr(sky_statistics.sum_aper_area, 'value', sky_statistics.sum_aper_area))

    photometry = aperture_photometry(image, aperture, error=error_image, mask=mask)
    aperture_area = np.asarray(aperture.area_overlap(image, mask=mask))
    flux = np.asarray(photometry['aperture_sum']) - sky_per_pixel * aperture_area
    if error_image is None:
        return flux, np.full(len(flux), np.nan)
    sky_level_error = sky_pixel_noise / np.sqrt(np.maximum(sky_pixel_count, 1))
    flux_error = np.sqrt(np.asarray(photometry['aperture_sum_err']) ** 2 + (aperture_area * sky_level_error) ** 2)
    return flux, flux_error


def peak_pixel_values(image, x_pixel, y_pixel):
    """ Highest pixel value (including sky) in a small box around each star: the level the detector actually saw."""
    peaks = np.full(len(x_pixel), np.nan)
    for star_index, (x, y) in enumerate(zip(x_pixel, y_pixel)):
        column, row = int(round(x)), int(round(y))
        box = image[max(row - PEAK_SEARCH_HALF_WIDTH, 0):row + PEAK_SEARCH_HALF_WIDTH + 1,
                    max(column - PEAK_SEARCH_HALF_WIDTH, 0):column + PEAK_SEARCH_HALF_WIDTH + 1]
        if box.size > 0:
            peaks[star_index] = np.max(box)
    return peaks


def measure_sky(image, bad_pixel_mask):
    """ Sky level [e-/pixel] and its robust pixel-to-pixel noise, from every 4th unflagged pixel.

    In a sparse star field the median is dominated by sky pixels; the noise is the MAD scaled to a Gaussian sigma.
    """
    sampled_pixels = image[::4, ::4]
    is_clean = (bad_pixel_mask[::4, ::4] == 0) & np.isfinite(sampled_pixels)
    sky_pixels = sampled_pixels[is_clean]
    sky_level = float(np.median(sky_pixels))
    sky_noise = float(1.4826 * np.median(np.abs(sky_pixels - sky_level)))
    return sky_level, sky_noise


def locate_stars(image, wcs, sky_level, fwhm_pixel, star_positions, edge_margin, max_centroid_shift):
    """ Pixel positions of the stars in one image: WCS prediction, re-centred with a centre-of-mass centroid in a
    box of about 2 FWHM around the prediction.

    Stars closer than edge_margin to the image edge are dropped. If the centroid moves by more than
    max_centroid_shift, the predicted position is kept and the star is marked as not re-centred.

    :return: star index into star_positions, x, y, is_recentred
    """
    x_predicted, y_predicted = wcs.world_to_pixel(star_positions)
    is_inside = (x_predicted > edge_margin) & (x_predicted < image.shape[1] - edge_margin)
    is_inside &= (y_predicted > edge_margin) & (y_predicted < image.shape[0] - edge_margin)
    star_id = np.arange(len(star_positions))[is_inside]
    x_predicted = x_predicted[is_inside]
    y_predicted = y_predicted[is_inside]

    box_size = int(2 * np.ceil(fwhm_pixel)) + 1  # odd number of pixels
    x_centroid, y_centroid = centroid_sources(image - sky_level, x_predicted, y_predicted, box_size=box_size,
                                              centroid_func=centroid_com)
    centroid_shift = np.hypot(x_centroid - x_predicted, y_centroid - y_predicted)
    is_recentred = np.isfinite(centroid_shift) & (centroid_shift < max_centroid_shift)
    x_pixel = np.where(is_recentred, x_centroid, x_predicted)
    y_pixel = np.where(is_recentred, y_centroid, y_predicted)
    return star_id, x_pixel, y_pixel, is_recentred


def fwhm_in_pixels(metadata):
    """ Seeing FWHM [pixels] from the BANZAI header; 4 pixels if unknown."""
    return metadata['fwhm_arcsec'] / metadata['pixscale'] if np.isfinite(metadata['fwhm_arcsec']) else 4.


def measure_image(filename, star_positions, args, wcs=None):
    """ Forced aperture photometry of all stars in one image.

    :param wcs: WCS of the image (e.g. from astrometric registration); the header WCS if None

    :return: dict of per-star arrays, image metadata dict
    """
    image, bad_pixel_mask, error_image, header, banzai_catalogue = read_image(filename)
    metadata = image_metadata(filename, header)

    # Sky level of the image, then the stars' positions in this image (kept well inside the image)
    sky_level, sky_noise = measure_sky(image, bad_pixel_mask)
    metadata['sky_level'] = sky_level
    metadata['sky_noise'] = sky_noise
    metadata['sky_rate'] = sky_level / metadata['exptime']
    star_id, x_pixel, y_pixel, is_recentred = locate_stars(image, wcs if wcs is not None else WCS(header), sky_level,
                                                           fwhm_in_pixels(metadata), star_positions, args.annulus[1] + 2,
                                                           args.max_centroid_shift)

    measurements = dict(star=star_id, x=x_pixel, y=y_pixel, is_recentred=is_recentred,
                        peak=peak_pixel_values(image, x_pixel, y_pixel))

    # Fixed-radius apertures; the local sky is the sigma-clipped median of an annulus around each star. Bad pixels and cosmic rays are
    # masked, which keeps them out of the sky estimate; a measurement with such pixels inside its own aperture is
    # rejected later (n_masked_r*). Saturated pixels are deliberately not masked: their flux is what we want to see.
    if error_image is None:
        error_image = np.sqrt(np.maximum(image, 0) + header.get('RDNOISE', 0) ** 2).astype(np.float32)
    is_masked_pixel = (bad_pixel_mask & (BPM_BAD_PIXEL | BPM_COSMIC_RAY)) > 0
    for radius in args.aperture_radii:
        flux, flux_error = aperture_photometry_with_local_sky(image, x_pixel, y_pixel, radius, args.annulus,
                                                              error_image=error_image, mask=is_masked_pixel)
        measurements[f'flux_r{radius:g}'] = flux
        measurements[f'flux_error_r{radius:g}'] = flux_error
        measurements[f'n_masked_r{radius:g}'] = count_flagged_pixels_in_aperture(
            bad_pixel_mask, BPM_BAD_PIXEL | BPM_COSMIC_RAY, x_pixel, y_pixel, radius)

    # Saturated pixels (above MAXLIN) inside the largest aperture
    largest_radius = max(args.aperture_radii)
    measurements['n_saturated'] = count_flagged_pixels_in_aperture(bad_pixel_mask, BPM_SATURATED, x_pixel, y_pixel, largest_radius)

    # BANZAI catalogue flux of the same star, for comparison with our own photometry. Matched on the sky if the catalogue
    # has ra / dec, otherwise (no astrometric solution) on pixel positions within 1 pixel; the catalogue's pixel
    # convention (0- or 1-based) is taken from whichever gives more matches.
    measurements['banzai_flux'] = np.full(len(star_id), np.nan)
    if banzai_catalogue is not None and len(banzai_catalogue) > 0 and len(star_id) > 0:
        if 'ra' in banzai_catalogue.colnames and 'dec' in banzai_catalogue.colnames:
            catalogue_positions = SkyCoord(banzai_catalogue['ra'], banzai_catalogue['dec'], unit='deg')
            catalogue_index, separation, _ = star_positions[star_id].match_to_catalog_sky(catalogue_positions)
            is_matched = separation.arcsec < 1.0
        else:
            catalogue_tree = cKDTree(np.c_[banzai_catalogue['x'], banzai_catalogue['y']])
            best_match = None
            for pixel_convention_offset in (0.0, 1.0):
                distance, catalogue_index = catalogue_tree.query(np.c_[x_pixel, y_pixel] + pixel_convention_offset)
                is_matched = distance < 1.0
                if best_match is None or is_matched.sum() > best_match[1].sum():
                    best_match = (catalogue_index, is_matched)
            catalogue_index, is_matched = best_match
        measurements['banzai_flux'][is_matched] = banzai_catalogue['flux'][catalogue_index[is_matched]]

    _log.info(f"{metadata['file']}: EXPTIME {metadata['exptime']:.1f} s, airmass {metadata['airmass']:.3f}, "
              f"FWHM {metadata['fwhm_arcsec']:.2f}\", {len(star_id)} stars measured")
    return measurements, metadata


def build_measurement_table(per_image_measurements, image_metadata_list, args):
    """ One row per (star, image) measurement, with count rates and quality flags."""
    columns = {}
    for image_index, (measurements, metadata) in enumerate(zip(per_image_measurements, image_metadata_list)):
        measurements = dict(measurements)
        measurements['image'] = np.full(len(measurements['star']), image_index)
        measurements['exptime'] = np.full(len(measurements['star']), metadata['exptime'])
        for name, values in measurements.items():
            columns.setdefault(name, []).append(values)
    table = {name: np.concatenate(values) for name, values in columns.items()}

    # Count rates [e-/s] per aperture, with a relative error floor
    for radius in args.aperture_radii:
        flux = table[f'flux_r{radius:g}']
        flux_error = table[f'flux_error_r{radius:g}']
        table[f'rate_r{radius:g}'] = flux / table['exptime']
        table[f'rate_error_r{radius:g}'] = np.sqrt((flux_error / table['exptime']) ** 2
                                                   + (RELATIVE_ERROR_FLOOR * table[f'rate_r{radius:g}']) ** 2)

    # Quality classes, per aperture: usable (clean, enough signal) and saturated (above MAXLIN; shown, never fitted)
    table['is_saturated'] = table['n_saturated'] > 0
    for radius in args.aperture_radii:
        flux = table[f'flux_r{radius:g}']
        with np.errstate(divide='ignore', invalid='ignore'):
            signal_to_noise = flux / table[f'flux_error_r{radius:g}']
        table[f'is_usable_r{radius:g}'] = (np.isfinite(flux) & (flux > 0) & (signal_to_noise > args.min_snr)
                                           & (table[f'n_masked_r{radius:g}'] == 0) & table['is_recentred'])
    table['is_usable'] = table[f'is_usable_r{args.main_aperture:g}']
    return table


# ----------------------------------------------------------------------------------------------------------------------
# Aperture optimisation
# ----------------------------------------------------------------------------------------------------------------------
# Radii [pixels] of the curve of growth; the curve is normalised at the largest radius
GROWTH_CURVE_RADII = np.arange(1, 41, dtype=float)
# Sky annulus [pixels] for the curve of growth, far outside the largest radius
GROWTH_CURVE_SKY_ANNULUS = (45., 60.)
# Half size [pixels] of the PSF model grid used for the neighbour contamination
PSF_GRID_HALF_SIZE = 128
# Minimum signal-to-noise of a star's total flux for its curve of growth; at least this many of the brightest stars
# are used if fewer qualify
GROWTH_CURVE_MIN_SNR = 300
GROWTH_CURVE_MIN_STARS = 20


def measure_growth_curves(filename, isolated_star_positions, args, wcs=None):
    """ Curve of growth (fraction of light inside radius r, normalised at the largest radius) of bright, unsaturated,
    isolated stars in one image.

    :return: dict with exposure time, FWHM, sky noise and the curves of growth (stars x radii)
    """
    image, bad_pixel_mask, error_image, header, _ = read_image(filename)
    metadata = image_metadata(filename, header)
    sky_level, sky_noise = measure_sky(image, bad_pixel_mask)
    if error_image is None:
        error_image = np.sqrt(np.maximum(image, 0) + header.get('RDNOISE', 0) ** 2).astype(np.float32)
    outer_radius = GROWTH_CURVE_SKY_ANNULUS[1]
    star_id, x_pixel, y_pixel, is_recentred = locate_stars(image, wcs if wcs is not None else WCS(header), sky_level,
                                                           fwhm_in_pixels(metadata), isolated_star_positions,
                                                           outer_radius + 2, args.max_centroid_shift)

    # Bright enough for a precise curve of growth, but below the linearity limit; no flagged pixel out to the sky annulus
    peak = peak_pixel_values(image, x_pixel, y_pixel)
    is_flagged = count_flagged_pixels_in_aperture(bad_pixel_mask, BPM_BAD_PIXEL | BPM_COSMIC_RAY | BPM_SATURATED,
                                                  x_pixel, y_pixel, outer_radius) > 0
    is_selected = is_recentred & ~is_flagged & (peak - sky_level > 50 * sky_noise) & (peak < 0.5 * metadata['maxlin'])
    x_pixel = x_pixel[is_selected]
    y_pixel = y_pixel[is_selected]

    if len(x_pixel) == 0:
        _log.warning(f"{metadata['file']}: no stars for the curve of growth")
        return dict(exptime=metadata['exptime'], fwhm_pixel=fwhm_in_pixels(metadata), sky_noise=sky_noise,
                    growth_curves=np.zeros((0, len(GROWTH_CURVE_RADII))), total_flux=np.zeros(0))

    # Sky per pixel once per star (sigma-clipped median of the far annulus), then the flux inside every radius
    positions = np.c_[x_pixel, y_pixel]
    sky_statistics = ApertureStats(image, CircularAnnulus(positions, *GROWTH_CURVE_SKY_ANNULUS),
                                   sigma_clip=SigmaClip(sigma=3.0, maxiters=5))
    sky_per_pixel = np.asarray(sky_statistics.median)
    apertures = [CircularAperture(positions, r=radius) for radius in GROWTH_CURVE_RADII]
    photometry = aperture_photometry(image, apertures, error=error_image)
    flux_inside_radius = np.zeros((len(x_pixel), len(GROWTH_CURVE_RADII)))
    for radius_index, aperture in enumerate(apertures):
        flux_inside_radius[:, radius_index] = (np.asarray(photometry[f'aperture_sum_{radius_index}'])
                                               - sky_per_pixel * np.asarray(aperture.area_overlap(image)))
    total_flux = flux_inside_radius[:, -1]
    total_flux_error = np.asarray(photometry[f'aperture_sum_err_{len(GROWTH_CURVE_RADII) - 1}'])

    # Only precise curves of growth: the noise of faint stars would otherwise look like star-to-star PSF scatter
    with np.errstate(divide='ignore', invalid='ignore'):
        total_flux_snr = np.where(total_flux > 0, total_flux / total_flux_error, 0)
    is_valid = total_flux_snr >= GROWTH_CURVE_MIN_SNR
    if is_valid.sum() < GROWTH_CURVE_MIN_STARS:
        brightest = np.argsort(total_flux_snr)[::-1][:GROWTH_CURVE_MIN_STARS]
        is_valid = np.zeros(len(total_flux), bool)
        is_valid[brightest[total_flux_snr[brightest] > 0]] = True
    growth_curves = flux_inside_radius[is_valid] / total_flux[is_valid, None]
    _log.info(f"{metadata['file']}: curve of growth from {is_valid.sum()} stars")
    return dict(exptime=metadata['exptime'], fwhm_pixel=fwhm_in_pixels(metadata), sky_noise=sky_noise,
                growth_curves=growth_curves, total_flux=total_flux[is_valid])


def psf_from_growth_curve(encircled_energy, radii, grid_half_size=PSF_GRID_HALF_SIZE):
    """ Circular PSF image (unit sum) built from a curve of growth:  PSF(rho) = dEE/drho / (2 pi rho)."""
    radii_with_zero = np.r_[0., radii]
    encircled_energy_with_zero = np.r_[0., encircled_energy]
    y_grid, x_grid = np.mgrid[-grid_half_size:grid_half_size + 1, -grid_half_size:grid_half_size + 1]
    rho = np.hypot(x_grid, y_grid)
    radial_derivative = np.gradient(encircled_energy_with_zero, radii_with_zero)
    psf = np.interp(np.maximum(rho, 0.5), radii_with_zero, radial_derivative, right=0.) / (2 * np.pi * np.maximum(rho, 0.5))
    psf[rho > radii[-1]] = 0
    return psf / psf.sum()


def neighbour_light_fraction_table(psf, aperture_radii):
    """ Fraction of a star's light that falls into an aperture of radius r centred at distance d from the star.

    :return: table [aperture radius index, distance in integer pixels]
    """
    grid_half_size = psf.shape[0] // 2
    y_grid, x_grid = np.mgrid[-grid_half_size:grid_half_size + 1, -grid_half_size:grid_half_size + 1]
    rho = np.hypot(x_grid, y_grid)
    table = np.zeros((len(aperture_radii), grid_half_size))
    for radius_index, radius in enumerate(aperture_radii):
        aperture_disk = (rho <= radius).astype(float)
        light_in_shifted_aperture = fftconvolve(psf, aperture_disk, mode='same')
        table[radius_index] = light_in_shifted_aperture[grid_half_size, grid_half_size:grid_half_size * 2]
    return np.clip(table, 0, None)


def neighbour_contamination(sources, light_fraction_table, min_peak):
    """ For every star brighter than min_peak: light from all other detected sources that falls into its aperture,
    relative to the star's own flux, per aperture radius.

    :return: array [star, aperture radius index]
    """
    max_distance = light_fraction_table.shape[1] - 1
    source_positions = np.c_[sources['x'], sources['y']]
    tree = cKDTree(source_positions)
    candidate_index = np.nonzero((sources['peak'] > min_peak) & (sources['flux'] > 0))[0]
    contamination = np.zeros((len(candidate_index), light_fraction_table.shape[0]))
    for row, star in enumerate(candidate_index):
        neighbours = [neighbour for neighbour in tree.query_ball_point(source_positions[star], max_distance) if neighbour != star]
        if not neighbours:
            continue
        distance = np.hypot(*(source_positions[neighbours] - source_positions[star]).T)
        flux_ratio = np.clip(sources['flux'][neighbours], 0, None) / sources['flux'][star]
        for radius_index in range(light_fraction_table.shape[0]):
            light_fraction = np.interp(distance, np.arange(max_distance + 1), light_fraction_table[radius_index])
            contamination[row, radius_index] = np.sum(flux_ratio * light_fraction)
    return contamination


def recommend_annulus(encircled_energy, radii, aperture_radius, wing_bias_tolerance, area_factor):
    """ Sky annulus for a given aperture: the smallest inner radius where the star's own wings bias the measured flux by
    less than wing_bias_tolerance, and a width that gives an annulus area of area_factor times the aperture area.

    Wing bias = (star light per pixel in the annulus) * (aperture area) / (star light in the aperture). Light beyond the
    largest growth curve radius is taken as zero.

    :return: inner radius, outer radius, inner radius grid and the wing bias on that grid
    """
    def encircled_energy_at(radius):
        return np.interp(radius, radii, encircled_energy, right=1.)

    aperture_area = np.pi * aperture_radius ** 2
    inner_radii = np.arange(np.ceil(aperture_radius) + 2, radii[-1] + 15, 1.)
    outer_radii = np.sqrt(inner_radii ** 2 + area_factor * aperture_radius ** 2)
    annulus_area = np.pi * (outer_radii ** 2 - inner_radii ** 2)
    wing_bias = ((encircled_energy_at(outer_radii) - encircled_energy_at(inner_radii)) / annulus_area * aperture_area
                 / encircled_energy_at(aperture_radius))
    is_acceptable = wing_bias < wing_bias_tolerance
    best = np.argmax(is_acceptable) if is_acceptable.any() else len(inner_radii) - 1
    return float(inner_radii[best]), float(np.ceil(outer_radii[best])), inner_radii, wing_bias


def analyse_apertures(files, sources, args, image_wcs):
    """ Optimum aperture radius and sky annulus: the balance between collecting most of the light and collecting light
    from neighbouring stars.

    For every candidate aperture radius r, the expected systematic error of a star is the quadrature sum of
      missed light   c * (1 - EE(r)): the light outside the aperture times its relative variability c. c is measured as
                     the frame-to-frame spread of the missed light divided by the missed light (median between 1.5 FWHM
                     and half the normalisation radius), i.e. how much the PSF wings change between images.
      neighbours     mean light of all detected neighbours inside the aperture, relative to the star's flux (the mean
                     includes the rare bright neighbours)
    The optimum radius minimises this total, between 1 FWHM and 75% of the normalisation radius (closer to the
    normalisation radius the curve of growth is 1 by construction). Photon/sky noise and the star-to-star scatter of the
    curve of growth are reported for information: noise averages down over many stars and does not set the choice.
    The sky annulus is then chosen for the optimum radius (see recommend_annulus).
    """
    isolated_star_positions, _ = select_isolated_stars(sources, GROWTH_CURVE_SKY_ANNULUS[1], args.isolation_flux_ratio)
    per_image = [measure_growth_curves(filename, isolated_star_positions, args, image_wcs[filename]) for filename in files]
    radii = GROWTH_CURVE_RADII

    # Light loss: median curve of growth per image, and the star-to-star scatter within each image
    median_curve_per_image = np.array([np.median(image['growth_curves'], axis=0) for image in per_image])
    star_scatter_per_image = np.array([1.4826 * np.median(np.abs(image['growth_curves'] - np.median(image['growth_curves'], axis=0)),
                                                          axis=0) for image in per_image])
    # Median over images, forced to be non-decreasing (light inside r cannot drop with r; residual noise can make it)
    median_encircled_energy = np.maximum.accumulate(np.median(median_curve_per_image, axis=0))
    light_loss_scatter = np.median(star_scatter_per_image, axis=0)
    frame_to_frame_spread = np.ptp(median_curve_per_image, axis=0)

    # Neighbour contamination with a PSF built from the median curve of growth
    psf = psf_from_growth_curve(median_encircled_energy, radii)
    light_fraction_table = neighbour_light_fraction_table(psf, radii)
    deepest_image = per_image[-1]
    contamination = neighbour_contamination(sources, light_fraction_table, min_peak=20 * sources['background_noise'])
    median_contamination = np.median(contamination, axis=0)
    contamination_90th_percentile = np.percentile(contamination, 90, axis=0)
    fraction_contaminated = np.mean(contamination > args.contamination_tolerance, axis=0)

    # Photon and sky noise of a typical star: peak at the logarithmic middle of the reference peak range, total flux
    # from the PSF (fraction of the light in the central pixel)
    typical_peak = float(np.sqrt(args.reference_peak_range[0] * args.reference_peak_range[1]))
    central_pixel_fraction = psf[psf.shape[0] // 2, psf.shape[1] // 2]
    typical_flux = typical_peak / central_pixel_fraction
    collected_flux = typical_flux * median_encircled_energy
    noise = np.sqrt(collected_flux + np.pi * radii ** 2 * deepest_image['sky_noise'] ** 2) / collected_flux

    # Missed light and its variability between images
    median_fwhm = np.median([image['fwhm_pixel'] for image in per_image])
    missed_light = 1 - median_encircled_energy
    is_variability_range = (radii >= 1.5 * median_fwhm) & (radii <= 0.5 * radii[-1]) & (missed_light > 0)
    missed_light_variability = float(np.clip(np.median(frame_to_frame_spread[is_variability_range]
                                                       / missed_light[is_variability_range]), 0.1, 1.0))
    missed_light_uncertainty = missed_light_variability * missed_light
    mean_contamination = np.mean(contamination, axis=0)

    total_error = np.sqrt(missed_light_uncertainty ** 2 + mean_contamination ** 2)
    is_candidate = (radii >= median_fwhm) & (radii <= 0.75 * radii[-1])
    optimum_index = np.nonzero(is_candidate)[0][np.argmin(total_error[is_candidate])]
    optimum_radius = float(radii[optimum_index])
    radius_99_percent = float(radii[np.argmax(median_encircled_energy >= 0.99)]) if (median_encircled_energy >= 0.99).any() else np.nan

    annulus_inner, annulus_outer, annulus_inner_grid, wing_bias = recommend_annulus(
        median_encircled_energy, radii, optimum_radius, args.annulus_wing_tolerance, args.annulus_area_factor)
    isolated_fraction_by_radius = np.array([np.mean(nearest_neighbour_distance(sources) > radius)
                                            for radius in np.arange(5, 81, 1.)])

    return dict(radii=radii, median_curve_per_image=median_curve_per_image, exposure_times=[image['exptime'] for image in per_image],
                growth_curves=[image['growth_curves'] for image in per_image], median_encircled_energy=median_encircled_energy,
                typical_peak=typical_peak,
                light_loss_scatter=light_loss_scatter, frame_to_frame_spread=frame_to_frame_spread,
                missed_light_variability=missed_light_variability, missed_light_uncertainty=missed_light_uncertainty,
                mean_contamination=mean_contamination,
                median_contamination=median_contamination, contamination_90th_percentile=contamination_90th_percentile,
                fraction_contaminated=fraction_contaminated, noise=noise, total_error=total_error, typical_flux=typical_flux,
                optimum_radius=optimum_radius, radius_99_percent=radius_99_percent, annulus_inner=annulus_inner,
                annulus_outer=annulus_outer, annulus_inner_grid=annulus_inner_grid, wing_bias=wing_bias,
                isolation_radius_grid=np.arange(5, 81, 1.), isolated_fraction_by_radius=isolated_fraction_by_radius,
                number_of_growth_curve_stars=[len(image['growth_curves']) for image in per_image])


def aperture_summary(aperture_analysis):
    """ JSON-friendly key numbers of the aperture analysis."""
    radii = aperture_analysis['radii']
    optimum_index = int(np.nonzero(radii == aperture_analysis['optimum_radius'])[0][0])
    return dict(optimum_radius_px=aperture_analysis['optimum_radius'],
                encircled_energy_at_optimum=float(aperture_analysis['median_encircled_energy'][optimum_index]),
                missed_light_variability=aperture_analysis['missed_light_variability'],
                missed_light_uncertainty_at_optimum=float(aperture_analysis['missed_light_uncertainty'][optimum_index]),
                mean_contamination_at_optimum=float(aperture_analysis['mean_contamination'][optimum_index]),
                light_loss_scatter_at_optimum=float(aperture_analysis['light_loss_scatter'][optimum_index]),
                frame_to_frame_spread_at_optimum=float(aperture_analysis['frame_to_frame_spread'][optimum_index]),
                median_contamination_at_optimum=float(aperture_analysis['median_contamination'][optimum_index]),
                noise_at_optimum=float(aperture_analysis['noise'][optimum_index]),
                total_error_at_optimum=float(aperture_analysis['total_error'][optimum_index]),
                radius_99_percent_px=aperture_analysis['radius_99_percent'],
                annulus_px=[aperture_analysis['annulus_inner'], aperture_analysis['annulus_outer']],
                growth_curve_stars_per_image=aperture_analysis['number_of_growth_curve_stars'],
                typical_star_peak=aperture_analysis['typical_peak'], typical_star_flux=aperture_analysis['typical_flux'])


# ----------------------------------------------------------------------------------------------------------------------
# Linearity analysis
# ----------------------------------------------------------------------------------------------------------------------
def weighted_mean_per_group(group_index, values, weights, number_of_groups):
    """ Weighted mean of values within each group (star or image); NaN for groups without weight."""
    # Zero-weight entries may hold NaN (e.g. from an excluded image); keep them out of the sum, since 0 * NaN = NaN
    values_with_weight = np.where(weights > 0, values, 0.)
    weighted_sum = np.bincount(group_index, weights * values_with_weight, number_of_groups)
    sum_of_weights = np.bincount(group_index, weights, number_of_groups)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(sum_of_weights > 0, weighted_sum / sum_of_weights, np.nan)


def fit_star_rates_and_image_throughput(table, radius, reference_peak_range, extra_selection=None, max_iterations=30,
                                        clip_sigma=5.):
    """ Fit  ln(rate_ij) = ln(R_i) + ln(Z_j)  using only reference measurements (peak inside the reference range).

    Alternating weighted least squares with robust sigma clipping. The throughputs Z_j are normalised to a weighted
    geometric mean of 1. Stars without any reference measurement get R_i = NaN. extra_selection (boolean per
    measurement) restricts the reference measurements further, e.g. to a high-quality subset.

    :return: dict with star rate R_i, image throughput Z_j, its relative uncertainty, and the mask of reference
             measurements used
    """
    rate = table[f'rate_r{radius:g}']
    rate_error = table[f'rate_error_r{radius:g}']
    star_index = table['star']
    image_index = table['image']
    number_of_stars = star_index.max() + 1
    number_of_images = image_index.max() + 1

    with np.errstate(divide='ignore', invalid='ignore'):
        log_rate = np.log(rate)
        log_rate_error = rate_error / rate
    is_reference = (table[f'is_usable_r{radius:g}'] & ~table['is_saturated'] & np.isfinite(log_rate)
                    & (table['peak'] >= reference_peak_range[0]) & (table['peak'] <= reference_peak_range[1]))
    if extra_selection is not None:
        is_reference &= extra_selection
    inverse_variance = np.where(is_reference, 1 / log_rate_error ** 2, 0.)

    is_used = is_reference.copy()
    log_throughput = np.zeros(number_of_images)
    safe_log_rate = np.where(is_reference, log_rate, 0.)
    for iteration in range(max_iterations):
        weights = inverse_variance * is_used
        log_star_rate = weighted_mean_per_group(star_index, safe_log_rate - log_throughput[image_index], weights, number_of_stars)
        log_throughput = weighted_mean_per_group(image_index, safe_log_rate - np.nan_to_num(log_star_rate[star_index]),
                                                 weights, number_of_images)
        weights_per_image = np.bincount(image_index, weights, number_of_images)
        log_throughput -= np.nansum(log_throughput * weights_per_image) / weights_per_image.sum()

        # Robust sigma clipping of the reference measurements
        normalised_residual = (log_rate - log_star_rate[star_index] - log_throughput[image_index]) / log_rate_error
        robust_scale = 1.4826 * np.nanmedian(np.abs(normalised_residual[is_used]))
        is_used_next = is_reference & (np.abs(normalised_residual) < clip_sigma * max(robust_scale, 1))
        if iteration > 3 and np.array_equal(is_used_next, is_used):
            break
        is_used = is_used_next

    # Relative uncertainty of Z_j from the weights of the reference measurements in each image
    sum_of_weights_per_image = np.bincount(image_index, inverse_variance * is_used, number_of_images)
    with np.errstate(divide='ignore'):
        relative_throughput_error = 1 / np.sqrt(sum_of_weights_per_image)

    return dict(star_rate=np.exp(log_star_rate), image_throughput=np.exp(log_throughput),
                relative_throughput_error=relative_throughput_error, is_reference_used=is_used)


def deviation_from_linearity(table, radius, normalisation):
    """ Relative deviation of every measurement from its linear expectation: rate_ij / (R_i * Z_j) - 1, and error."""
    expected_rate = normalisation['star_rate'][table['star']] * normalisation['image_throughput'][table['image']]
    deviation = table[f'rate_r{radius:g}'] / expected_rate - 1
    deviation_error = table[f'rate_error_r{radius:g}'] / expected_rate
    return deviation, deviation_error, expected_rate


def fit_nonlinearity_polynomial(peak, deviation, deviation_error, polynomial_degree, level_scale=1e4):
    """ Weighted polynomial fit  deviation = sum_k c_k * (peak / level_scale)^k, k = 0 .. polynomial_degree.

    :return: dict with coefficients, covariance (scaled by the reduced chi^2 if > 1), reduced chi^2, evaluation function
    """
    scaled_peak = peak / level_scale
    design = np.vander(scaled_peak, polynomial_degree + 1, increasing=True)
    weight = 1 / deviation_error
    coefficients, *_ = np.linalg.lstsq(design * weight[:, None], deviation * weight, rcond=None)
    chi2 = np.sum(((deviation - design @ coefficients) * weight) ** 2)
    degrees_of_freedom = max(len(deviation) - polynomial_degree - 1, 1)
    reduced_chi2 = chi2 / degrees_of_freedom
    covariance = np.linalg.inv((design * weight[:, None]).T @ (design * weight[:, None])) * max(reduced_chi2, 1)

    def evaluate(level):
        level_design = np.vander(np.atleast_1d(level) / level_scale, polynomial_degree + 1, increasing=True)
        value = level_design @ coefficients
        error = np.sqrt(np.einsum('ij,jk,ik->i', level_design, covariance, level_design))
        return value, error

    return dict(coefficients=coefficients, covariance=covariance, reduced_chi2=reduced_chi2, level_scale=level_scale,
                evaluate=evaluate)


def binned_median(x_values, y_values, bin_edges, min_count=5):
    """ Median of y in bins of x, with the standard error of the median (1.2533 * sigma / sqrt(N))."""
    bin_index = np.digitize(x_values, bin_edges) - 1
    number_of_bins = len(bin_edges) - 1
    medians = np.full(number_of_bins, np.nan)
    median_errors = np.full(number_of_bins, np.nan)
    counts = np.zeros(number_of_bins, int)
    for bin_number in range(number_of_bins):
        is_in_bin = bin_index == bin_number
        counts[bin_number] = is_in_bin.sum()
        if counts[bin_number] >= min_count:
            medians[bin_number] = np.median(y_values[is_in_bin])
            median_errors[bin_number] = 1.2533 * np.std(y_values[is_in_bin]) / np.sqrt(counts[bin_number])
    bin_centres = np.sqrt(bin_edges[1:] * bin_edges[:-1])  # geometric centres for log-spaced bins
    return bin_centres, medians, median_errors, counts


def pairwise_rate_ratios(table, radius, reference_peak_range):
    """ Per star: count rate in each image divided by its count rate in the shortest exposure.

    Each image's ratios are normalised by the median ratio of stars whose peak lies in the reference range in both
    images, which removes transparency and shutter differences. A linear detector gives 1 at all levels.
    Note that bright stars may already be non-linear in the shortest exposure.
    """
    rate = table[f'rate_r{radius:g}']
    is_good = table['is_usable']
    shortest_image = table['image'][np.argmin(table['exptime'])]
    is_in_shortest = is_good & (table['image'] == shortest_image)
    shortest_rate_by_star = dict(zip(table['star'][is_in_shortest], rate[is_in_shortest]))
    shortest_peak_by_star = dict(zip(table['star'][is_in_shortest], table['peak'][is_in_shortest]))

    ratios = {}
    for image_number in np.unique(table['image']):
        if image_number == shortest_image:
            continue
        is_in_image = is_good & (table['image'] == image_number) & np.isin(table['star'], list(shortest_rate_by_star))
        stars = table['star'][is_in_image]
        ratio = rate[is_in_image] / np.array([shortest_rate_by_star[star] for star in stars])
        peak = table['peak'][is_in_image]
        peak_in_shortest = np.array([shortest_peak_by_star[star] for star in stars])
        is_both_reference = ((peak >= reference_peak_range[0]) & (peak <= reference_peak_range[1])
                             & (peak_in_shortest >= reference_peak_range[0]) & (peak_in_shortest <= reference_peak_range[1]))
        if is_both_reference.sum() < 3:
            continue
        ratio /= np.median(ratio[is_both_reference])
        ratios[int(image_number)] = dict(peak=peak, ratio=ratio, is_saturated=table['is_saturated'][is_in_image])
    return ratios, int(shortest_image)


def explain_image_throughput(exposure_times, airmasses, image_throughput):
    """ Two alternative explanations of the per-image throughput Z_j, each fitted on its own:

    shutter     Z_j = A * (1 + delta_t / t_j)       an exposure time error delta_t changes the apparent rate by delta_t / t
    extinction  Z_j = B * 10^(-0.4 * k * X_j)        atmospheric extinction k [mag / airmass]

    If airmass and 1 / t are strongly correlated (e.g. a monotonic sequence during rising or setting), the two cannot be
    separated and the image sequence should be repeated in a shuffled exposure time order.
    """
    if len(np.unique(exposure_times)) < 3:
        return None
    shutter_design = np.c_[np.ones_like(exposure_times), 1 / exposure_times]
    (shutter_amplitude, amplitude_times_delta), *_ = np.linalg.lstsq(shutter_design, image_throughput, rcond=None)

    # ln Z = ln B - 0.4 ln(10) k X
    extinction_design = np.c_[np.ones_like(airmasses), -0.4 * np.log(10) * airmasses]
    (log_extinction_amplitude, extinction_coefficient), *_ = np.linalg.lstsq(extinction_design, np.log(image_throughput),
                                                                             rcond=None)
    correlation = float(np.corrcoef(airmasses, 1 / exposure_times)[0, 1]) if np.ptp(airmasses) > 0 else 0.
    shutter_residual = image_throughput - shutter_design @ [shutter_amplitude, amplitude_times_delta]
    extinction_residual = image_throughput - np.exp(extinction_design @ [log_extinction_amplitude, extinction_coefficient])
    return dict(delta_t=float(amplitude_times_delta / shutter_amplitude), shutter_amplitude=float(shutter_amplitude),
                shutter_rms=float(np.std(shutter_residual)),
                extinction_k=float(extinction_coefficient), extinction_amplitude=float(np.exp(log_extinction_amplitude)),
                extinction_rms=float(np.std(extinction_residual)),
                airmass_exptime_correlation=correlation, is_degenerate=bool(abs(correlation) > 0.9))


def weighted_linear_fit(design, values, value_errors):
    """ Weighted linear least squares. Parameter covariance is scaled by the reduced chi^2 if that exceeds 1.

    :return: parameters, covariance, reduced chi^2
    """
    weight = 1 / value_errors
    weighted_design = design * weight[:, None]
    parameters, *_ = np.linalg.lstsq(weighted_design, values * weight, rcond=None)
    degrees_of_freedom = len(values) - design.shape[1]
    chi2 = np.sum(((values - design @ parameters) * weight) ** 2)
    reduced_chi2 = chi2 / degrees_of_freedom if degrees_of_freedom > 0 else np.nan
    covariance = np.linalg.pinv(weighted_design.T @ weighted_design)
    if degrees_of_freedom > 0:
        covariance *= max(reduced_chi2, 1)
    return parameters, covariance, reduced_chi2


def test_residual_shutter_offset(exposure_times, airmasses, image_throughput, relative_throughput_error,
                                 extinction_coefficient, is_extinction_given, significance_threshold=3.):
    """ Is a shutter (exposure time) offset still needed after the extinction correction?

    The model is linearised for small offsets:  ln Z_j = a - 0.4 ln(10) k X_j + delta_t / t_j.

    fixed k   k as used for the cloud test; fit a + delta_t / t_j to the extinction-corrected throughput, and compare
              with a constant (no shutter offset)
    joint     k and delta_t fitted together (needs >= 4 images). If k was fitted from the same images without a shutter
              term, it may have absorbed part of the shutter effect; the joint fit and the correlation coefficient
              between k and delta_t show whether the image sequence can separate the two at all.

    :return: dict with the fixed-k and the joint results, or None if fewer than 3 distinct exposure times
    """
    if len(np.unique(np.round(exposure_times, 1))) < 3:
        return None
    log_throughput = np.log(image_throughput)
    extinction_slope = -0.4 * np.log(10) * airmasses
    inverse_exposure_time = 1 / exposure_times

    # Fixed k: extinction-corrected throughput vs 1 / t
    log_corrected_throughput = log_throughput - extinction_coefficient * extinction_slope
    (offset, delta_t), covariance, reduced_chi2 = weighted_linear_fit(
        np.c_[np.ones_like(exposure_times), inverse_exposure_time], log_corrected_throughput, relative_throughput_error)
    delta_t_error = float(np.sqrt(covariance[1, 1]))
    _, _, reduced_chi2_constant = weighted_linear_fit(np.ones((len(exposure_times), 1)), log_corrected_throughput,
                                                      relative_throughput_error)
    fixed_k = dict(extinction_coefficient=float(extinction_coefficient), delta_t=float(delta_t), delta_t_error=delta_t_error,
                   significance=float(abs(delta_t) / delta_t_error), amplitude=float(np.exp(offset)),
                   reduced_chi2_with_shutter=float(reduced_chi2), reduced_chi2_without_shutter=float(reduced_chi2_constant),
                   effect_on_shortest_exposure_percent=float(100 * delta_t / exposure_times.min()),
                   is_significant=bool(abs(delta_t) / delta_t_error > significance_threshold))

    # Joint fit of k and delta_t
    joint = None
    if len(exposure_times) >= 4:
        (joint_offset, joint_k, joint_delta_t), joint_covariance, joint_reduced_chi2 = weighted_linear_fit(
            np.c_[np.ones_like(exposure_times), extinction_slope, inverse_exposure_time], log_throughput,
            relative_throughput_error)
        joint_errors = np.sqrt(np.diag(joint_covariance))
        correlation = float(joint_covariance[1, 2] / (joint_errors[1] * joint_errors[2]))
        joint = dict(extinction_coefficient=float(joint_k), extinction_coefficient_error=float(joint_errors[1]),
                     delta_t=float(joint_delta_t), delta_t_error=float(joint_errors[2]),
                     significance=float(abs(joint_delta_t) / joint_errors[2]), k_delta_t_correlation=correlation,
                     reduced_chi2=float(joint_reduced_chi2), degrees_of_freedom=int(len(exposure_times) - 3),
                     is_significant=bool(abs(joint_delta_t) / joint_errors[2] > significance_threshold),
                     is_degenerate=bool(abs(correlation) > 0.9))
    return dict(fixed_k=fixed_k, joint=joint, is_extinction_given=bool(is_extinction_given),
                significance_threshold=significance_threshold)


def fit_extinction_coefficient(airmasses, image_throughput):
    """ Extinction coefficient k [mag / airmass] from  ln Z_j = ln B - 0.4 ln(10) k X_j. None without airmass spread."""
    if len(airmasses) < 2 or np.ptp(airmasses) < 0.01:
        return None
    design = np.c_[np.ones_like(airmasses), -0.4 * np.log(10) * airmasses]
    (_, extinction_coefficient), *_ = np.linalg.lstsq(design, np.log(image_throughput), rcond=None)
    return float(extinction_coefficient)


def group_by_exposure_time(exposure_times, relative_tolerance=0.01):
    """ Indices of images grouped by (nearly) equal exposure time, e.g. 120.0 s and 120.001 s belong together."""
    groups = []
    for image_number in np.argsort(exposure_times):
        if groups and np.isclose(exposure_times[image_number], exposure_times[groups[-1][0]], rtol=relative_tolerance):
            groups[-1].append(int(image_number))
        else:
            groups.append([int(image_number)])
    return groups


def check_against_airmass_law(airmasses, image_throughput, relative_throughput_error, extinction_coefficient, tolerance,
                              is_reference_image):
    """ Test every image, regardless of its exposure time, against the extinction law.

    Each throughput is corrected for extinction,  Z_j * 10^(0.4 k X_j);  the reference level is the weighted mean of the
    corrected throughput of the reference (unflagged) images. Clouds can only remove light, so an image is flagged when
    its deficit below the reference exceeds both 3 sigma (image and reference errors combined) and the tolerance.

    Unlike the repeated exposure time test, this test is not immune to a shutter offset: an offset that makes short
    exposures appear fainter also shows up as a deficit.

    :return: dict with reference level, per-image deficit, rejection threshold and flags
    """
    if not is_reference_image.any():
        is_reference_image = np.ones(len(image_throughput), bool)
    corrected_throughput = image_throughput * 10 ** (0.4 * extinction_coefficient * airmasses)
    reference_weight = 1 / relative_throughput_error[is_reference_image] ** 2
    reference_level = np.sum(corrected_throughput[is_reference_image] * reference_weight) / np.sum(reference_weight)
    reference_relative_error = 1 / np.sqrt(np.sum(reference_weight))

    deficit = 1 - corrected_throughput / reference_level
    threshold = np.maximum(3 * np.sqrt(relative_throughput_error ** 2 + reference_relative_error ** 2), tolerance)
    return dict(reference_level=float(reference_level), reference_relative_error=float(reference_relative_error),
                deficit=deficit, threshold=threshold, is_flagged=deficit > threshold)


def check_repeated_exposure_times(exposure_times, airmasses, image_throughput, relative_throughput_error,
                                  extinction_coefficient, tolerance):
    """ Test whether images with the same exposure time differ in throughput only by the airmass term.

    Within such a group a shutter or exposure time error is identical, so after correcting each throughput for
    extinction,  Z_j * 10^(0.4 k X_j),  all values must agree. Clouds can only remove light: the image with the highest
    corrected throughput is taken as the clearest, and every other image whose deficit exceeds both 3 sigma (combined
    throughput errors) and the tolerance is flagged as questionable.

    :return: list of per-group dicts, boolean array 'is_questionable' per image, cloud deficit per image (NaN if untested)
    """
    is_questionable = np.zeros(len(exposure_times), bool)
    cloud_deficit = np.full(len(exposure_times), np.nan)
    group_reports = []
    for group in group_by_exposure_time(exposure_times):
        if len(group) < 2:
            continue
        group = np.array(group)
        airmass_corrected_throughput = image_throughput[group] * 10 ** (0.4 * extinction_coefficient * airmasses[group])
        clearest = group[np.argmax(airmass_corrected_throughput)]
        clearest_throughput = airmass_corrected_throughput.max()

        deficit = 1 - airmass_corrected_throughput / clearest_throughput
        combined_error = np.sqrt(relative_throughput_error[group] ** 2 + relative_throughput_error[clearest] ** 2)
        threshold = np.maximum(3 * combined_error, tolerance)
        is_flagged = (deficit > threshold) & (group != clearest)

        is_questionable[group] = is_flagged
        cloud_deficit[group] = deficit
        group_reports.append(dict(exptime=float(np.median(exposure_times[group])), images=group.tolist(),
                                  clearest_image=int(clearest), deficit_percent=(deficit * 100).tolist(),
                                  threshold_percent=(threshold * 100).tolist(), questionable_images=group[is_flagged].tolist(),
                                  is_consistent_with_airmass=bool(not is_flagged.any())))
    return group_reports, is_questionable, cloud_deficit


# ----------------------------------------------------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------------------------------------------------
def plot_diagnostics(table, image_metadata_list, results, normalisations, throughput_normalisation, deviations,
                     polynomial_fit, pairwise, args,
                     output_filename):
    """ Nine-panel diagnostic figure of the linearity analysis (main aperture unless stated otherwise)."""
    main_radius = args.main_aperture
    deviation, deviation_error, expected_rate = deviations[main_radius]
    normalisation = normalisations[main_radius]
    exposure_times = np.array([metadata['exptime'] for metadata in image_metadata_list])
    saturation_level = image_metadata_list[0]['saturate']
    maxlin_level = image_metadata_list[0]['maxlin']

    is_plotted = table['is_usable'] & np.isfinite(deviation)
    is_unsaturated = is_plotted & ~table['is_saturated']
    is_saturated = is_plotted & table['is_saturated']
    image_colours = plt.cm.viridis(np.linspace(0, 0.9, len(image_metadata_list)))
    peak_bin_edges = np.logspace(np.log10(max(np.nanmin(table['peak'][is_plotted]), 10)), np.log10(1e5), 25)

    fig, axes = plt.subplots(3, 3, figsize=(21, 17))

    def mark_detector_levels(axis):
        axis.axvspan(*args.reference_peak_range, color='green', alpha=0.08, label='reference range')
        axis.axvline(maxlin_level, color='orange', ls='--', lw=1, label=f'MAXLIN {maxlin_level:.0f}')
        axis.axvline(saturation_level, color='red', ls='--', lw=1, label=f'SATURATE {saturation_level:.0f}')

    # Panel 1: deviation from linearity vs peak pixel level
    axis = axes[0, 0]
    for image_number, metadata in enumerate(image_metadata_list):
        is_this_image = is_unsaturated & (table['image'] == image_number)
        axis.errorbar(table['peak'][is_this_image], deviation[is_this_image] * 100, deviation_error[is_this_image] * 100,
                      fmt='.', ms=4, alpha=0.5, color=image_colours[image_number], elinewidth=0.5,
                      label=f"{metadata['exptime']:.0f} s")
    axis.plot(table['peak'][is_saturated], deviation[is_saturated] * 100, 'rx', ms=5, label='saturated pixels in aperture')
    bin_centres, medians, median_errors, _ = binned_median(table['peak'][is_unsaturated], deviation[is_unsaturated],
                                                           peak_bin_edges)
    axis.errorbar(bin_centres, medians * 100, median_errors * 100, fmt='ks', ms=5, label='binned median')
    level_grid = np.logspace(np.log10(peak_bin_edges[0]), np.log10(maxlin_level), 200)
    fitted_deviation, fitted_error = polynomial_fit['evaluate'](level_grid)
    axis.plot(level_grid, fitted_deviation * 100, 'k-', lw=1.5, label=f'polynomial fit (degree {args.polynomial_degree})')
    axis.fill_between(level_grid, (fitted_deviation - fitted_error) * 100, (fitted_deviation + fitted_error) * 100,
                      color='k', alpha=0.15)
    mark_detector_levels(axis)
    # Sky level of each image (in the image's colour): the floor of the peak pixel values
    sky_levels = [metadata['sky_level'] for metadata in image_metadata_list]
    for image_number, sky_level in enumerate(sky_levels):
        axis.axvline(sky_level, color=image_colours[image_number], ls='-.', lw=1.2,
                     label=f"sky level per image ({min(sky_levels):.0f}–{max(sky_levels):.0f} e-)" if image_number == 0 else None)
    axis.axhline(0, color='k', lw=0.5)
    axis.set_xscale('log')
    axis.set_xlim(left=min(axis.get_xlim()[0], 0.7 * min(sky_levels)))
    axis.set_ylim(-args.plot_range_percent, args.plot_range_percent)
    axis.set_xlabel('peak pixel value incl. sky [e-]')
    axis.set_ylabel('rate / (R_i · Z_j) − 1  [%]')
    axis.set_title(f'deviation from linearity, aperture r = {main_radius:g} px')
    axis.legend(fontsize=7, ncol=2)
    axis.grid(alpha=0.3)

    # Panel 2: deviation vs total counts in the aperture
    axis = axes[0, 1]
    total_counts = table[f'flux_r{main_radius:g}']
    for image_number, metadata in enumerate(image_metadata_list):
        is_this_image = is_unsaturated & (table['image'] == image_number)
        axis.plot(total_counts[is_this_image], deviation[is_this_image] * 100, '.', ms=4, alpha=0.5,
                  color=image_colours[image_number], label=f"{metadata['exptime']:.0f} s")
    axis.plot(total_counts[is_saturated], deviation[is_saturated] * 100, 'rx', ms=5, label='saturated')
    count_bin_edges = np.logspace(np.log10(np.nanmin(total_counts[is_unsaturated])),
                                  np.log10(np.nanmax(total_counts[is_unsaturated])), 25)
    bin_centres, medians, median_errors, _ = binned_median(total_counts[is_unsaturated], deviation[is_unsaturated],
                                                           count_bin_edges)
    axis.errorbar(bin_centres, medians * 100, median_errors * 100, fmt='ks', ms=5, label='binned median')
    axis.axhline(0, color='k', lw=0.5)
    axis.set_xscale('log')
    axis.set_ylim(-args.plot_range_percent, args.plot_range_percent)
    axis.set_xlabel(f'total counts in aperture r = {main_radius:g} px [e-]')
    axis.set_ylabel('deviation [%]')
    axis.set_title('deviation vs integrated signal (faint-end / charge-transfer effects)')
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3)

    # Panel 3: measured vs expected counts (classic linearity curve)
    axis = axes[0, 2]
    expected_counts = expected_rate * table['exptime']
    for image_number, metadata in enumerate(image_metadata_list):
        is_this_image = is_plotted & (table['image'] == image_number)
        axis.loglog(expected_counts[is_this_image], total_counts[is_this_image], '.', ms=4, alpha=0.6,
                    color=image_colours[image_number], label=f"{metadata['exptime']:.0f} s")
    axis.loglog(expected_counts[is_saturated], total_counts[is_saturated], 'rx', ms=5, label='saturated')
    count_limits = [np.nanmin(expected_counts[is_plotted]), np.nanmax(expected_counts[is_plotted])]
    axis.loglog(count_limits, count_limits, 'k--', lw=1, label='1:1')
    axis.set_xlabel('expected counts R_i · Z_j · t_j [e-]')
    axis.set_ylabel('measured counts [e-]')
    axis.set_title('measured vs expected counts')
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3, which='both')

    # Panel 4: pairwise ratios relative to the shortest exposure
    axis = axes[1, 0]
    ratios, shortest_image = pairwise
    for image_number, ratio in ratios.items():
        axis.plot(ratio['peak'][~ratio['is_saturated']], (ratio['ratio'][~ratio['is_saturated']] - 1) * 100, '.', ms=4,
                  alpha=0.6, color=image_colours[image_number],
                  label=f"{image_metadata_list[image_number]['exptime']:.0f} s / "
                        f"{image_metadata_list[shortest_image]['exptime']:.0f} s")
        axis.plot(ratio['peak'][ratio['is_saturated']], (ratio['ratio'][ratio['is_saturated']] - 1) * 100, 'x', ms=5,
                  color=image_colours[image_number])
    mark_detector_levels(axis)
    axis.axhline(0, color='k', lw=0.5)
    axis.set_xscale('log')
    axis.set_ylim(-args.plot_range_percent, args.plot_range_percent)
    axis.set_xlabel('peak pixel value in the longer exposure [e-]')
    axis.set_ylabel('normalised rate ratio − 1 [%]')
    axis.set_title('model-free check: rate ratio to the shortest exposure (x = saturated)')
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3)

    # Panel 5: aperture dependence of the binned deviation (PSF changes with level, e.g. brighter-fatter)
    axis = axes[1, 1]
    for radius in args.aperture_radii:
        radius_deviation, _, _ = deviations[radius]
        is_radius_good = table[f'is_usable_r{radius:g}'] & ~table['is_saturated'] & np.isfinite(radius_deviation)
        bin_centres, medians, median_errors, _ = binned_median(table['peak'][is_radius_good], radius_deviation[is_radius_good],
                                                               peak_bin_edges)
        axis.errorbar(bin_centres, medians * 100, median_errors * 100, fmt='o-', ms=4, capsize=2, label=f'r = {radius:g} px')
    mark_detector_levels(axis)
    axis.axhline(0, color='k', lw=0.5)
    axis.set_xscale('log')
    axis.set_ylim(-args.plot_range_percent, args.plot_range_percent)
    axis.set_xlabel('peak pixel value incl. sky [e-]')
    axis.set_ylabel('binned median deviation [%]')
    axis.set_title('aperture dependence of the deviation')
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3)

    # Panel 6: per-star tracks: normalised rate vs exposure time, coloured by peak level in the longest exposure
    axis = axes[1, 2]
    longest_image = int(np.argmax(exposure_times))
    is_in_longest = table['is_usable'] & (table['image'] == longest_image)
    peak_in_longest = dict(zip(table['star'][is_in_longest], table['peak'][is_in_longest]))
    colour_norm = matplotlib.colors.LogNorm(vmin=max(min(peak_in_longest.values(), default=100), 100), vmax=7e4)
    for star in np.unique(table['star'][is_plotted]):
        if star not in peak_in_longest:
            continue
        is_this_star = is_plotted & (table['star'] == star)
        order = np.argsort(table['exptime'][is_this_star])
        axis.plot(table['exptime'][is_this_star][order], deviation[is_this_star][order] * 100, '-', lw=0.6, alpha=0.6,
                  color=plt.cm.plasma(colour_norm(peak_in_longest[star])))
    scalar_mappable = plt.cm.ScalarMappable(norm=colour_norm, cmap='plasma')
    plt.colorbar(scalar_mappable, ax=axis, fraction=0.046, label=f'peak in the {exposure_times[longest_image]:.0f} s image [e-]')
    axis.axhline(0, color='k', lw=0.5)
    axis.set_ylim(-args.plot_range_percent, args.plot_range_percent)
    axis.set_xlabel('exposure time [s]')
    axis.set_ylabel('deviation [%]')
    axis.set_title('per-star deviation vs exposure time')
    axis.grid(alpha=0.3)

    # Panel 7: image throughput Z_j (quality subset) vs exposure time, explained either by a shutter offset or extinction
    axis = axes[2, 0]
    image_throughput = throughput_normalisation['image_throughput']
    airmasses = np.array([metadata['airmass'] for metadata in image_metadata_list])
    axis.plot(exposure_times, image_throughput, 'ko', ms=7, zorder=3,
              label=f"Z_j (quality stars: ≥ {args.throughput_detection_sigma:g} σ detection, S/N > {args.throughput_min_snr:g})")
    for image_number, metadata in enumerate(image_metadata_list):
        axis.annotate(f"X={metadata['airmass']:.3f}", (exposure_times[image_number], image_throughput[image_number]),
                      textcoords='offset points', xytext=(5, 5), fontsize=7)
    is_questionable_image = np.array([image['questionable'] for image in results['images']])
    if is_questionable_image.any():
        axis.plot(exposure_times[is_questionable_image], image_throughput[is_questionable_image], 'o', ms=14, mfc='none',
                  mec='red', mew=2, label='questionable: not explained by airmass (clouds?)', zorder=4)
    throughput_models = results['throughput_models']
    title = 'per-image throughput (transparency, airmass, shutter)'
    if throughput_models is not None:
        time_grid = np.linspace(exposure_times.min() * 0.8, exposure_times.max() * 1.1, 100)
        axis.plot(time_grid, throughput_models['shutter_amplitude'] * (1 + throughput_models['delta_t'] / time_grid), 'b-',
                  label=f"shutter: Δt = {throughput_models['delta_t'] * 1000:+.0f} ms "
                        f"(rms {throughput_models['shutter_rms'] * 100:.2f}%)")
        extinction_prediction = throughput_models['extinction_amplitude'] * 10 ** (-0.4 * throughput_models['extinction_k'] * airmasses)
        order = np.argsort(exposure_times)
        axis.plot(exposure_times[order], extinction_prediction[order], 'g^--',
                  label=f"extinction: k = {throughput_models['extinction_k']:.2f} mag/airmass "
                        f"(rms {throughput_models['extinction_rms'] * 100:.2f}%)")
        if throughput_models['is_degenerate']:
            title = (f"throughput: shutter and extinction degenerate "
                     f"(corr(X, 1/t) = {throughput_models['airmass_exptime_correlation']:+.2f})")
    axis.set_xlabel('exposure time [s]')
    axis.set_ylabel('image throughput Z_j')
    axis.set_title(title)
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3)

    # Panel 8: BANZAI catalogue flux vs our aperture photometry
    axis = axes[2, 1]
    banzai_ratio = table['banzai_flux'] / table[f'flux_r{max(args.aperture_radii):g}']
    for image_number, metadata in enumerate(image_metadata_list):
        is_this_image = is_plotted & (table['image'] == image_number) & np.isfinite(banzai_ratio)
        axis.plot(table['peak'][is_this_image], banzai_ratio[is_this_image], '.', ms=4, alpha=0.6,
                  color=image_colours[image_number], label=f"{metadata['exptime']:.0f} s")
    mark_detector_levels(axis)
    axis.axhline(1, color='k', lw=0.5)
    axis.set_xscale('log')
    axis.set_ylim(0.8, 1.2)
    axis.set_xlabel('peak pixel value incl. sky [e-]')
    axis.set_ylabel(f'BANZAI CAT flux / own flux (r = {max(args.aperture_radii):g} px)')
    axis.set_title('BANZAI catalogue vs own photometry (stars missing from CAT are not shown)')
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3)

    # Panel 9: normalised residuals of the reference measurements (validates the error model). Only stars with at least
    # three reference measurements are used; residuals are scaled by sqrt(n / (n - 1)) because R_i is fitted from the
    # same n measurements. A star with a single reference measurement would have a residual of exactly zero.
    axis = axes[2, 2]
    is_reference_used = normalisation['is_reference_used']
    reference_count_per_star = np.bincount(table['star'][is_reference_used], minlength=table['star'].max() + 1)
    reference_count = reference_count_per_star[table['star']]
    is_checked = is_reference_used & (reference_count >= 3)
    dof_correction = np.sqrt(reference_count[is_checked] / (reference_count[is_checked] - 1))
    normalised_residual = deviation[is_checked] / deviation_error[is_checked] * dof_correction
    axis.hist(normalised_residual, bins=np.linspace(-6, 6, 61), density=True, color='grey', label='reference measurements')
    gaussian_grid = np.linspace(-6, 6, 200)
    axis.plot(gaussian_grid, np.exp(-gaussian_grid ** 2 / 2) / np.sqrt(2 * np.pi), 'k-', label='unit Gaussian')
    axis.set_xlabel('deviation / error')
    axis.set_title(f"reference residuals: robust σ = {1.4826 * np.median(np.abs(normalised_residual)):.2f}")
    axis.legend(fontsize=7)

    first_image = image_metadata_list[0]
    fig.suptitle(f"{first_image['instrument']} {first_image['filter']}  {first_image['object']}  {first_image['dateobs'][:10]}: "
                 f"{results['number_of_stars']} stars, {len(image_metadata_list)} images "
                 f"({exposure_times.min():.0f}–{exposure_times.max():.0f} s)", fontsize=14)
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    plt.savefig(output_filename, dpi=80)
    plt.close()


def plot_airmass_dependence(image_metadata_list, image_throughput, relative_throughput_error, extinction_coefficient,
                            extinction_source, airmass_law_check, repeated_exposure_check, is_questionable_image,
                            residual_shutter_test, output_filename):
    """ Throughput vs airmass, and the cloud tests.

    Left: Z_j vs airmass. The grey dotted line is the extinction law of the clear images, the grey dashed line the
    rejection threshold of the airmass law test (applies to all images). For every group of equal exposure times, the
    coloured solid line is the extinction law through the group's clearest image, and the coloured dashed line the
    rejection threshold of the repeated exposure time test.
    Middle: deficit of every image below the airmass law, in observing time order, with each image's threshold.
    Right: extinction-corrected throughput vs 1 / exposure time. A shutter offset delta_t is a straight line with slope
    delta_t (in ln Z); no residual shutter offset means a horizontal line.
    """
    exposure_times = np.array([metadata['exptime'] for metadata in image_metadata_list])
    airmasses = np.array([metadata['airmass'] for metadata in image_metadata_list])
    observing_minutes = (np.array([metadata['mjd'] for metadata in image_metadata_list])
                         - np.nanmin([metadata['mjd'] for metadata in image_metadata_list])) * 24 * 60
    images_in_repeated_groups = sorted(image for group in repeated_exposure_check for image in group['images'])
    is_in_repeated_group = np.isin(np.arange(len(image_metadata_list)), images_in_repeated_groups)

    fig, (airmass_axis, deficit_axis, shutter_axis) = plt.subplots(1, 3, figsize=(24, 6.5))
    airmass_grid = np.linspace(airmasses.min() - 0.02, airmasses.max() + 0.02, 100)
    group_colours = plt.cm.tab10(np.arange(max(len(repeated_exposure_check), 1)) % 10)

    # Left: throughput vs airmass; extinction law of the clear images and the airmass law rejection threshold
    clear_law = airmass_law_check['reference_level'] * 10 ** (-0.4 * extinction_coefficient * airmass_grid)
    airmass_rejection_fraction = float(np.median(airmass_law_check['threshold']))
    airmass_axis.plot(airmass_grid, clear_law, ':', color='grey', lw=2.5,
                      label=f'extinction law of clear images, k = {extinction_coefficient:.2f} mag/airmass ({extinction_source})')
    airmass_axis.plot(airmass_grid, clear_law * (1 - airmass_rejection_fraction), '--', color='grey', lw=1.5,
                      label=f'airmass law test: rejection threshold (−{airmass_rejection_fraction * 100:.2f}%, median)')
    airmass_axis.fill_between(airmass_grid, clear_law * (1 - airmass_rejection_fraction), 0, color='grey', alpha=0.08)

    # Per group of repeated exposure times: extinction law through the clearest image, and the rejection threshold
    for group_number, group in enumerate(repeated_exposure_check):
        colour = group_colours[group_number]
        clearest = group['clearest_image']
        clearest_corrected_throughput = image_throughput[clearest] * 10 ** (0.4 * extinction_coefficient * airmasses[clearest])
        group_law = clearest_corrected_throughput * 10 ** (-0.4 * extinction_coefficient * airmass_grid)
        rejection_fraction = np.max(group['threshold_percent']) / 100
        airmass_axis.plot(airmass_grid, group_law, '-', color=colour, lw=1.5,
                          label=f"t = {group['exptime']:.0f} s: airmass law through clearest image")
        airmass_axis.plot(airmass_grid, group_law * (1 - rejection_fraction), '--', color=colour, lw=1.5,
                          label=f"t = {group['exptime']:.0f} s: repeated exposure test threshold (−{rejection_fraction * 100:.2f}%)")
        airmass_axis.fill_between(airmass_grid, group_law * (1 - rejection_fraction), 0, color=colour, alpha=0.06)
        group_images = np.array(group['images'])
        airmass_axis.errorbar(airmasses[group_images], image_throughput[group_images],
                              image_throughput[group_images] * relative_throughput_error[group_images], fmt='o', ms=8,
                              color=colour, zorder=3)

    # Images with a unique exposure time: tested against the airmass law only
    if (~is_in_repeated_group).any():
        airmass_axis.errorbar(airmasses[~is_in_repeated_group], image_throughput[~is_in_repeated_group],
                              image_throughput[~is_in_repeated_group] * relative_throughput_error[~is_in_repeated_group],
                              fmt='o', ms=7, mfc='white', color='k', zorder=3, label='unique exposure time (airmass law test only)')
    if is_questionable_image.any():
        airmass_axis.plot(airmasses[is_questionable_image], image_throughput[is_questionable_image], 'o', ms=16, mfc='none',
                          mec='red', mew=2, zorder=4, label='questionable (clouds?)')
    for image_number in range(len(image_metadata_list)):
        airmass_axis.annotate(f"{exposure_times[image_number]:.0f} s", (airmasses[image_number], image_throughput[image_number]),
                              textcoords='offset points', xytext=(8, -3), fontsize=8)
    throughput_range = np.ptp(image_throughput)
    airmass_axis.set_ylim(image_throughput.min() - max(0.5 * throughput_range, 0.01),
                          image_throughput.max() + max(0.3 * throughput_range, 0.005))
    airmass_axis.set_xlabel('airmass')
    airmass_axis.set_ylabel('image throughput Z_j')
    airmass_axis.set_title('throughput vs airmass')
    airmass_axis.legend(fontsize=7, loc='lower left')
    airmass_axis.grid(alpha=0.3)

    # Middle: deficit of every image below the airmass law of the clear images, in time order
    image_colours = plt.cm.viridis(np.linspace(0, 0.9, len(image_metadata_list)))
    deficit_percent = airmass_law_check['deficit'] * 100
    threshold_percent = airmass_law_check['threshold'] * 100
    combined_error_percent = 100 * np.sqrt(relative_throughput_error ** 2 + airmass_law_check['reference_relative_error'] ** 2)
    for image_number in range(len(image_metadata_list)):
        deficit_axis.errorbar(observing_minutes[image_number], deficit_percent[image_number], combined_error_percent[image_number],
                              fmt='o', ms=8, color=image_colours[image_number])
        deficit_axis.plot(observing_minutes[image_number], threshold_percent[image_number], '_', ms=25, mew=2.5,
                          color=image_colours[image_number])
        deficit_axis.annotate(f"{exposure_times[image_number]:.0f} s", (observing_minutes[image_number], deficit_percent[image_number]),
                              textcoords='offset points', xytext=(8, -3), fontsize=8)
    if airmass_law_check['is_flagged'].any():
        deficit_axis.plot(observing_minutes[airmass_law_check['is_flagged']], deficit_percent[airmass_law_check['is_flagged']],
                          'o', ms=16, mfc='none', mec='red', mew=2, label='questionable (clouds?)')
    deficit_axis.plot([], [], '_', ms=15, mew=2.5, color='k', label='rejection threshold per image: max(3σ, tolerance)')
    deficit_axis.axhline(0, color='k', lw=0.5)
    deficit_axis.invert_yaxis()  # light lost to clouds plotted downwards, as in the left panel
    deficit_axis.set_xlabel('time since first image [min]')
    deficit_axis.set_ylabel('airmass-corrected deficit vs clear images [%]')
    deficit_axis.set_title('cloud test: all images vs airmass law')
    deficit_axis.legend(fontsize=7)
    deficit_axis.grid(alpha=0.3)

    # Right: residual shutter offset, extinction-corrected throughput vs 1 / t (cloudy images marked, not fitted)
    inverse_exposure_time = 1 / exposure_times
    corrected_throughput = image_throughput * 10 ** (0.4 * extinction_coefficient * airmasses)
    corrected_throughput /= np.median(corrected_throughput[~is_questionable_image])
    shutter_axis.errorbar(inverse_exposure_time[~is_questionable_image], corrected_throughput[~is_questionable_image],
                          corrected_throughput[~is_questionable_image] * relative_throughput_error[~is_questionable_image],
                          fmt='ko', ms=7, label='extinction-corrected throughput')
    if is_questionable_image.any():
        shutter_axis.plot(inverse_exposure_time[is_questionable_image], corrected_throughput[is_questionable_image], 'o',
                          ms=8, color='grey')
        shutter_axis.plot(inverse_exposure_time[is_questionable_image], corrected_throughput[is_questionable_image], 'o',
                          ms=16, mfc='none', mec='red', mew=2, label='questionable (clouds?), not fitted')
    for image_number in range(len(image_metadata_list)):
        shutter_axis.annotate(f"{exposure_times[image_number]:.0f} s",
                              (inverse_exposure_time[image_number], corrected_throughput[image_number]),
                              textcoords='offset points', xytext=(8, -3), fontsize=8)
    if residual_shutter_test is not None:
        fixed_k = residual_shutter_test['fixed_k']
        inverse_time_grid = np.linspace(0, inverse_exposure_time.max() * 1.1, 100)
        clear_median = np.median((image_throughput * 10 ** (0.4 * extinction_coefficient * airmasses))[~is_questionable_image])
        shutter_axis.plot(inverse_time_grid, fixed_k['amplitude'] * np.exp(fixed_k['delta_t'] * inverse_time_grid) / clear_median,
                          'b-', label=f"Δt = {fixed_k['delta_t'] * 1000:+.0f} ± {fixed_k['delta_t_error'] * 1000:.0f} ms "
                                      f"({fixed_k['significance']:.1f} σ), k = {fixed_k['extinction_coefficient']:.2f} fixed")
        joint = residual_shutter_test['joint']
        if joint is not None:
            shutter_axis.plot([], [], ' ', label=f"joint fit: Δt = {joint['delta_t'] * 1000:+.0f} ± {joint['delta_t_error'] * 1000:.0f} ms, "
                                                 f"k = {joint['extinction_coefficient']:.2f} ± {joint['extinction_coefficient_error']:.2f}, "
                                                 f"corr = {joint['k_delta_t_correlation']:+.2f}")
        shutter_axis.set_xlim(0, inverse_exposure_time.max() * 1.1)
    shutter_axis.axhline(1, color='k', lw=0.5)
    shutter_axis.set_xlabel('1 / exposure time [1/s]')
    shutter_axis.set_ylabel('extinction-corrected throughput (normalised)')
    shutter_axis.set_title('residual shutter offset after extinction correction')
    shutter_axis.legend(fontsize=7)
    shutter_axis.grid(alpha=0.3)

    first_image = image_metadata_list[0]
    fig.suptitle(f"{first_image['instrument']} {first_image['filter']}  {first_image['object']}  "
                 f"{first_image['dateobs'][:10]}: airmass dependence of the throughput")
    plt.tight_layout()
    plt.savefig(output_filename, dpi=90)
    plt.close()


def plot_aperture_analysis(aperture_analysis, args, used_aperture, used_annulus, title, output_filename):
    """ Four panels: curves of growth, error terms vs aperture radius, sky annulus wing bias, isolation vs radius."""
    radii = aperture_analysis['radii']
    optimum_radius = aperture_analysis['optimum_radius']
    annulus = (aperture_analysis['annulus_inner'], aperture_analysis['annulus_outer'])
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    def mark_choices(axis):
        axis.axvline(optimum_radius, color='green', lw=2, label=f'optimum aperture r = {optimum_radius:g} px')
        axis.axvspan(*annulus, color='orange', alpha=0.15, label=f'recommended sky annulus {annulus[0]:g}–{annulus[1]:g} px')
        if used_aperture != optimum_radius:
            axis.axvline(used_aperture, color='grey', ls='--', lw=1.5, label=f'aperture used r = {used_aperture:g} px')

    # Panel 1: curve of growth per image, with the star-to-star spread of the deepest image
    axis = axes[0, 0]
    image_colours = plt.cm.viridis(np.linspace(0, 0.9, len(aperture_analysis['exposure_times'])))
    deepest_curves = aperture_analysis['growth_curves'][-1]
    axis.fill_between(radii, np.percentile(deepest_curves, 16, axis=0) * 100, np.percentile(deepest_curves, 84, axis=0) * 100,
                      color=image_colours[-1], alpha=0.2, label='16–84% of stars, deepest image')
    for image_number, exposure_time in enumerate(aperture_analysis['exposure_times']):
        axis.plot(radii, aperture_analysis['median_curve_per_image'][image_number] * 100, '-', color=image_colours[image_number],
                  label=f"{exposure_time:.0f} s ({aperture_analysis['number_of_growth_curve_stars'][image_number]} stars)")
    axis.axhline(99, color='k', ls=':', lw=1, label='99%')
    mark_choices(axis)
    axis.set_ylim(50, 101)
    axis.set_xlabel('aperture radius [px]')
    axis.set_ylabel(f'collected light [% of light within {radii[-1]:g} px]')
    axis.set_title('curve of growth')
    axis.legend(fontsize=7, loc='lower right')
    axis.grid(alpha=0.3)

    # Panel 2: error terms vs aperture radius
    axis = axes[0, 1]
    floor = 1e-4  # keeps zeros visible on the log axis
    axis.semilogy(radii, (1 - aperture_analysis['median_encircled_energy']) * 100 + floor, 'k:', label='light outside aperture (1 − EE)')
    axis.semilogy(radii, aperture_analysis['missed_light_uncertainty'] * 100 + floor, 'b-', lw=2,
                  label=f"missed light uncertainty = {aperture_analysis['missed_light_variability']:.2f} × (1 − EE)")
    axis.semilogy(radii, aperture_analysis['frame_to_frame_spread'] * 100 + floor, 'b:', label='missed light: frame-to-frame spread')
    axis.semilogy(radii, aperture_analysis['mean_contamination'] * 100 + floor, 'r-', lw=2, label='neighbour light: mean over stars')
    axis.semilogy(radii, aperture_analysis['median_contamination'] * 100 + floor, 'r--', lw=1, label='neighbour light: median star')
    axis.semilogy(radii, aperture_analysis['contamination_90th_percentile'] * 100 + floor, 'r:', lw=1, label='neighbour light: 90th percentile')
    axis.semilogy(radii, aperture_analysis['total_error'] * 100, 'k-', lw=2.5, label='total = quadrature sum of the thick lines')
    axis.semilogy(radii, aperture_analysis['noise'] * 100, 'm-', lw=1,
                  label=f"for information: noise, star with peak {aperture_analysis['typical_peak']:.0f}")
    axis.semilogy(radii, aperture_analysis['light_loss_scatter'] * 100 + floor, '-', color='lightblue', lw=1,
                  label='for information: star-to-star scatter of EE')
    axis.axvspan(0.75 * radii[-1], radii[-1], color='grey', alpha=0.1, label='not considered (near normalisation radius)')
    mark_choices(axis)
    axis.set_ylim(1e-3, 100)
    axis.set_xlabel('aperture radius [px]')
    axis.set_ylabel('relative error [%]')
    axis.set_title('aperture trade-off: missed light vs neighbour light')
    axis.legend(fontsize=6.5)
    axis.grid(alpha=0.3, which='both')

    # Panel 3: wing bias of the sky annulus vs its inner radius, for the optimum aperture
    axis = axes[1, 0]
    axis.semilogy(aperture_analysis['annulus_inner_grid'], aperture_analysis['wing_bias'] * 100 + 1e-5, 'k-',
                  label=f'star light in the annulus, as flux bias (aperture r = {optimum_radius:g} px)')
    axis.axhline(args.annulus_wing_tolerance * 100, color='red', ls='--', label=f'tolerance {args.annulus_wing_tolerance * 100:g}%')
    axis.axvline(annulus[0], color='orange', lw=2, label=f'recommended inner radius {annulus[0]:g} px')
    axis.axvline(used_annulus[0], color='grey', ls='--', label=f'inner radius used {used_annulus[0]:g} px')
    axis.set_xlabel(f'annulus inner radius [px] (area = {args.annulus_area_factor:g} × aperture area)')
    axis.set_ylabel('flux bias from star wings in the annulus [%]')
    axis.set_title('sky annulus: contamination by the star itself')
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3, which='both')

    # Panel 4: stars kept by the isolation cut vs annulus outer radius, and stars with significant neighbour light
    axis = axes[1, 1]
    axis.plot(aperture_analysis['isolation_radius_grid'], aperture_analysis['isolated_fraction_by_radius'] * 100, 'k-',
              label='stars with no neighbour within radius (isolation cut at annulus outer radius)')
    axis.plot(radii, aperture_analysis['fraction_contaminated'] * 100, 'r-',
              label=f'stars with neighbour light > {args.contamination_tolerance * 100:g}% inside aperture radius')
    axis.axvline(annulus[1], color='orange', lw=2, label=f'recommended annulus outer radius {annulus[1]:g} px')
    axis.axvline(used_annulus[1], color='grey', ls='--', label=f'annulus outer radius used {used_annulus[1]:g} px')
    axis.axvline(optimum_radius, color='green', lw=2)
    axis.set_ylim(0, 101)
    axis.set_xlabel('radius [px]')
    axis.set_ylabel('stars [%]')
    axis.set_title('crowding: isolation and contaminated apertures')
    axis.legend(fontsize=7)
    axis.grid(alpha=0.3)

    fig.suptitle(f'{title}: aperture optimisation')
    plt.tight_layout()
    plt.savefig(output_filename, dpi=80)
    plt.close()


def least_squares_with_scatter(design, values):
    """ Unweighted linear least squares; parameter errors from the scatter of the residuals.

    :return: parameters, parameter errors (NaN without degrees of freedom), residuals, degrees of freedom
    """
    parameters, *_ = np.linalg.lstsq(design, values, rcond=None)
    residuals = values - design @ parameters
    degrees_of_freedom = len(values) - design.shape[1]
    if degrees_of_freedom > 0:
        residual_variance = np.sum(residuals ** 2) / degrees_of_freedom
        parameter_errors = np.sqrt(np.diag(np.linalg.pinv(design.T @ design)) * residual_variance)
    else:
        parameter_errors = np.full(design.shape[1], np.nan)
    return parameters, parameter_errors, residuals, degrees_of_freedom


def test_sky_scaling(image_metadata_list, significance_threshold=3.):
    """ Does the sky level scale with exposure time?

    For a linear detector and a constant sky, the sky level is proportional to the exposure time, S_j = r * t_j. Fitted:
      proportional   S = r t
      offset         S = r t + b             b: additive level present in every image, independent of the exposure time
                                             (e.g. an error of the subtracted bias level, or a dark residual)
      offset + drift S = (r0 + r1 tau) t + b sky brightness changing with time tau [h] since the first image; needs
                                             repeated exposure times or enough images to separate drift and offset
      curvature      S = r t + b + q t^2     a significant q means the sky signal is not linear in t (non-linearity at sky
                                             levels, or a sky change correlated with the exposure time)
    Parameter errors come from the scatter of the residuals; |value| / error > significance_threshold is significant.

    :return: dict of fit results, or None if fewer than 3 images or fewer than 2 distinct exposure times
    """
    exposure_times = np.array([metadata['exptime'] for metadata in image_metadata_list])
    sky_levels = np.array([metadata['sky_level'] for metadata in image_metadata_list])
    mjd = np.array([metadata['mjd'] for metadata in image_metadata_list], float)
    hours_since_first = (mjd - np.nanmin(mjd)) * 24 if np.isfinite(mjd).all() else np.zeros(len(mjd))
    if len(exposure_times) < 3 or len(np.unique(np.round(exposure_times, 1))) < 2:
        return None

    def describe(name, design, parameter_names):
        parameters, errors, residuals, degrees_of_freedom = least_squares_with_scatter(design, sky_levels)
        return dict(model=name, degrees_of_freedom=int(degrees_of_freedom), rms_residual=float(np.std(residuals)),
                    parameters={parameter: dict(value=float(value), error=float(error),
                                                significance=float(abs(value) / error) if error > 0 else np.nan)
                                for parameter, value, error in zip(parameter_names, parameters, errors)},
                    residuals=residuals.tolist())

    models = {
        'proportional': describe('S = r t', exposure_times[:, None], ['r']),
        'offset': describe('S = r t + b', np.c_[exposure_times, np.ones_like(exposure_times)], ['r', 'b']),
        'curvature': describe('S = r t + b + q t^2', np.c_[exposure_times, np.ones_like(exposure_times), exposure_times ** 2],
                              ['r', 'b', 'q']) if len(np.unique(np.round(exposure_times, 1))) >= 3 and len(exposure_times) >= 4 else None,
        'offset_drift': describe('S = (r0 + r1 tau) t + b', np.c_[exposure_times, exposure_times * hours_since_first,
                                                                   np.ones_like(exposure_times)], ['r0', 'r1', 'b'])
        if len(exposure_times) >= 4 and np.ptp(hours_since_first) > 0 else None,
    }
    offset = models['offset']['parameters']['b']
    curvature = models['curvature']['parameters']['q'] if models['curvature'] else None
    return dict(exposure_times=exposure_times.tolist(), sky_levels=sky_levels.tolist(),
                hours_since_first=hours_since_first.tolist(), models=models,
                is_offset_significant=bool(offset['significance'] > significance_threshold),
                is_curvature_significant=bool(curvature is not None and curvature['significance'] > significance_threshold),
                significance_threshold=significance_threshold)


def plot_sky_scaling(sky_test, title, output_filename):
    """ Left: sky level vs exposure time with the proportional and the offset model.
    Right: sky rate (sky level / t) vs time, raw and with the fitted offset removed."""
    exposure_times = np.array(sky_test['exposure_times'])
    sky_levels = np.array(sky_test['sky_levels'])
    hours = np.array(sky_test['hours_since_first'])
    offset_model = sky_test['models']['offset']['parameters']
    proportional_rate = sky_test['models']['proportional']['parameters']['r']['value']
    rate, offset = offset_model['r']['value'], offset_model['b']['value']

    fig, (level_axis, rate_axis) = plt.subplots(1, 2, figsize=(14, 5))
    time_grid = np.linspace(0, exposure_times.max() * 1.05, 100)
    level_axis.plot(exposure_times, sky_levels, 'ko', ms=6, label='sky level (median of the image)')
    level_axis.plot(time_grid, proportional_rate * time_grid, 'b--', lw=1, label=f'S = r t, r = {proportional_rate:.3f} e-/px/s')
    level_axis.plot(time_grid, rate * time_grid + offset, 'r-', lw=1.5,
                    label=f"S = r t + b, r = {rate:.3f} e-/px/s, b = {offset:+.1f} ± {offset_model['b']['error']:.1f} e-")
    level_axis.axhline(0, color='k', lw=0.5)
    level_axis.set_xlabel('exposure time [s]')
    level_axis.set_ylabel('sky level [e-/px]')
    level_axis.set_title('sky level vs exposure time')
    level_axis.legend(fontsize=8)
    level_axis.grid(alpha=0.3)

    # Group by exposure time rounded to 0.1 s (nominally equal exposure times differ in the third decimal)
    rounded_exposure_times = np.round(exposure_times, 1)
    image_colours = plt.cm.viridis(np.linspace(0, 0.9, len(np.unique(rounded_exposure_times))))
    colour_by_time = dict(zip(np.unique(rounded_exposure_times), image_colours))
    for exposure_time in np.unique(rounded_exposure_times):
        is_this_time = rounded_exposure_times == exposure_time
        rate_axis.plot(hours[is_this_time] * 60, sky_levels[is_this_time] / exposure_times[is_this_time], 'o', mfc='none',
                       color=colour_by_time[exposure_time])
        rate_axis.plot(hours[is_this_time] * 60, (sky_levels[is_this_time] - offset) / exposure_times[is_this_time], 'o',
                       color=colour_by_time[exposure_time], label=f'{exposure_time:.0f} s')
    rate_axis.plot([], [], 'ko', mfc='none', label='S / t (raw)')
    rate_axis.plot([], [], 'ko', label='(S - b) / t (offset removed)')
    rate_axis.set_xlabel('time since first image [min]')
    rate_axis.set_ylabel('sky rate [e-/px/s]')
    rate_axis.set_title('sky rate: should depend on time, not on exposure time')
    rate_axis.legend(fontsize=8, ncol=2)
    rate_axis.grid(alpha=0.3)

    fig.suptitle(f'{title}: does the sky brightness scale with exposure time?')
    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()


def register_images(files, reference_sources, args):
    """ WCS of every image, registered against the stars of the reference (last, deepest) image where needed.

    :return: dict filename -> WCS, without the images that could not be registered
    """
    reference_file = files[-1]
    is_reference_solved = frame_has_valid_astrometry(reference_file, fits.getheader(reference_file, 'SCI'))
    image_wcs = {}
    for filename in files:
        header = fits.getheader(filename, 'SCI')
        header_wcs = WCS(header)
        if filename == reference_file or (is_reference_solved and frame_has_valid_astrometry(filename, header)):
            image_wcs[filename] = header_wcs
            continue
        image_sources = detect_sources(filename, args)
        registration = register_frame(image_sources['x'], image_sources['y'], header_wcs, reference_sources['sky_positions'])
        if registration is None:
            _log.warning(f"{os.path.basename(filename)}: registration against the reference image failed; image dropped")
            continue
        image_wcs[filename] = registration['wcs']
        _log.info(f"{os.path.basename(filename)}: registered with {registration['number_of_matches']} stars, "
                  f"rms {registration['rms_arcsec']:.2f}\", offset E {registration['pointing_offset_arcsec'][0]:+.1f}\" "
                  f"N {registration['pointing_offset_arcsec'][1]:+.1f}\"")
    if not is_reference_solved:
        _log.info(f"reference image {os.path.basename(reference_file)} has no astrometric solution: images registered "
                  f"relative to it")
    return image_wcs


# ----------------------------------------------------------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------------------------------------------------------
def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='On-sky linearity test from a series of BANZAI images of a star field with increasing exposure times',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('files', nargs='+', help='BANZAI e91 images of one field, same filter, varying EXPTIME')
    parser.add_argument('--outdir', default='plots/linearity', help='output directory for plots and tables')
    parser.add_argument('--aperture-radii', type=float, nargs='+', default=[8, 15, 25],
                        help='fixed aperture radii [pixels]')
    parser.add_argument('--main-aperture', type=float, default=15,
                        help='aperture radius [pixels] used for the headline result; must be one of --aperture-radii')
    parser.add_argument('--annulus', type=float, nargs=2, default=[20, 25],
                        help='inner and outer radius of the local sky annulus [pixels]')
    parser.add_argument('--reference-peak-range', type=float, nargs=2, default=[100, 20000],
                        help='peak pixel range [e-] that defines linear behaviour (deviation = 0)')
    parser.add_argument('--polynomial-degree', type=int, default=2, help='degree of the non-linearity polynomial')
    parser.add_argument('--detection-threshold', type=float, default=5, help='star detection threshold [sigma]')
    parser.add_argument('--isolation-flux-ratio', type=float, default=0.01,
                        help='neighbours fainter than this fraction of a star\'s flux do not count in the isolation cut; '
                             'their light is the same fraction of the star in every image and cancels in the rate comparison')
    parser.add_argument('--min-snr', type=float, default=20, help='minimum signal-to-noise of a measurement')
    parser.add_argument('--max-centroid-shift', type=float, default=2,
                        help='maximum re-centring shift [pixels] before falling back to the WCS position')
    parser.add_argument('--plot-range-percent', type=float, default=5, help='y-axis range of deviation plots [%%]')
    parser.add_argument('--extinction-coefficient', type=float, default=None,
                        help='extinction k [mag/airmass] for the repeated exposure time check; fitted from all images if not given')
    parser.add_argument('--cloud-tolerance', type=float, default=0.005,
                        help='throughput deficit (fraction) beyond the airmass term tolerated before an image is questionable')
    parser.add_argument('--throughput-detection-sigma', type=float, default=50,
                        help='minimum detection significance [sigma] of stars used for the airmass, cloud and shutter tests')
    parser.add_argument('--throughput-min-snr', type=float, default=50,
                        help='minimum signal-to-noise of measurements used for the airmass, cloud and shutter tests')
    parser.add_argument('--exclude-questionable', action='store_true',
                        help='exclude questionable images (possible clouds) from the normalisation and the fit')
    parser.add_argument('--auto-aperture', action='store_true',
                        help='use the optimum aperture and the recommended sky annulus from the aperture analysis')
    parser.add_argument('--contamination-tolerance', type=float, default=0.001,
                        help='neighbour light inside the aperture (fraction of the star flux) counted as contaminated')
    parser.add_argument('--annulus-wing-tolerance', type=float, default=0.001,
                        help='maximum flux bias (fraction) from the star\'s own light in the sky annulus')
    parser.add_argument('--annulus-area-factor', type=float, default=5,
                        help='sky annulus area in units of the aperture area')
    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    if args.main_aperture not in args.aperture_radii:
        parser.error('--main-aperture must be one of --aperture-radii')
    return args


def main():
    args = parseCommandLine()
    os.makedirs(args.outdir, exist_ok=True)

    # Sort images by exposure time; detect stars on the deepest image, the astrometric reference
    exposure_time_by_file = {filename: float(fits.getval(filename, 'EXPTIME', extname='SCI')) for filename in args.files}
    files = sorted(args.files, key=lambda filename: exposure_time_by_file[filename])
    filters = {fits.getval(filename, 'FILTER', extname='SCI') for filename in files}
    if len(filters) > 1:
        _log.warning(f"images have different filters {filters}; colour-dependent throughput will add scatter")
    sources = detect_sources(files[-1], args)

    # WCS of every image: images without a BANZAI astrometric solution are registered against the stars of the
    # reference image. If the reference itself is unsolved, all images are registered relative to it, which keeps the
    # star positions consistent between images (all that the forced photometry needs). Unregistrable images are dropped.
    image_wcs = register_images(files, sources, args)
    files = [filename for filename in files if filename in image_wcs]
    if len(files) < 2:
        _log.error(f"only {len(files)} image(s) with a usable WCS; at least 2 are needed")
        sys.exit(1)

    # Aperture optimisation; with --auto-aperture its choice is used for the photometry
    aperture_analysis = analyse_apertures(files, sources, args, image_wcs)
    if args.auto_aperture:
        args.main_aperture = aperture_analysis['optimum_radius']
        if args.main_aperture not in args.aperture_radii:
            args.aperture_radii = sorted(args.aperture_radii + [args.main_aperture])
        args.annulus = [aperture_analysis['annulus_inner'], aperture_analysis['annulus_outer']]
        _log.info(f"using optimum aperture r = {args.main_aperture:g} px, sky annulus {args.annulus[0]:g}-{args.annulus[1]:g} px")
    star_positions, star_detection_significance = select_isolated_stars(sources, args.annulus[1], args.isolation_flux_ratio)

    # Forced photometry on every image
    per_image_measurements, image_metadata_list = [], []
    for filename in files:
        measurements, metadata = measure_image(filename, star_positions, args, image_wcs[filename])
        per_image_measurements.append(measurements)
        image_metadata_list.append(metadata)
    table = build_measurement_table(per_image_measurements, image_metadata_list, args)
    table['detection_significance'] = star_detection_significance[table['star']]

    # Normalisation (R_i, Z_j) and deviations for every aperture
    normalisations, deviations = {}, {}
    for radius in args.aperture_radii:
        normalisations[radius] = fit_star_rates_and_image_throughput(table, radius, args.reference_peak_range)
        deviations[radius] = deviation_from_linearity(table, radius, normalisations[radius])

    # Image throughput for the airmass, cloud and shutter tests: from a separate high-quality subset of strongly detected
    # stars with precise measurements, so that faint-star noise does not enter these tests
    with np.errstate(divide='ignore', invalid='ignore'):
        main_signal_to_noise = table[f'flux_r{args.main_aperture:g}'] / table[f'flux_error_r{args.main_aperture:g}']
    table['is_throughput_quality'] = ((table['detection_significance'] >= args.throughput_detection_sigma)
                                      & (main_signal_to_noise > args.throughput_min_snr))
    throughput_normalisation = fit_star_rates_and_image_throughput(table, args.main_aperture, args.reference_peak_range,
                                                                   extra_selection=table['is_throughput_quality'])
    quality_measurements_per_image = np.bincount(table['image'][throughput_normalisation['is_reference_used']],
                                                 minlength=len(image_metadata_list))
    _log.info(f"throughput quality stars (detection >= {args.throughput_detection_sigma:g} sigma, "
              f"S/N > {args.throughput_min_snr:g}) per image: {quality_measurements_per_image.tolist()}")
    if (quality_measurements_per_image < 5).any():
        _log.warning("fewer than 5 throughput quality stars in some images; airmass and shutter tests are unreliable")

    # Cloud tests: is the throughput variation explained by airmass alone? Two tests, an image failing either one is
    # questionable:
    #   airmass law          every image against the extinction law of the clear images
    #   repeated exposures   images with the same exposure time against each other (immune to a shutter offset)
    # A fitted k and the reference level are refitted without the questionable images until the flags no longer
    # change, so that a cloudy image does not bias the law it is tested against.
    exposure_times = np.array([metadata['exptime'] for metadata in image_metadata_list])
    airmasses = np.array([metadata['airmass'] for metadata in image_metadata_list])
    is_questionable_image = np.zeros(len(image_metadata_list), bool)
    for check_iteration in range(10):
        if args.extinction_coefficient is not None:
            extinction_coefficient, extinction_source = args.extinction_coefficient, 'given'
        else:
            is_clear_image = ~is_questionable_image
            extinction_coefficient = fit_extinction_coefficient(airmasses[is_clear_image],
                                                                throughput_normalisation['image_throughput'][is_clear_image])
            extinction_source = f'fitted to {is_clear_image.sum()} unflagged images'
            if extinction_coefficient is None:
                extinction_coefficient, extinction_source = 0., 'no airmass spread, k = 0 assumed'
        airmass_law_check = check_against_airmass_law(
            airmasses, throughput_normalisation['image_throughput'], throughput_normalisation['relative_throughput_error'],
            extinction_coefficient, args.cloud_tolerance, is_reference_image=~is_questionable_image)
        repeated_exposure_check, is_flagged_by_repeated_exposures, repeated_exposure_deficit = check_repeated_exposure_times(
            exposure_times, airmasses, throughput_normalisation['image_throughput'],
            throughput_normalisation['relative_throughput_error'], extinction_coefficient, args.cloud_tolerance)
        is_questionable_next = airmass_law_check['is_flagged'] | is_flagged_by_repeated_exposures
        if np.array_equal(is_questionable_next, is_questionable_image):
            break
        is_questionable_image = is_questionable_next
    questionable_reasons = [[reason for reason, is_flagged in (('airmass law', airmass_law_check['is_flagged'][image]),
                                                               ('repeated exposure time', is_flagged_by_repeated_exposures[image]))
                             if is_flagged]
                            for image in range(len(image_metadata_list))]
    table['is_questionable_image'] = is_questionable_image[table['image']]

    if args.exclude_questionable and is_questionable_image.any():
        _log.info(f"excluding questionable images {np.nonzero(is_questionable_image)[0].tolist()} and refitting")
        for radius in args.aperture_radii:
            table[f'is_usable_r{radius:g}'] &= ~table['is_questionable_image']
            normalisations[radius] = fit_star_rates_and_image_throughput(table, radius, args.reference_peak_range)
            deviations[radius] = deviation_from_linearity(table, radius, normalisations[radius])
        table['is_usable'] = table[f'is_usable_r{args.main_aperture:g}']

    # Non-linearity polynomial (main aperture), fitted to clean measurements up to MAXLIN
    deviation, deviation_error, _ = deviations[args.main_aperture]
    maxlin_level = image_metadata_list[0]['maxlin']
    is_fitted = (table['is_usable'] & ~table['is_saturated'] & np.isfinite(deviation)
                 & (table['peak'] <= maxlin_level))
    polynomial_fit = fit_nonlinearity_polynomial(table['peak'][is_fitted], deviation[is_fitted],
                                                 deviation_error[is_fitted], args.polynomial_degree)
    pairwise = pairwise_rate_ratios(table, args.main_aperture, args.reference_peak_range)

    # Residual shutter offset after extinction correction, on the unflagged images (throughputs as in the cloud test)
    is_clear_image = ~is_questionable_image
    residual_shutter_test = test_residual_shutter_offset(
        exposure_times[is_clear_image], airmasses[is_clear_image], throughput_normalisation['image_throughput'][is_clear_image],
        throughput_normalisation['relative_throughput_error'][is_clear_image], extinction_coefficient,
        is_extinction_given=args.extinction_coefficient is not None)

    # Does the sky level scale with exposure time?
    sky_scaling_test = test_sky_scaling(image_metadata_list)

    # Summary: image throughput from the linearity normalisation, and from the quality subset used for the tests
    image_throughput = normalisations[args.main_aperture]['image_throughput']
    quality_throughput = throughput_normalisation['image_throughput']
    has_star_rate = np.isfinite(normalisations[args.main_aperture]['star_rate'])
    max_fitted_peak = float(np.max(table['peak'][is_fitted]))
    report_levels = [level for level in REPORT_LEVELS if level <= max_fitted_peak]
    deviation_at_levels, deviation_error_at_levels = polynomial_fit['evaluate'](report_levels)
    results = dict(
        images=[dict(metadata, throughput=float(throughput), throughput_quality=float(quality),
                     throughput_quality_stars=int(quality_count), questionable=bool(questionable),
                     questionable_reasons=reasons,
                     airmass_law_deficit_percent=float(airmass_deficit * 100),
                     airmass_law_threshold_percent=float(airmass_threshold * 100),
                     repeated_exposure_deficit_percent=float(repeated_deficit * 100))
                for metadata, throughput, quality, quality_count, questionable, reasons, airmass_deficit, airmass_threshold,
                    repeated_deficit
                in zip(image_metadata_list, image_throughput, quality_throughput, quality_measurements_per_image,
                       is_questionable_image, questionable_reasons, airmass_law_check['deficit'],
                       airmass_law_check['threshold'], repeated_exposure_deficit)],
        airmass_law_reference_level=airmass_law_check['reference_level'],
        throughput_quality_selection=dict(detection_sigma=args.throughput_detection_sigma, min_snr=args.throughput_min_snr),
        extinction_coefficient=dict(value=extinction_coefficient, source=extinction_source),
        repeated_exposure_check=repeated_exposure_check, questionable_images_excluded=bool(args.exclude_questionable),
        number_of_stars=int(has_star_rate.sum()), number_of_fitted_measurements=int(is_fitted.sum()),
        main_aperture_px=args.main_aperture, reference_peak_range=args.reference_peak_range,
        max_fitted_peak=max_fitted_peak,
        polynomial_coefficients=polynomial_fit['coefficients'].tolist(),
        polynomial_level_scale=polynomial_fit['level_scale'], polynomial_reduced_chi2=float(polynomial_fit['reduced_chi2']),
        deviation_percent_at_level={f'{level:d}': dict(value=float(value * 100), error=float(error * 100))
                                    for level, value, error in zip(report_levels, deviation_at_levels, deviation_error_at_levels)},
        throughput_models=explain_image_throughput(exposure_times[~is_questionable_image], airmasses[~is_questionable_image],
                                                   quality_throughput[~is_questionable_image]),
        residual_shutter_test=residual_shutter_test,
        sky_scaling=sky_scaling_test,
        aperture_analysis=dict(aperture_summary(aperture_analysis), aperture_used_px=args.main_aperture,
                               annulus_used_px=list(args.annulus)))

    first_image = image_metadata_list[0]
    output_basename = os.path.join(args.outdir, f"linearity_{first_image['instrument']}_{first_image['filter']}_"
                                                f"{first_image['dateobs'][:10].replace('-', '')}")
    plot_diagnostics(table, image_metadata_list, results, normalisations, throughput_normalisation, deviations,
                     polynomial_fit, pairwise, args,
                     f'{output_basename}.pdf')
    # The cloud test was done on the throughputs before any exclusion; plot exactly those
    plot_airmass_dependence(image_metadata_list, throughput_normalisation['image_throughput'],
                            throughput_normalisation['relative_throughput_error'], extinction_coefficient, extinction_source,
                            airmass_law_check, repeated_exposure_check, is_questionable_image, residual_shutter_test,
                            f'{output_basename}_airmass.pdf')
    plot_aperture_analysis(aperture_analysis, args, args.main_aperture, args.annulus,
                           f"{first_image['instrument']} {first_image['filter']}  {first_image['object']}  {first_image['dateobs'][:10]}",
                           f'{output_basename}_aperture.pdf')
    if sky_scaling_test is not None:
        plot_sky_scaling(sky_scaling_test, f"{first_image['instrument']} {first_image['filter']}  {first_image['object']}  "
                                           f"{first_image['dateobs'][:10]}", f'{output_basename}_sky.pdf')

    # Measurement table and summary for further analysis
    output_table = Table({name: values for name, values in table.items()})
    output_table['deviation'] = deviation
    output_table['deviation_error'] = deviation_error
    output_table.write(f'{output_basename}_measurements.ecsv', overwrite=True)
    with open(f'{output_basename}_summary.json', 'w') as summary_file:
        json.dump(results, summary_file, indent=1, default=float)

    # Report
    print(f"\n{first_image['instrument']} {first_image['filter']} {first_image['object']}: {results['number_of_stars']} stars, "
          f"{results['number_of_fitted_measurements']} measurements fitted (aperture r = {args.main_aperture:g} px)")
    for image in results['images']:
        questionable_note = f"  QUESTIONABLE (clouds? {', '.join(image['questionable_reasons'])})" if image['questionable'] else ''
        print(f"  {image['file']}  t = {image['exptime']:7.2f} s  airmass {image['airmass']:.3f}  "
              f"FWHM {image['fwhm_arcsec']:.2f}\"  sky {image['sky_level']:8.1f} e-/px ({image['sky_rate']:.3f} e-/px/s, "
              f"noise {image['sky_noise']:.1f} e-)  Z = {image['throughput']:.4f}, "
              f"Z_quality = {image['throughput_quality']:.4f} ({image['throughput_quality_stars']} stars){questionable_note}")
    print(f"  airmass law check, all images (k = {extinction_coefficient:.2f} mag/airmass, {extinction_source}; "
          f"deficit vs the clear images, rejected above the threshold):")
    for image in results['images']:
        verdict = 'NOT explained by airmass - clouds?' if 'airmass law' in image['questionable_reasons'] else 'consistent with airmass'
        print(f"    {image['file'][-20:-12]}  t = {image['exptime']:6.1f} s  X = {image['airmass']:.3f}  deficit "
              f"{image['airmass_law_deficit_percent']:+.2f}% (threshold {image['airmass_law_threshold_percent']:.2f}%) -> {verdict}")
    if repeated_exposure_check:
        print(f"  repeated exposure time check (k = {extinction_coefficient:.2f} mag/airmass, {extinction_source}):")
        for group in repeated_exposure_check:
            deficits = ', '.join(f"{results['images'][image]['file'][-20:-12]} {deficit:+.2f}%"
                                 for image, deficit in zip(group['images'], group['deficit_percent']))
            verdict = 'consistent with airmass' if group['is_consistent_with_airmass'] else 'NOT explained by airmass - clouds?'
            print(f"    t = {group['exptime']:.1f} s: deficit vs clearest image: {deficits} -> {verdict}")
    else:
        print("  repeated exposure time check: no exposure time was repeated, not tested")
    if is_questionable_image.any():
        action = 'excluded from the fit' if args.exclude_questionable else 'kept in the fit (use --exclude-questionable)'
        print(f"  WARNING: {is_questionable_image.sum()} questionable image(s), {action}")
        is_short_exposure_flagged_by_law_only = (airmass_law_check['is_flagged'] & ~is_flagged_by_repeated_exposures
                                                 & (exposure_times < np.median(exposure_times)))
        if (residual_shutter_test is not None and residual_shutter_test['fixed_k']['is_significant']
                and is_short_exposure_flagged_by_law_only.any()):
            print("  WARNING: the residual shutter offset is significant; airmass law flags on short exposures may be "
                  "caused by the shutter rather than by clouds")
    throughput_models = results['throughput_models']
    if throughput_models is not None:
        print(f"  throughput explained by a shutter offset of {throughput_models['delta_t'] * 1000:+.0f} ms "
              f"(rms {throughput_models['shutter_rms'] * 100:.2f}%) or by extinction k = "
              f"{throughput_models['extinction_k']:.2f} mag/airmass (rms {throughput_models['extinction_rms'] * 100:.2f}%)")
        if throughput_models['is_degenerate']:
            print(f"  WARNING: airmass and 1/exposure time are correlated ({throughput_models['airmass_exptime_correlation']:+.2f}); "
                  f"shutter offset and extinction cannot be separated. Shuffle the exposure time order.")
    if residual_shutter_test is None:
        print("  residual shutter offset test: fewer than 3 distinct exposure times, not tested")
    else:
        fixed_k = residual_shutter_test['fixed_k']
        verdict = 'SIGNIFICANT' if fixed_k['is_significant'] else 'not significant'
        print(f"  residual shutter offset after extinction correction (k = {fixed_k['extinction_coefficient']:.2f}, "
              f"{extinction_source}): Δt = {fixed_k['delta_t'] * 1000:+.0f} ± {fixed_k['delta_t_error'] * 1000:.0f} ms "
              f"({fixed_k['significance']:.1f} σ, {verdict}; {fixed_k['effect_on_shortest_exposure_percent']:+.2f}% at the "
              f"shortest exposure); reduced χ² {fixed_k['reduced_chi2_without_shutter']:.2f} without, "
              f"{fixed_k['reduced_chi2_with_shutter']:.2f} with shutter term")
        joint = residual_shutter_test['joint']
        if joint is not None:
            verdict = 'SIGNIFICANT' if joint['is_significant'] else 'not significant'
            print(f"  joint fit: k = {joint['extinction_coefficient']:.2f} ± {joint['extinction_coefficient_error']:.2f} mag/airmass, "
                  f"Δt = {joint['delta_t'] * 1000:+.0f} ± {joint['delta_t_error'] * 1000:.0f} ms ({joint['significance']:.1f} σ, "
                  f"{verdict}), correlation(k, Δt) = {joint['k_delta_t_correlation']:+.2f}, {joint['degrees_of_freedom']} dof")
            if joint['is_degenerate']:
                print("  WARNING: k and Δt are strongly correlated in the joint fit; the image sequence cannot separate "
                      "extinction from a shutter offset. Shuffle the exposure time order or give --extinction-coefficient.")
    print(f"  deviation from linearity (relative to {args.reference_peak_range[0]:.0f}-{args.reference_peak_range[1]:.0f} e- peak):")
    for level, value in results['deviation_percent_at_level'].items():
        print(f"    peak {int(level):6d} e-: {value['value']:+.2f} ± {value['error']:.2f} %")
    if sky_scaling_test is None:
        print("  sky scaling test: fewer than 3 images or 2 distinct exposure times, not tested")
    else:
        models = sky_scaling_test['models']
        offset = models['offset']['parameters']
        print(f"  sky scaling with exposure time: S = r t + b with r = {offset['r']['value']:.3f} e-/px/s, "
              f"b = {offset['b']['value']:+.1f} ± {offset['b']['error']:.1f} e- ({offset['b']['significance']:.1f} σ, "
              f"{'SIGNIFICANT additive offset' if sky_scaling_test['is_offset_significant'] else 'not significant'}); "
              f"rms residual {models['offset']['rms_residual']:.1f} e- (proportional model: {models['proportional']['rms_residual']:.1f} e-)")
        if models['offset_drift'] is not None:
            drift = models['offset_drift']['parameters']
            print(f"    with sky drift: r0 = {drift['r0']['value']:.3f} e-/px/s, drift r1 = {drift['r1']['value']:+.3f} ± "
                  f"{drift['r1']['error']:.3f} e-/px/s per hour, b = {drift['b']['value']:+.1f} ± {drift['b']['error']:.1f} e-, "
                  f"rms residual {models['offset_drift']['rms_residual']:.1f} e-")
        if models['curvature'] is not None:
            curvature = models['curvature']['parameters']['q']
            print(f"    curvature q = {curvature['value']:+.2e} ± {curvature['error']:.1e} e-/s^2 ({curvature['significance']:.1f} σ, "
                  f"{'SIGNIFICANT: sky not linear in exposure time' if sky_scaling_test['is_curvature_significant'] else 'not significant'})")
    aperture = results['aperture_analysis']
    print(f"  aperture analysis: optimum r = {aperture['optimum_radius_px']:g} px, collects "
          f"{aperture['encircled_energy_at_optimum'] * 100:.1f}% of the light within {GROWTH_CURVE_RADII[-1]:g} px; "
          f"missed light uncertainty {aperture['missed_light_uncertainty_at_optimum'] * 100:.2f}% "
          f"(variability {aperture['missed_light_variability']:.2f} of the missed light), "
          f"mean neighbour light {aperture['mean_contamination_at_optimum'] * 100:.2f}%, "
          f"total {aperture['total_error_at_optimum'] * 100:.2f}%")
    print(f"    99% of the light within r = {aperture['radius_99_percent_px']:g} px; recommended sky annulus "
          f"{aperture['annulus_px'][0]:g}-{aperture['annulus_px'][1]:g} px")
    if not args.auto_aperture and (args.main_aperture != aperture['optimum_radius_px']
                                   or list(args.annulus) != aperture['annulus_px']):
        print(f"  note: aperture r = {args.main_aperture:g} px and annulus {args.annulus[0]:g}-{args.annulus[1]:g} px were used; "
              f"--auto-aperture applies the optimum")
    print(f"  plots and tables: {output_basename}*")


if __name__ == '__main__':
    main()
