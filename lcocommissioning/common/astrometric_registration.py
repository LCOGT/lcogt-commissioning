"""
Astrometric registration of frames without a valid WCS (e.g. BANZAI frames with WCSERR != 0) against a reference star
list: the pointing offset from a 2D histogram of star-to-star offsets, then a linear WCS (shift, rotation, scale) fitted
to the matched stars.

Used by onsky_linearity.py and experimental/pupilghost/starflat_photometry.py.
"""
import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.wcs.utils import fit_wcs_from_points

REGISTRATION_SEARCH_RADIUS_ARCSEC = 60.  # maximum expected pointing error of the nominal WCS
REGISTRATION_BIN_ARCSEC = 1.0            # bin size of the offset histogram
REGISTRATION_MIN_MATCHES = 10            # minimum number of matched stars for a WCS fit


def frame_has_valid_astrometry(filename, header):
    """ True if BANZAI solved the astrometry: WCSERR == 0, and the catalogue (if present) has ra / dec columns."""
    if header.get('WCSERR', 0) != 0:
        return False
    with fits.open(filename) as hdulist:
        if 'CAT' not in hdulist:
            return True
        catalogue_columns = [name.lower() for name in hdulist['CAT'].columns.names]
    return 'ra' in catalogue_columns and 'dec' in catalogue_columns


def find_pointing_offset(frame_positions, reference_positions, search_radius_arcsec=REGISTRATION_SEARCH_RADIUS_ARCSEC,
                         bin_size_arcsec=REGISTRATION_BIN_ARCSEC):
    """ Most common offset (east, north [arcsec]) from the frame's stars to the reference stars.

    All pairs within the search radius vote in a 2D histogram of offsets; the true pointing error collects the votes of
    all real matches in one bin, random pairs spread evenly. The peak is refined with the median of the pairs near it.
    """
    frame_index, reference_index, _, _ = reference_positions.search_around_sky(frame_positions, search_radius_arcsec * u.arcsec)
    if len(frame_index) == 0:
        return None
    east_offset, north_offset = frame_positions[frame_index].spherical_offsets_to(reference_positions[reference_index])
    east_offset = east_offset.to_value(u.arcsec)
    north_offset = north_offset.to_value(u.arcsec)
    bin_edges = np.arange(-search_radius_arcsec, search_radius_arcsec + bin_size_arcsec, bin_size_arcsec)
    votes, _, _ = np.histogram2d(east_offset, north_offset, bins=[bin_edges, bin_edges])
    peak_east_bin, peak_north_bin = np.unravel_index(np.argmax(votes), votes.shape)
    peak_east = 0.5 * (bin_edges[peak_east_bin] + bin_edges[peak_east_bin + 1])
    peak_north = 0.5 * (bin_edges[peak_north_bin] + bin_edges[peak_north_bin + 1])
    is_near_peak = np.hypot(east_offset - peak_east, north_offset - peak_north) < 2 * bin_size_arcsec
    return float(np.median(east_offset[is_near_peak])), float(np.median(north_offset[is_near_peak])), int(votes.max())


def register_frame(x_pixel, y_pixel, nominal_wcs, reference_positions):
    """ WCS of a frame from its star pixel positions and a reference star list.

    Start: nominal WCS shifted by the pointing offset. Then a linear WCS (shift, rotation, scale) is fitted to the
    stars matched within 3", and refitted to the stars matched within 1" with the improved WCS.

    :return: dict with the fitted WCS, number of matched stars, rms match residual ["] and the pointing offset, or None
    """
    frame_positions = nominal_wcs.pixel_to_world(x_pixel, y_pixel)
    pointing_offset = find_pointing_offset(frame_positions, reference_positions)
    if pointing_offset is None:
        return None
    east_offset, north_offset, offset_votes = pointing_offset
    corrected_positions = frame_positions.spherical_offsets_by(east_offset * u.arcsec, north_offset * u.arcsec)

    fitted_wcs = None
    for match_radius_arcsec in (3.0, 1.0):
        reference_index, separation, _ = corrected_positions.match_to_catalog_sky(reference_positions)
        is_matched = separation.arcsec < match_radius_arcsec
        if is_matched.sum() < REGISTRATION_MIN_MATCHES:
            return None
        fitted_wcs = fit_wcs_from_points((x_pixel[is_matched], y_pixel[is_matched]), reference_positions[reference_index[is_matched]],
                                         proj_point='center', projection='TAN')
        corrected_positions = fitted_wcs.pixel_to_world(x_pixel, y_pixel)

    reference_index, separation, _ = corrected_positions.match_to_catalog_sky(reference_positions)
    is_matched = separation.arcsec < 1.0
    return dict(wcs=fitted_wcs, number_of_matches=int(is_matched.sum()),
                rms_arcsec=float(np.sqrt(np.mean(separation.arcsec[is_matched] ** 2))),
                pointing_offset_arcsec=(east_offset, north_offset), offset_votes=offset_votes)
