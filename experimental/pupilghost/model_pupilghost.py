"""
Empirical pupil-ghost model for twilight sky flats (developed on elp 1m0-08 / ep60).

Model
-----
A sky flat is  F_f(x,y) = P_f(x,y) * (1 + G_f(r))  with P_f the true response (pixel QE x illumination)
and G_f the additive pupil ghost, expressed as a fraction of the sky level.  G_f is assumed
axisymmetric about a per-filter ghost centre.

  1. Bin every master flat (BPM-masked median), ratio each against a ghost-free reference
     filter (default rp). This cancels pixel-to-pixel response, dust donuts and vignetting.
  2. Divide out a 2D quadratic (twilight gradient, chromatic illumination) fitted in the
     normalisation zone (an annulus around the detector centre), which defines G = 0 there.
     -> G is only defined up to a constant; a uniform ghost floor is degenerate with the flat
     normalisation and does not affect relative photometry.
  3. Fit the ghost centre (maximise azimuthal symmetry) and a least-squares cubic spline G_f(r),
     forced to zero beyond the inner radius of the normalisation zone.

Non-axisymmetric chromatic structure (QE patches, blue stripe pattern) stays in the residual and
deliberately remains in the flat field.

Correcting a flat:   F_corr = F / (1 + G_f(r));  photometric error of an uncorrected flat at radius r is ~ -G_f(r).

Usage:
    python model_pupilghost.py --flat-dir flatcache --outdir plots/pupilghost [--full-resolution]
"""
import argparse
import glob
import json
import logging
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.interpolate import make_lsq_spline
from scipy.optimize import minimize

_log = logging.getLogger(__name__)

# Flats are binned BINNING x BINNING pixels (median) before modelling; the ghost varies on arcminute scales
BINNING = 16

# Half width [binned pixels] of the central box used to normalise each binned flat to 1
NORMALISATION_BOX_HALF_WIDTH = 16

# The ghost centre may not wander further than this [binned pixels] from the detector centre. This keeps the fit
# on the chip for filters without a detectable ghost.
MAX_CENTRE_OFFSET_BINNED = 40

# Spline knot spacing of the radial ghost profile [arcmin]
SPLINE_KNOT_SPACING_ARCMIN = 0.5


# ----------------------------------------------------------------------------------------------------------------------
# Binned flat fields
# ----------------------------------------------------------------------------------------------------------------------
def bin_flat(filename, binning=BINNING):
    """ Read a master flat, mask bad pixels, median-bin it and normalise it to 1 in the detector centre."""
    with fits.open(filename) as hdulist:
        flat = hdulist['SCI'].data.astype(np.float32)
        flat[hdulist['BPM'].data > 0] = np.nan
        header = hdulist['SCI'].header

    # Median of each binning x binning block, ignoring bad pixels
    binned_size = flat.shape[0] // binning
    pixel_blocks = flat.reshape(binned_size, binning, binned_size, binning).transpose(0, 2, 1, 3)
    binned_flat = np.nanmedian(pixel_blocks.reshape(binned_size, binned_size, binning * binning), axis=2)

    centre = binned_size // 2
    central_box = binned_flat[centre - NORMALISATION_BOX_HALF_WIDTH:centre + NORMALISATION_BOX_HALF_WIDTH,
                              centre - NORMALISATION_BOX_HALF_WIDTH:centre + NORMALISATION_BOX_HALF_WIDTH]
    return binned_flat / np.nanmedian(central_box), header


def load_binned_flats(flat_dir, cache_filename):
    """ Binned versions of all full-frame (bin1x1) master flats in flat_dir, cached in an .npz file.

    :return: dict flat name -> binned flat, dict flat name -> metadata (filter, pixel scale, full image shape)
    """
    if os.path.exists(cache_filename):
        cache = np.load(cache_filename)
        binned_flats = {name: cache[name] for name in cache.files if name != '__meta__'}
        return binned_flats, json.loads(str(cache['__meta__']))

    binned_flats, flat_metadata = {}, {}
    for filename in sorted(glob.glob(os.path.join(flat_dir, '*bin1x1-*.fits.fz'))):
        flat_name = os.path.basename(filename).replace('.fits.fz', '')
        binned_flats[flat_name], header = bin_flat(filename)
        flat_metadata[flat_name] = dict(filter=header['FILTER'], pixscale=header['PIXSCALE'],
                                        shape=[header['NAXIS2'], header['NAXIS1']])
        _log.info(f"binned {flat_name}")
    np.savez(cache_filename, __meta__=json.dumps(flat_metadata), **binned_flats)
    return binned_flats, flat_metadata


# ----------------------------------------------------------------------------------------------------------------------
# Ghost model of one filter
# ----------------------------------------------------------------------------------------------------------------------
def quadratic_surface_design(x_grid, y_grid, centre_x, centre_y, image_size):
    """ Design matrix of a 2D quadratic surface (1, x, y, x^2, xy, y^2) in normalised coordinates."""
    x_normalised = ((x_grid - centre_x) / (image_size / 2)).ravel()
    y_normalised = ((y_grid - centre_y) / (image_size / 2)).ravel()
    return np.c_[np.ones_like(x_normalised), x_normalised, y_normalised,
                 x_normalised * x_normalised, x_normalised * y_normalised, y_normalised * y_normalised]


def fit_radial_spline(radius_arcmin, ghost_fraction, max_radius_arcmin, knot_spacing=SPLINE_KNOT_SPACING_ARCMIN):
    """ Least-squares cubic spline G(r) through the data inside max_radius_arcmin, pinned to zero at max_radius_arcmin."""
    is_valid = np.isfinite(ghost_fraction) & (radius_arcmin < max_radius_arcmin)
    sort_order = np.argsort(radius_arcmin[is_valid])
    sorted_radius = radius_arcmin[is_valid][sort_order]
    sorted_ghost = ghost_fraction[is_valid][sort_order]

    # Pin the profile to zero at max_radius_arcmin by appending zero-valued points just outside
    number_of_pinning_points = 50
    sorted_radius = np.r_[sorted_radius, np.linspace(max_radius_arcmin, max_radius_arcmin + 1, number_of_pinning_points)]
    sorted_ghost = np.r_[sorted_ghost, np.zeros(number_of_pinning_points)]

    # Knots every knot_spacing; boundary knots repeated 4 times as required for a cubic B-spline
    knots = np.r_[[0] * 4, np.arange(knot_spacing, max_radius_arcmin + 1, knot_spacing), [max_radius_arcmin + 1] * 4]
    return make_lsq_spline(sorted_radius, sorted_ghost, knots, k=3)


def fit_ghost_centre(ghost_image, x_grid, y_grid, arcmin_per_pixel, max_radius_arcmin, start_centre):
    """ Find the centre about which the ghost image is most azimuthally symmetric.

    For a trial centre, a radial spline is fitted; the centre with the smallest sum of squared residuals around that
    spline wins. Five starting points (the start and +/-15 pixel offsets) guard against local minima.
    """

    def sum_of_squared_residuals(trial_centre):
        if np.hypot(trial_centre[0] - start_centre[0], trial_centre[1] - start_centre[1]) > MAX_CENTRE_OFFSET_BINNED:
            return np.inf
        radius_arcmin = np.hypot(x_grid - trial_centre[0], y_grid - trial_centre[1]) * arcmin_per_pixel
        is_valid = np.isfinite(ghost_image) & (radius_arcmin < max_radius_arcmin)
        radial_spline = fit_radial_spline(radius_arcmin[is_valid], ghost_image[is_valid], max_radius_arcmin)
        return np.sum((ghost_image[is_valid] - radial_spline(radius_arcmin[is_valid])) ** 2)

    starting_points = [start_centre,
                       (start_centre[0] - 15, start_centre[1]), (start_centre[0] + 15, start_centre[1]),
                       (start_centre[0], start_centre[1] - 15), (start_centre[0], start_centre[1] + 15)]
    best_result = None
    for starting_point in starting_points:
        result = minimize(sum_of_squared_residuals, starting_point, method='Nelder-Mead', options=dict(xatol=0.05, fatol=1e-12))
        if best_result is None or result.fun < best_result.fun:
            best_result = result
    return best_result.x


def model_filter(binned_flats, reference_flat, arcmin_per_pixel, normalisation_zone):
    """ Ghost model of one filter from its binned flats (one or more nights) and the reference filter's flat.

    :return: dict with the averaged ratio image (= ghost image), per-flat ratios, ghost centre [binned px],
             radial spline, 2D model, residual and the radius map [arcmin]
    """
    image_size = reference_flat.shape[0]
    y_grid, x_grid = np.mgrid[:image_size, :image_size].astype(float)
    detector_centre = ((image_size - 1) / 2, (image_size - 1) / 2)
    radius_from_detector_centre = np.hypot(x_grid - detector_centre[0], y_grid - detector_centre[1]) * arcmin_per_pixel

    # Quadratic surface fitted only in the normalisation zone, where the ghost is defined to be zero
    surface_design = quadratic_surface_design(x_grid, y_grid, *detector_centre, image_size)
    is_in_normalisation_zone = ((radius_from_detector_centre.ravel() > normalisation_zone[0])
                                & (radius_from_detector_centre.ravel() < normalisation_zone[1]))

    # Ratio of each flat to the reference, divided by its own quadratic surface (twilight gradient differs per night)
    ghost_images = []
    for binned_flat in binned_flats:
        ratio = (binned_flat / reference_flat).ravel()
        is_fit_pixel = is_in_normalisation_zone & np.isfinite(ratio)
        surface_coefficients = np.linalg.lstsq(surface_design[is_fit_pixel], ratio[is_fit_pixel], rcond=None)[0]
        ghost_images.append((ratio / (surface_design @ surface_coefficients)).reshape(image_size, image_size) - 1)
    mean_ghost_image = np.nanmean(ghost_images, axis=0)

    # Ghost centre and radial profile G(r); the model is zero outside the inner radius of the normalisation zone
    ghost_centre_x, ghost_centre_y = fit_ghost_centre(mean_ghost_image, x_grid, y_grid, arcmin_per_pixel,
                                                      normalisation_zone[0], detector_centre)
    radius_from_ghost_centre = np.hypot(x_grid - ghost_centre_x, y_grid - ghost_centre_y) * arcmin_per_pixel
    radial_spline = fit_radial_spline(radius_from_ghost_centre.ravel(), mean_ghost_image.ravel(), normalisation_zone[0])
    ghost_model_image = evaluate_ghost_profile(radial_spline, radius_from_ghost_centre, normalisation_zone[0])

    return dict(ratio=mean_ghost_image, ratios=ghost_images, cx=ghost_centre_x, cy=ghost_centre_y,
                spline=radial_spline, model=ghost_model_image, resid=mean_ghost_image - ghost_model_image,
                r=radius_from_ghost_centre)


def evaluate_ghost_profile(radial_spline, radius_arcmin, max_radius_arcmin):
    """ G(r) from the spline inside max_radius_arcmin, zero outside."""
    return np.where(radius_arcmin < max_radius_arcmin, radial_spline(np.minimum(radius_arcmin, max_radius_arcmin)), 0)


# ----------------------------------------------------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------------------------------------------------
def write_full_resolution_model(filter_name, radial_spline, centre_x_pixel, centre_y_pixel, full_image_shape,
                                pixel_scale_arcsec, max_radius_arcmin, reference_filter, outdir):
    """ Write 1 + G at full detector resolution; a flat is corrected by dividing it by this image."""
    y_grid, x_grid = np.mgrid[:full_image_shape[0], :full_image_shape[1]]
    radius_arcmin = np.hypot(x_grid - centre_x_pixel, y_grid - centre_y_pixel) * pixel_scale_arcsec / 60
    ghost_fraction = evaluate_ghost_profile(radial_spline, radius_arcmin, max_radius_arcmin).astype(np.float32)
    header = fits.Header(dict(FILTER=filter_name, GHOSTCX=centre_x_pixel, GHOSTCY=centre_y_pixel, REFFILT=reference_filter,
                              COMMENT='1 + pupil ghost fraction; corrected flat = flat / this'))
    fits.CompImageHDU(1 + ghost_fraction, header=header).writeto(
        os.path.join(outdir, f'pupilghost_{filter_name}.fits.fz'), overwrite=True)


def plot_radial_profiles(filter_models, args):
    """ Azimuthal median of each filter's ghost image (dots) and the spline model (line)."""
    fig, axis = plt.subplots(figsize=(10, 6))
    radius_grid = np.arange(0, 17, 0.1)
    profile_bin_edges = np.arange(0, 17.01, 0.25)
    for filter_name, filter_model in filter_models.items():
        radius_arcmin = filter_model['r'].ravel()
        ghost_image = filter_model['ratio'].ravel()
        azimuthal_median = [np.nanmedian(ghost_image[(radius_arcmin >= inner) & (radius_arcmin < outer)])
                            for inner, outer in zip(profile_bin_edges[:-1], profile_bin_edges[1:])]
        data_line, = axis.plot(profile_bin_edges[:-1] + 0.125, np.array(azimuthal_median) * 100, '.', ms=3)
        model_profile = evaluate_ghost_profile(filter_model['spline'], radius_grid, args.normalisation_zone[0])
        axis.plot(radius_grid, model_profile * 100, '-', color=data_line.get_color(), label=filter_name)
    axis.axvspan(*args.normalisation_zone, color='grey', alpha=0.15, label='normalisation zone')
    axis.axhline(0, color='k', lw=0.5)
    axis.set_xlabel('radius from ghost centre [arcmin]')
    axis.set_ylabel(f'pupil ghost G (flat/{args.reference} - 1) [%]')
    axis.grid()
    axis.legend(ncol=2)
    plt.savefig(os.path.join(args.outdir, 'pupilghost_profiles.png'), dpi=100, bbox_inches='tight')
    plt.close()


def plot_model_images(filter_models, args):
    """ Per filter: ratio to the reference (data), radial ghost model, and residual."""
    rows = [('ratio', 'flat/ref', 2), ('model', 'ghost model', 2), ('resid', 'residual', 0.5)]  # key, title, colour limit [%]
    fig, axes = plt.subplots(len(rows), len(filter_models), figsize=(3.2 * len(filter_models), 9.5), squeeze=False)
    for column, (filter_name, filter_model) in enumerate(filter_models.items()):
        for row, (image_key, title, colour_limit_percent) in enumerate(rows):
            axis = axes[row, column]
            image = axis.imshow(filter_model[image_key] * 100, origin='lower', cmap='RdBu_r',
                                vmin=-colour_limit_percent, vmax=colour_limit_percent)
            axis.set_title(f'{filter_name} {title} [%]', fontsize=9)
            axis.set_xticks([])
            axis.set_yticks([])
            if column == len(filter_models) - 1:
                plt.colorbar(image, ax=axis, fraction=0.046)
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, 'pupilghost_2d.png'), dpi=80)
    plt.close()


# ----------------------------------------------------------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------------------------------------------------------
def parseCommandLine():
    parser = argparse.ArgumentParser(description='Empirical pupil-ghost model from twilight master flats',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--flat-dir', default='flatcache', help='directory with BANZAI master sky flats')
    parser.add_argument('--outdir', default='plots/pupilghost')
    parser.add_argument('--reference', default='rp', help='ghost-free reference filter')
    parser.add_argument('--filters', nargs='+', default=['V', 'R', 'I', 'ip', 'zs', 'Y', 'w', 'gp', 'B'],
                        help='filters to model')
    parser.add_argument('--normalisation-zone', nargs=2, type=float, default=[12.5, 16.5],
                        help='radial range [arcmin] defined as ghost-free; G = 0 there')
    parser.add_argument('--full-resolution', action='store_true', help='also write full-resolution 1+G FITS per filter')
    parser.add_argument('--loglevel', default='INFO', choices=['DEBUG', 'INFO', 'WARN'])
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.loglevel), format='%(asctime)s %(levelname)s %(message)s')
    return args


def main():
    args = parseCommandLine()
    os.makedirs(args.outdir, exist_ok=True)

    # Binned flats, grouped by filter (several nights per filter are averaged later)
    binned_flats, flat_metadata = load_binned_flats(args.flat_dir, os.path.join(args.outdir, f'binned{BINNING}.npz'))
    binned_flats_by_filter = {}
    for flat_name, binned_flat in binned_flats.items():
        binned_flats_by_filter.setdefault(flat_metadata[flat_name]['filter'], []).append(binned_flat)

    any_flat_metadata = next(iter(flat_metadata.values()))
    pixel_scale_arcsec = any_flat_metadata['pixscale']
    full_image_shape = any_flat_metadata['shape']
    arcmin_per_binned_pixel = BINNING * pixel_scale_arcsec / 60
    reference_flat = np.nanmean(binned_flats_by_filter[args.reference], axis=0)
    max_ghost_radius_arcmin = args.normalisation_zone[0]

    filter_models, summary = {}, {}
    for filter_name in args.filters:
        if filter_name not in binned_flats_by_filter:
            _log.warning(f"no flats for filter {filter_name}")
            continue
        filter_model = model_filter(binned_flats_by_filter[filter_name], reference_flat, arcmin_per_binned_pixel,
                                    args.normalisation_zone)
        filter_models[filter_name] = filter_model

        # Ghost centre in unbinned pixel coordinates (0-based): centre of the binned pixel
        centre_x_pixel = float(filter_model['cx'] * BINNING + (BINNING - 1) / 2)
        centre_y_pixel = float(filter_model['cy'] * BINNING + (BINNING - 1) / 2)

        # Key numbers of the radial profile
        profile_radius = np.arange(0, max_ghost_radius_arcmin + 0.01, 0.1)
        profile_ghost = filter_model['spline'](profile_radius)
        is_inside_ghost = filter_model['r'] < max_ghost_radius_arcmin
        summary[filter_name] = dict(
            center_x_px=centre_x_pixel, center_y_px=centre_y_pixel,
            peak_to_peak_pct=float((profile_ghost.max() - profile_ghost.min()) * 100),
            min_pct=float(profile_ghost.min() * 100), r_min_arcmin=float(profile_radius[profile_ghost.argmin()]),
            max_pct=float(profile_ghost.max() * 100), r_max_arcmin=float(profile_radius[profile_ghost.argmax()]),
            resid_rms_pct=float(np.nanstd(filter_model['resid'][is_inside_ghost]) * 100),
            nflats=len(binned_flats_by_filter[filter_name]))
        _log.info(f"{filter_name}: {summary[filter_name]}")

        # Radial profile table, read by starflat_photometry.py
        np.savetxt(os.path.join(args.outdir, f'ghost_profile_{filter_name}.dat'), np.c_[profile_radius, profile_ghost],
                   header=f"r[arcmin] G (fraction of sky); centre x,y [unbinned px, 0-based] = "
                          f"{centre_x_pixel:.1f} {centre_y_pixel:.1f}")

        if args.full_resolution:
            write_full_resolution_model(filter_name, filter_model['spline'], centre_x_pixel, centre_y_pixel,
                                        full_image_shape, pixel_scale_arcsec, max_ghost_radius_arcmin, args.reference,
                                        args.outdir)

    with open(os.path.join(args.outdir, 'pupilghost_summary.json'), 'w') as summary_file:
        json.dump(summary, summary_file, indent=1)

    plot_radial_profiles(filter_models, args)
    plot_model_images(filter_models, args)


if __name__ == '__main__':
    main()
