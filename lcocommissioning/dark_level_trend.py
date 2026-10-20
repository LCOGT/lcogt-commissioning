"""
Trend of the dark and bias offset levels of a camera over the last days, correlated with time, CCD temperature and
ambient conditions, optionally split into epochs (e.g. before / after a change of the camera readout).

Darks: for every processed (BANZAI, RLEVEL 91) individual dark frame in the selected days
  - header data from the archive OpenSearch index: DATE-OBS, EXPTIME, CCDATEMP, WMSTEMP (weather station air
    temperature), TUBETEMP (telescope tube), FOCTEMP (focuser), BIASLVL (bias level removed by BANZAI, i.e. the level of
    the master bias used), SATURATE, enclosure state and lights
  - dark level from the image: sigma-clipped median of the frame (100 pixel border excluded) and of a 3 x 3 grid of
    regions, in e-/s (BANZAI stores processed darks normalised to one second)
  - derived: local solar time and altitude of the Sun at the site
Biases: the bias level (BIASLVL) of every processed individual bias, from the archive index only (no download).

Bias offset model, per epoch: without overscan, BANZAI removes the bias level by subtracting the master bias. If the
bias level drifts with temperature, every frame keeps an additive offset
    bias offset = bias level at the time of the frame - bias level of the master bias used.
The bias level is modelled as BIASLVL = a + b * TUBETEMP from the daily biases; the predicted offset of each dark
(divided by its exposure time) is compared with the measured dark level, and removed:
    measured dark level = alpha * predicted offset + dark current.

Frames are downloaded one at a time (in parallel threads), measured and deleted; the measurements are kept in a cache
table so that a rerun only processes new frames.

Output: correlation and bias-model summary on stdout, a table (.ecsv) and two figures (.pdf): the dark level trend and
the bias offset analysis.

Usage:
    export ARCHIVE_API_TOKEN=xxxxxxxx
    darktrend --camera ep60 --ndays 11 --epoch-boundary 2026-10-19T22:00 --outdir plots/darktrend
"""
import argparse
import datetime as dt
import logging
import os
from concurrent.futures import ThreadPoolExecutor

import matplotlib

matplotlib.use('Agg')
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, get_sun
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, join, vstack
from astropy.time import Time
from scipy import stats

import lcocommissioning.common.lco_archive_utilities as lco_archive_utilities
from lcocommissioning.common.common import lco_site_alt, lco_site_lonlat

_log = logging.getLogger(__name__)

# Header keywords taken from the archive index for every frame
HEADER_KEYWORDS = ['DATE-OBS', 'MJD-OBS', 'EXPTIME', 'CCDATEMP', 'CCDSTEMP', 'WMSTEMP', 'TUBETEMP', 'FOCTEMP',
                   'BIASLVL', 'L1IDBIAS', 'SATURATE', 'ENC1STAT', 'ENCWLIGT', 'ENCRLIGT', 'SITEID']

# Header keywords holding text; all others are numbers (missing values become NaN)
TEXT_KEYWORDS = {'DATE-OBS', 'L1IDBIAS', 'ENC1STAT', 'ENCWLIGT', 'ENCRLIGT', 'SITEID'}

# Columns of the measurement cache; header data are always taken fresh from the archive index
MEASUREMENT_COLUMNS = ['filename', 'dark_level', 'dark_noise', 'region_levels']

# Quantities the dark level is correlated with: (column, axis label)
CORRELATED_QUANTITIES = [('days_since_start', 'time since first dark [days]'),
                         ('local_solar_hour', 'local solar time [h]'),
                         ('sun_altitude', 'altitude of the Sun [deg]'),
                         ('CCDATEMP', 'CCD temperature CCDATEMP [C]'),
                         ('WMSTEMP', 'outside air temperature WMSTEMP [C]'),
                         ('TUBETEMP', 'telescope tube temperature TUBETEMP [C]'),
                         ('FOCTEMP', 'focuser temperature FOCTEMP [C]')]

# Border [pixels] excluded from the dark level (edge effects), and the grid of regions for the spatial pattern
BORDER = 100
GRID_SIZE = 3

# Plot style per epoch
EPOCH_COLOURS = ['C0', 'C3', 'C2', 'C1', 'C4']
EPOCH_MARKERS = ['o', '^', 's', 'D', 'v']


# ----------------------------------------------------------------------------------------------------------------------
# Archive
# ----------------------------------------------------------------------------------------------------------------------
def query_frames(camera, dates, readmode, obstype, frame_code):
    """ Archive index records of the individual processed frames of one OBSTYPE (frame_code e.g. 'd91', 'b91').

    Master frames have no sequence number and therefore no '-<frame_code>' in their name; they are skipped.
    """
    rows = []
    for dayobs in dates:
        search = lco_archive_utilities.make_opensearch(
            'fitsheaders', [{'DAY-OBS': dayobs}, {'RLEVEL': 91}, {'INSTRUME': camera}], None,
            terms_filters=[{'OBSTYPE': [obstype]}, {'CONFMODE': [readmode]}])
        for hit in lco_archive_utilities.scan_opensearch(search, f'{camera} {obstype} {dayobs}'):
            record = hit.to_dict()
            if f'-{frame_code}' not in record['filename']:
                continue
            row = {'filename': record['filename'], 'frameid': int(record['frameid']), 'DAY-OBS': dayobs}
            for keyword in HEADER_KEYWORDS:
                value = record.get(keyword)
                if keyword in TEXT_KEYWORDS:
                    row[keyword] = '' if value is None else str(value)
                else:
                    row[keyword] = np.nan if value in (None, '', 'N/A', 'UNKNOWN') else float(value)
            rows.append(row)
    _log.info(f"{len(rows)} {obstype} frames of {camera} in {dates[0]} - {dates[-1]}")
    return rows


# ----------------------------------------------------------------------------------------------------------------------
# Dark level measurement
# ----------------------------------------------------------------------------------------------------------------------
def measure_dark_level(frameid):
    """ Download one dark, return its sigma-clipped median level [e-/s] and the median of each grid region."""
    hdulist = lco_archive_utilities.download_from_archive(frameid)
    try:
        image = hdulist['SCI'].data.astype(np.float32)[BORDER:-BORDER, BORDER:-BORDER]
    finally:
        hdulist.close()
    _, median_level, robust_noise = sigma_clipped_stats(image[::2, ::2], sigma=3, maxiters=5)

    region_rows = np.array_split(np.arange(image.shape[0]), GRID_SIZE)
    region_columns = np.array_split(np.arange(image.shape[1]), GRID_SIZE)
    region_levels = np.zeros((GRID_SIZE, GRID_SIZE))
    for row_index, rows in enumerate(region_rows):
        for column_index, columns in enumerate(region_columns):
            region = image[rows[0]:rows[-1] + 1:2, columns[0]:columns[-1] + 1:2]
            region_levels[row_index, column_index] = sigma_clipped_stats(region, sigma=3, maxiters=5)[1]
    return float(median_level), float(robust_noise), region_levels.ravel()


def measure_all(records, cache_filename, number_of_threads):
    """ Dark level of every record, reusing cached measurements.

    :return: table with the header data of the records joined with their measurements
    """
    cache = Table.read(cache_filename)[MEASUREMENT_COLUMNS] if os.path.exists(cache_filename) else None
    measured_filenames = set(cache['filename']) if cache is not None else set()
    to_measure = [record for record in records if record['filename'] not in measured_filenames]
    _log.info(f"{len(records) - len(to_measure)} darks already measured, {len(to_measure)} to download")

    def measure(record):
        try:
            level, noise, region_levels = measure_dark_level(record['frameid'])
        except Exception as error:
            _log.warning(f"{record['filename']}: measurement failed ({error})")
            return None
        _log.info(f"{record['filename']}: {level:.4f} e-/s")
        return dict(filename=record['filename'], dark_level=level, dark_noise=noise, region_levels=region_levels)

    with ThreadPoolExecutor(max_workers=number_of_threads) as executor:
        for row in executor.map(measure, to_measure):
            if row is not None:
                # Save after every frame, so that an interrupted run keeps its progress
                cache = vstack([cache, Table(rows=[row])]) if cache is not None else Table(rows=[row])
                cache.write(cache_filename, overwrite=True)

    if cache is None:
        return None
    return join(Table(rows=records), cache, keys='filename', join_type='inner')


# ----------------------------------------------------------------------------------------------------------------------
# Derived quantities, epochs and correlations
# ----------------------------------------------------------------------------------------------------------------------
def add_time_and_sun(table, reference_time):
    """ Mid-exposure time, days since reference_time, local solar time and the altitude of the Sun at the site."""
    site = table['SITEID'][0]
    longitude, latitude = lco_site_lonlat[site]
    location = EarthLocation(lon=longitude * u.deg, lat=latitude * u.deg, height=lco_site_alt.get(site, 2000) * u.m)
    mid_exposure = Time(np.asarray(table['DATE-OBS'], dtype=str), format='isot', scale='utc') + np.asarray(table['EXPTIME'], float) / 2 * u.s
    table['days_since_start'] = mid_exposure.mjd - reference_time.mjd
    table['local_solar_hour'] = ((mid_exposure.datetime64.astype('datetime64[s]').astype(float) / 3600 + longitude / 15) % 24)
    table['sun_altitude'] = get_sun(mid_exposure).transform_to(AltAz(obstime=mid_exposure, location=location)).alt.deg
    return table


def epoch_labels(boundaries):
    """ Human-readable label of every epoch defined by the sorted boundary times."""
    if not boundaries:
        return ['all data']
    labels = [f'before {boundaries[0].isot[:16]}']
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        labels.append(f'{start.isot[:16]} - {end.isot[:16]}')
    labels.append(f'after {boundaries[-1].isot[:16]}')
    return labels


def assign_epochs(table, boundaries):
    """ Epoch index of every frame: number of boundaries at or before its DATE-OBS."""
    frame_time = Time(np.asarray(table['DATE-OBS'], dtype=str), format='isot', scale='utc')
    table['epoch'] = np.array([int(np.sum([frame_time[index] >= boundary for boundary in boundaries]))
                               for index in range(len(table))], dtype=int)
    return table


def correlations(table):
    """ Pearson r and Spearman rho (with p-value) of the dark level with every correlated quantity."""
    results = {}
    dark_level = np.asarray(table['dark_level'], float)
    for column, label in CORRELATED_QUANTITIES:
        values = np.asarray(table[column], float)
        is_valid = np.isfinite(values) & np.isfinite(dark_level)
        if is_valid.sum() < 3 or np.ptp(values[is_valid]) == 0:
            results[column] = dict(label=label, pearson=np.nan, spearman=np.nan, spearman_p=np.nan,
                                   note='constant' if is_valid.sum() >= 3 else 'too few values')
            continue
        pearson = stats.pearsonr(values[is_valid], dark_level[is_valid])[0]
        spearman, spearman_p = stats.spearmanr(values[is_valid], dark_level[is_valid])
        results[column] = dict(label=label, pearson=float(pearson), spearman=float(spearman), spearman_p=float(spearman_p),
                               note=f'range {np.min(values[is_valid]):.2f} .. {np.max(values[is_valid]):.2f}')
    return results


# ----------------------------------------------------------------------------------------------------------------------
# Bias offset model
# ----------------------------------------------------------------------------------------------------------------------
def fit_bias_level_model(biases):
    """ BIASLVL = a + b * TUBETEMP from the daily biases of one epoch; None with fewer than 3 usable biases."""
    tube_temperature = np.asarray(biases['TUBETEMP'], float)
    bias_level = np.asarray(biases['BIASLVL'], float)
    is_valid = np.isfinite(tube_temperature) & np.isfinite(bias_level)
    if is_valid.sum() < 3 or np.ptp(tube_temperature[is_valid]) == 0:
        return None
    slope, intercept = np.polyfit(tube_temperature[is_valid], bias_level[is_valid], 1)
    residuals = bias_level[is_valid] - (intercept + slope * tube_temperature[is_valid])
    return dict(intercept=float(intercept), slope=float(slope), rms=float(np.std(residuals)),
                correlation=float(np.corrcoef(tube_temperature[is_valid], bias_level[is_valid])[0, 1]),
                number_of_biases=int(is_valid.sum()), bias_level_range=[float(bias_level[is_valid].min()),
                                                                        float(bias_level[is_valid].max())])


def apply_bias_offset_model(darks, biases, number_of_epochs):
    """ Per epoch: bias level model, predicted bias offset of every dark, and the dark level with that offset removed.

    Adds columns predicted_bias_offset [e-/s] and dark_level_corrected [e-/s] (NaN where an epoch has no bias model).
    :return: dict epoch -> bias model and dark comparison results
    """
    darks['predicted_bias_offset'] = np.full(len(darks), np.nan)
    darks['dark_level_corrected'] = np.full(len(darks), np.nan)
    results = {}
    for epoch in range(number_of_epochs):
        model = fit_bias_level_model(biases[biases['epoch'] == epoch])
        is_epoch_dark = darks['epoch'] == epoch
        results[epoch] = dict(bias_model=model, dark_comparison=None)
        if model is None or is_epoch_dark.sum() == 0:
            continue
        # offset = bias level at the dark's tube temperature - level of the master bias used (the dark's BIASLVL)
        predicted = ((model['intercept'] + model['slope'] * np.asarray(darks['TUBETEMP'][is_epoch_dark], float)
                      - np.asarray(darks['BIASLVL'][is_epoch_dark], float)) / np.asarray(darks['EXPTIME'][is_epoch_dark], float))
        measured = np.asarray(darks['dark_level'][is_epoch_dark], float)
        darks['predicted_bias_offset'][is_epoch_dark] = predicted
        is_valid = np.isfinite(predicted) & np.isfinite(measured)
        if is_valid.sum() < 3:
            continue
        scale, dark_current = np.polyfit(predicted[is_valid], measured[is_valid], 1)
        corrected = measured - scale * predicted
        darks['dark_level_corrected'][is_epoch_dark] = corrected
        results[epoch]['dark_comparison'] = dict(
            scale=float(scale), dark_current=float(dark_current),
            correlation=float(np.corrcoef(predicted[is_valid], measured[is_valid])[0, 1]),
            rms_before=float(np.std(measured[is_valid])), rms_after=float(np.std(corrected[is_valid])),
            number_of_darks=int(is_valid.sum()))
    return results


# ----------------------------------------------------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------------------------------------------------
def mark_epoch_boundaries(axis, boundaries, reference_time):
    for boundary in boundaries:
        axis.axvline(boundary.mjd - reference_time.mjd, color='k', ls='--', lw=1)


def plot_dark_trend(table, correlation_results, labels, boundaries, reference_time, camera, output_filename):
    """ Dark level vs time, vs each correlated quantity, the temperatures vs time, and the spatial pattern; epochs in
    different colours and markers."""
    dark_level = np.asarray(table['dark_level'], float)
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))

    def scatter_by_epoch(axis, x_values, y_values, with_labels=False):
        for epoch, label in enumerate(labels):
            is_epoch = table['epoch'] == epoch
            if is_epoch.any():
                axis.plot(np.asarray(x_values, float)[is_epoch], np.asarray(y_values, float)[is_epoch],
                          EPOCH_MARKERS[epoch % len(EPOCH_MARKERS)], ms=5, color=EPOCH_COLOURS[epoch % len(EPOCH_COLOURS)],
                          label=f'{label} ({is_epoch.sum()} darks)' if with_labels else None)

    # Panel 1: dark level vs time
    axis = axes[0, 0]
    scatter_by_epoch(axis, table['days_since_start'], dark_level, with_labels=True)
    mark_epoch_boundaries(axis, boundaries, reference_time)
    axis.set_xlabel(f"days since {reference_time.isot[:16]} UTC")
    axis.set_ylabel('dark level [e-/s]')
    axis.set_title(f'{camera}: dark level of {len(table)} darks')
    axis.legend(fontsize=8)
    axis.grid(alpha=0.3)

    # Panels 2-7: dark level vs each quantity (except time, shown in panel 1); correlations per epoch in the title
    panel_axes = [axes[0, 1], axes[0, 2], axes[1, 0], axes[1, 1], axes[1, 2], axes[2, 0]]
    for axis, (column, label) in zip(panel_axes, CORRELATED_QUANTITIES[1:]):
        scatter_by_epoch(axis, table[column], dark_level)
        title_parts = []
        for epoch in range(len(labels)):
            result = correlation_results[epoch][column] if epoch in correlation_results else None
            if result is not None:
                title_parts.append(f"r = {result['pearson']:+.2f}" if np.isfinite(result['pearson']) else result['note'])
        axis.set_title(' | '.join(title_parts) + ('  (per epoch)' if len(labels) > 1 else ''))
        axis.set_xlabel(label)
        axis.set_ylabel('dark level [e-/s]')
        axis.grid(alpha=0.3)

    # Panel 8: temperatures vs time
    axis = axes[2, 1]
    for column, colour in (('WMSTEMP', 'C0'), ('TUBETEMP', 'C1'), ('FOCTEMP', 'C2')):
        axis.plot(table['days_since_start'], np.asarray(table[column], float), 'o', ms=3, color=colour, label=column)
    mark_epoch_boundaries(axis, boundaries, reference_time)
    axis.set_xlabel(f"days since {reference_time.isot[:16]} UTC")
    axis.set_ylabel('temperature [C]')
    axis.legend(fontsize=8)
    axis.set_title('ambient temperatures at the dark exposures')
    axis.grid(alpha=0.3)

    # Panel 9: spatial pattern of the excess: mean region level of the brightest third minus the faintest third
    axis = axes[2, 2]
    region_levels = np.vstack(table['region_levels'])
    order = np.argsort(dark_level)
    third = max(len(order) // 3, 1)
    excess = region_levels[order[-third:]].mean(axis=0) - region_levels[order[:third]].mean(axis=0)
    image = axis.imshow(excess.reshape(GRID_SIZE, GRID_SIZE), origin='lower', cmap='viridis')
    for index, value in enumerate(excess):
        axis.text(index % GRID_SIZE, index // GRID_SIZE, f'{value:.4f}', ha='center', va='center', color='w', fontsize=9)
    plt.colorbar(image, ax=axis, label='e-/s')
    axis.set_xticks(range(GRID_SIZE))
    axis.set_yticks(range(GRID_SIZE))
    axis.set_xlabel('region column (x)')
    axis.set_ylabel('region row (y)')
    axis.set_title('excess of the brightest third of darks\nover the faintest third, per region (all epochs)')

    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()


def plot_bias_offset(darks, biases, bias_results, labels, boundaries, reference_time, camera, output_filename):
    """ Bias level vs tube temperature and vs time, measured dark level vs predicted bias offset, and the dark level
    before / after removing the bias offset; per epoch."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    for epoch, label in enumerate(labels):
        colour = EPOCH_COLOURS[epoch % len(EPOCH_COLOURS)]
        marker = EPOCH_MARKERS[epoch % len(EPOCH_MARKERS)]
        epoch_biases = biases[biases['epoch'] == epoch]
        epoch_darks = darks[darks['epoch'] == epoch]
        model = bias_results[epoch]['bias_model']
        comparison = bias_results[epoch]['dark_comparison']

        # Bias level vs tube temperature, with the epoch's model
        axis = axes[0, 0]
        axis.plot(epoch_biases['TUBETEMP'], epoch_biases['BIASLVL'], marker, ms=3, color=colour, alpha=0.6,
                  label=f'{label}: {len(epoch_biases)} biases')
        if model is not None:
            temperature_grid = np.linspace(np.nanmin(epoch_biases['TUBETEMP']), np.nanmax(epoch_biases['TUBETEMP']), 10)
            axis.plot(temperature_grid, model['intercept'] + model['slope'] * temperature_grid, '-', color=colour,
                      label=f"BIASLVL = {model['intercept']:.1f} + {model['slope']:.2f} T  (r = {model['correlation']:.2f}, rms {model['rms']:.1f})")

        # Bias level vs time: daily biases, and the master bias levels used by the darks
        axis = axes[0, 1]
        axis.plot(epoch_biases['days_since_start'], epoch_biases['BIASLVL'], marker, ms=3, color=colour, alpha=0.6,
                  label=f'{label}: daily biases')
        axis.plot(epoch_darks['days_since_start'], epoch_darks['BIASLVL'], 'x', ms=6, color=colour,
                  label=f'{label}: master bias level used by the darks')

        # Measured dark level vs predicted bias offset
        axis = axes[1, 0]
        axis.plot(epoch_darks['predicted_bias_offset'], epoch_darks['dark_level'], marker, ms=5, color=colour,
                  label=label if comparison is None else
                  f"{label}: slope {comparison['scale']:.2f}, dark current {comparison['dark_current']:.4f} e-/s, r = {comparison['correlation']:.3f}")
        if comparison is not None:
            offset_grid = np.linspace(np.nanmin(epoch_darks['predicted_bias_offset']), np.nanmax(epoch_darks['predicted_bias_offset']), 10)
            axis.plot(offset_grid, comparison['scale'] * offset_grid + comparison['dark_current'], '-', color=colour)

        # Dark level before (open) and after (filled) removing the bias offset
        axis = axes[1, 1]
        axis.plot(epoch_darks['days_since_start'], epoch_darks['dark_level'], marker, ms=5, mfc='none', color=colour,
                  label=f'{label}: measured' + (f" (rms {comparison['rms_before']:.4f})" if comparison else ''))
        axis.plot(epoch_darks['days_since_start'], epoch_darks['dark_level_corrected'], marker, ms=5, color=colour,
                  label=f'{label}: bias offset removed' + (f" (rms {comparison['rms_after']:.4f})" if comparison else ''))

    axes[0, 0].set_xlabel('telescope tube temperature TUBETEMP [C]')
    axes[0, 0].set_ylabel('bias level BIASLVL')
    axes[0, 0].set_title(f'{camera}: bias level vs tube temperature')
    axes[0, 1].set_xlabel(f"days since {reference_time.isot[:16]} UTC")
    axes[0, 1].set_ylabel('bias level BIASLVL')
    axes[0, 1].set_title('bias level vs time')
    axes[1, 0].set_xlabel('predicted bias offset / EXPTIME [e-/s]\n(bias level at the dark\'s tube temperature - master bias level)')
    axes[1, 0].set_ylabel('measured dark level [e-/s]')
    axes[1, 0].set_title('darks: measured level vs predicted bias offset')
    axes[1, 1].set_xlabel(f"days since {reference_time.isot[:16]} UTC")
    axes[1, 1].set_ylabel('dark level [e-/s]')
    axes[1, 1].set_title('dark level before / after removing the bias offset')
    for axis in (axes[0, 1], axes[1, 1]):
        mark_epoch_boundaries(axis, boundaries, reference_time)
    for axis in axes.flat:
        axis.legend(fontsize=7)
        axis.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()


# ----------------------------------------------------------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------------------------------------------------------
def parseCommandLine():
    parser = argparse.ArgumentParser(description='Dark and bias offset level trend of a camera, correlated with time, '
                                                 'CCD and ambient temperature',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--camera', required=True, help='instrument, e.g. ep60')
    parser.add_argument('--ndays', type=int, default=10, help='number of DAY-OBS dates, ending on --dayobs')
    parser.add_argument('--dayobs', default=None, help='last DAY-OBS, YYYYMMDD (default: today UTC)')
    parser.add_argument('--epoch-boundary', nargs='*', default=[],
                        help='UTC time(s), e.g. 2026-10-19T22:00, that split the data into epochs (e.g. a change of the '
                             'camera readout); frames are labelled and analysed per epoch')
    parser.add_argument('--readmode', default='full_frame')
    parser.add_argument('--outdir', default='plots/darktrend')
    parser.add_argument('--threads', type=int, default=4, help='parallel downloads')
    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'])
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    args.epoch_boundary = sorted(Time(boundary, format='isot', scale='utc') for boundary in args.epoch_boundary)
    return args


def print_summary(darks, biases, correlation_results, bias_results, labels, dates, camera, basename):
    print(f"\n{camera}: {len(darks)} darks and {len(biases)} biases, {dates[0]} - {dates[-1]}")
    for epoch, label in enumerate(labels):
        epoch_darks = darks[darks['epoch'] == epoch]
        epoch_biases = biases[biases['epoch'] == epoch]
        print(f"\n  epoch {label}: {len(epoch_darks)} darks, {len(epoch_biases)} biases")
        if len(epoch_darks) > 0:
            saturation_levels = sorted({float(level) for level in epoch_darks['SATURATE'] if np.isfinite(level)})
            print(f"    dark level {np.min(epoch_darks['dark_level']):.4f} .. {np.max(epoch_darks['dark_level']):.4f} e-/s "
                  f"(median {np.median(epoch_darks['dark_level']):.4f}); SATURATE in darks: {saturation_levels}")
        if epoch in correlation_results:
            print(f"    {'quantity':45s} {'Pearson r':>10s} {'Spearman rho':>13s} {'p-value':>10s}")
            for column, result in correlation_results[epoch].items():
                print(f"    {result['label']:45s} {result['pearson']:+10.2f} {result['spearman']:+13.2f} "
                      f"{result['spearman_p']:10.1e}  {result['note']}")
        model = bias_results[epoch]['bias_model']
        if model is None:
            print("    bias level model: not enough biases")
        else:
            print(f"    bias level: BIASLVL = {model['intercept']:.1f} + {model['slope']:.2f} x TUBETEMP  (r = {model['correlation']:.3f}, "
                  f"rms {model['rms']:.1f}, {model['number_of_biases']} biases, range {model['bias_level_range'][0]:.1f} .. "
                  f"{model['bias_level_range'][1]:.1f})")
        comparison = bias_results[epoch]['dark_comparison']
        if comparison is not None:
            print(f"    darks: measured = {comparison['scale']:.3f} x predicted bias offset + {comparison['dark_current']:.4f} e-/s "
                  f"(r = {comparison['correlation']:.3f}); rms {comparison['rms_before']:.4f} -> {comparison['rms_after']:.4f} e-/s "
                  f"after removing the bias offset")
    print(f"\n  plots and tables: {basename}*")


def main():
    args = parseCommandLine()
    os.makedirs(args.outdir, exist_ok=True)
    last_day = dt.datetime.strptime(args.dayobs, '%Y%m%d').date() if args.dayobs else dt.datetime.now(dt.timezone.utc).date()
    dates = [(last_day - dt.timedelta(days=offset)).strftime('%Y%m%d') for offset in range(args.ndays - 1, -1, -1)]
    labels = epoch_labels(args.epoch_boundary)

    # Darks: header data from the index, dark level measured from the images
    dark_records = query_frames(args.camera, dates, args.readmode, 'DARK', 'd91')
    if not dark_records:
        _log.error("no darks found")
        return
    darks = measure_all(dark_records, os.path.join(args.outdir, f'darklevel_{args.camera}_cache.ecsv'), args.threads)
    if darks is None or len(darks) == 0:
        _log.error("no dark could be measured")
        return
    # Biases: header data only
    biases = Table(rows=query_frames(args.camera, dates, args.readmode, 'BIAS', 'b91'))

    # Common time axis, epochs, derived quantities
    reference_time = Time(min(np.asarray(darks['DATE-OBS'], dtype=str)), format='isot', scale='utc')
    darks = assign_epochs(add_time_and_sun(darks, reference_time), args.epoch_boundary)
    biases = assign_epochs(add_time_and_sun(biases, reference_time), args.epoch_boundary)
    darks.sort('days_since_start')
    correlation_results = {epoch: correlations(darks[darks['epoch'] == epoch])
                           for epoch in range(len(labels)) if (darks['epoch'] == epoch).sum() >= 3}
    bias_results = apply_bias_offset_model(darks, biases, len(labels))

    basename = os.path.join(args.outdir, f'darklevel_{args.camera}_{dates[0]}-{dates[-1]}')
    darks.write(f'{basename}.ecsv', overwrite=True)
    biases.write(f'{basename}_biases.ecsv', overwrite=True)
    plot_dark_trend(darks, correlation_results, labels, args.epoch_boundary, reference_time, args.camera, f'{basename}.pdf')
    plot_bias_offset(darks, biases, bias_results, labels, args.epoch_boundary, reference_time, args.camera,
                     f'{basename}_bias_offset.pdf')
    print_summary(darks, biases, correlation_results, bias_results, labels, dates, args.camera, basename)


if __name__ == '__main__':
    main()
