"""
Compare the daily CCD calibration frames of one camera with the current master (super-calibration) frames.

For a given DAY-OBS:
  1. Find the most recent master flat (per readout mode and filter), master bias and master dark (per readout mode) in
     the archive within the --ndays days ending on DAY-OBS, download them and plot each full-frame master.
  2. Download the individual processed calibration frames taken on DAY-OBS and compare each with its master:
       flat   ratio       daily flat / master flat of the same filter   (dimensionless, around 1)
       bias   difference  daily bias - master bias                      [e-]
       dark   difference  daily dark - master dark                      [e-/s]
     Daily flats are only compared for the filters given with --filter.
  3. Per calibration type, a summary of the daily comparison over the day: median and robust rms of every residual.

BANZAI conventions used: processed daily biases are overscan-corrected but not master-bias subtracted; processed darks
and master darks are bias-subtracted and normalised to one second of exposure time (e-/s).

Usage:
    export ARCHIVE_API_TOKEN=xxxxxxxx
    flatcomp --camera ep60 --dayobs 20261017 --filter V --caltypes flat bias dark
"""
import argparse
import datetime as dt
import logging
import os
import re

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import lcocommissioning.common.lco_archive_utilities as lco_archive_utilities
from lcocommissioning.common.logging_config import setup_logging

_log = logging.getLogger(__name__)

# Archive OBSTYPE of each calibration type
OBSTYPE = {'flat': 'SKYFLAT', 'bias': 'BIAS', 'dark': 'DARK'}

# How a daily frame is compared with its master: ratio for the multiplicative flat, difference for additive bias / dark
COMPARISON = {'flat': 'ratio', 'bias': 'difference', 'dark': 'difference'}

# Physical unit of the image values
UNIT = {'flat': 'relative', 'bias': 'e-', 'dark': 'e-/s'}

# Readout modes of the masters; daily frames are compared in full-frame mode only
MASTER_READOUT_MODES = ['central_2k_2x2', 'full_frame']
DAILY_READOUT_MODES = ['full_frame']

# Fixed display range of flats and flat ratios
FLAT_DISPLAY_RANGE = (0.97, 1.03)

# Half width of the display range of bias and dark images, in units of the robust standard deviation
DISPLAY_RANGE_SIGMA = 5

# Individual frames carry a 4-digit sequence number and a frame-type code, e.g. ...-20261017-0002-b91.fits.fz;
# masters do not (e.g. ...-20261017-bias-bin1x1.fits.fz)
INDIVIDUAL_FRAME_PATTERN = re.compile(r'-(\d{4})-[a-z]\d{2}\.fits')


# --------------------------------------------------------------------------
# DAY-OBS handling
# --------------------------------------------------------------------------
def parse_dayobs(dayobs):
    """Parse a DAY-OBS string (YYYYMMDD or YYYY-MM-DD) into a date."""
    for fmt in ("%Y%m%d", "%Y-%m-%d"):
        try:
            return dt.datetime.strptime(dayobs, fmt).date()
        except ValueError:
            pass
    raise argparse.ArgumentTypeError(f"Invalid DAY-OBS '{dayobs}', expected YYYYMMDD")


def default_dayobs():
    """Default DAY-OBS: yesterday (UTC), i.e. the most recently completed night."""
    return (dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=1)).date()


def dates_ending_on(dayobs, number_of_days):
    """The number_of_days DAY-OBS strings (YYYYMMDD) ending on dayobs, oldest first."""
    return [(dayobs - dt.timedelta(days=offset)).strftime('%Y%m%d') for offset in range(number_of_days - 1, -1, -1)]


# --------------------------------------------------------------------------
# Archive access
# --------------------------------------------------------------------------
def is_individual_frame(filename):
    """True for an individual processed frame, False for a master (super-calibration) frame."""
    return INDIVIDUAL_FRAME_PATTERN.search(filename) is not None


def query_calibrations(dayobs, camera, calibration_type, is_master, site=None, filters=None, readmodes=DAILY_READOUT_MODES):
    """Archive records of processed (RLEVEL 91) calibration frames of one type taken on dayobs.

    The archive query with is_master=False also returns masters; they are removed here, so that the daily list holds
    individual frames only.
    """
    frames = lco_archive_utilities.get_frames_for_noisegainanalysis(
        dayobs=dayobs, camera=camera, obstype=[OBSTYPE[calibration_type]], rlevel=91, site=site, isMaster=is_master,
        filter=filters if calibration_type == 'flat' else None, readmode=readmodes)
    if not is_master:
        frames = [frame for frame in frames if is_individual_frame(frame['FILENAME'])]
    _log.info(f"{dayobs}: {len(frames)} {'master' if is_master else 'daily'} {calibration_type} frames")
    return list(frames)


def download_to_cache(frame, cachedir):
    """Download a frame into the cache directory unless it is there already; return the local path."""
    cached_path = os.path.join(cachedir, frame['FILENAME'])
    if not os.path.exists(cached_path):
        hdulist = lco_archive_utilities.download_from_archive(frame['frameid'])
        hdulist.writeto(cached_path, overwrite=True)
        hdulist.close()
    return cached_path


def master_key(frame, calibration_type):
    """Masters are kept per readout mode, flats additionally per filter."""
    if calibration_type == 'flat':
        return frame['CONFMODE'], frame['FILTER']
    return frame['CONFMODE'], None


def find_masters(args, calibration_type):
    """Most recent master per readout mode (and filter, for flats) within the --ndays days ending on DAY-OBS."""
    masters = {}
    for date in dates_ending_on(args.dayobs, args.ndays):  # oldest first, so later masters replace earlier ones
        for frame in query_calibrations(date, args.camera, calibration_type, is_master=True, site=args.site,
                                        filters=args.filter, readmodes=MASTER_READOUT_MODES):
            masters[master_key(frame, calibration_type)] = frame
    return masters


# --------------------------------------------------------------------------
# Image comparison and plots
# --------------------------------------------------------------------------
def read_science_image(path):
    """Science image (float32) and header of a BANZAI frame."""
    with fits.open(path) as hdulist:
        return hdulist['SCI'].data.astype(np.float32), hdulist['SCI'].header


def robust_statistics(image):
    """Median and robust standard deviation (MAD scaled to a Gaussian sigma) of an image, from every 4th pixel."""
    sampled = image[::4, ::4]
    sampled = sampled[np.isfinite(sampled)]
    median = float(np.median(sampled))
    return median, float(1.4826 * np.median(np.abs(sampled - median)))


def plot_calibration(frame, calibration_type, args, master=None, cut_width=200):
    """Plot a calibration frame, or its comparison with a master, with a horizontal cut through the centre.

    flat: the image (or the ratio to the master) on a fixed +-3% scale; the cut is normalised to its median.
    bias / dark: the image (or the difference to the master) on a scale of median +- 5 robust sigma.

    :return: dict with the robust statistics of the plotted image
    """
    image, header = read_science_image(os.path.join(args.cachedir, frame['FILENAME']))
    comparison_text = ''
    if master is not None:
        master_image, _ = read_science_image(os.path.join(args.cachedir, master['FILENAME']))
        if COMPARISON[calibration_type] == 'ratio':
            image = image / master_image
            comparison_text = f"divided by {master['FILENAME']}"
        else:
            image = image - master_image
            comparison_text = f"minus {master['FILENAME']}"

    median, robust_std = robust_statistics(image)
    if calibration_type == 'flat':
        display_range = FLAT_DISPLAY_RANGE
    else:
        half_width = DISPLAY_RANGE_SIGMA * max(robust_std, 1e-6)
        display_range = (median - half_width, median + half_width)
    outlier_fraction = float(np.mean(np.abs(image[::4, ::4] - median) > DISPLAY_RANGE_SIGMA * max(robust_std, 1e-6)))

    plt.figure(figsize=(7.5, 10))
    grid = plt.GridSpec(2, 1, wspace=0.1, hspace=0.1, height_ratios=[4, 1])

    # Image
    plt.subplot(grid[0])
    plt.imshow(image, origin='lower', cmap='gray', vmin=display_range[0], vmax=display_range[1])
    colourbar = plt.colorbar()
    colourbar.set_label(UNIT[calibration_type] if COMPARISON[calibration_type] == 'difference' or master is None else 'ratio')
    description = f"Filter: {frame['FILTER']:<5s}" if calibration_type == 'flat' else f"{calibration_type}, EXPTIME {header.get('EXPTIME', 0):.0f} s"
    plt.title(f"{frame['FILENAME']}\n{comparison_text}\n{frame['INSTRUME']:<4s}   {frame['CONFMODE']:<15s}  {description}   "
              f"DAY-OBS: {frame['DAY-OBS']}\nmedian {median:.4g}  robust rms {robust_std:.3g} {UNIT[calibration_type]}  "
              f"beyond {DISPLAY_RANGE_SIGMA} sigma: {outlier_fraction * 100:.2f}%", fontsize=9)
    plt.xlabel("Pixel X coordinate")
    plt.ylabel("Pixel Y coordinate")
    centre_row = image.shape[0] // 2
    plt.axhline(y=centre_row - cut_width // 2, color='red', linestyle='--')
    plt.axhline(y=centre_row + cut_width // 2, color='red', linestyle='--')

    # Horizontal cut: median over the band of rows around the centre
    plt.subplot(grid[1])
    horizontal_cut = np.median(image[centre_row - cut_width // 2:centre_row + cut_width // 2, :], axis=0)
    if calibration_type == 'flat':
        horizontal_cut /= np.median(horizontal_cut)
        plt.axhline(y=1, color='red', linestyle='--')
        plt.ylim(*FLAT_DISPLAY_RANGE)
        plt.ylabel("Median pixel value\nin cut (normalised)")
    else:
        cut_median, cut_std = np.median(horizontal_cut), 1.4826 * np.median(np.abs(horizontal_cut - np.median(horizontal_cut)))
        plt.axhline(y=cut_median, color='red', linestyle='--')
        plt.ylim(cut_median - 10 * max(cut_std, 1e-6), cut_median + 10 * max(cut_std, 1e-6))
        plt.ylabel(f"Median pixel value\nin cut [{UNIT[calibration_type]}]")
    plt.plot(horizontal_cut)
    plt.xlabel("Pixel X coordinate")

    label = frame['FILTER'] if calibration_type == 'flat' else calibration_type
    comparison_suffix = f"_{COMPARISON[calibration_type]}_to_master" if master is not None else ''
    figure_name = f"{frame['FILENAME'].replace('.fz', '').replace('.fits', '')}_{frame['CONFMODE']}_{label}{comparison_suffix}.pdf"
    plt.savefig(os.path.join(args.cachedir, figure_name), dpi=150, bbox_inches='tight')
    plt.close()
    return dict(filename=frame['FILENAME'], median=median, robust_std=robust_std, outlier_fraction=outlier_fraction,
                exptime=float(header.get('EXPTIME', 0)), filter=frame['FILTER'])


def plot_daily_summary(comparison_statistics, calibration_type, args):
    """Median and robust rms of every daily comparison over the day, in frame sequence order."""
    sequence_numbers = [int(INDIVIDUAL_FRAME_PATTERN.search(entry['filename']).group(1)) for entry in comparison_statistics]
    order = np.argsort(sequence_numbers)
    sequence_numbers = np.array(sequence_numbers)[order]
    medians = np.array([entry['median'] for entry in comparison_statistics])[order]
    robust_stds = np.array([entry['robust_std'] for entry in comparison_statistics])[order]
    labels = np.array([entry['filter'] if calibration_type == 'flat' else f"{entry['exptime']:.0f}s"
                       for entry in comparison_statistics])[order]

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    comparison_name = 'ratio to master' if COMPARISON[calibration_type] == 'ratio' else 'difference to master'
    axes[0].plot(sequence_numbers, medians, 'o-')
    axes[0].axhline(1 if COMPARISON[calibration_type] == 'ratio' else 0, color='k', lw=0.5)
    axes[0].set_ylabel(f'median of {comparison_name}\n[{UNIT[calibration_type]}]')
    axes[1].plot(sequence_numbers, robust_stds, 'o-')
    axes[1].set_ylabel(f'robust rms of {comparison_name}\n[{UNIT[calibration_type]}]')
    axes[1].set_xlabel('frame sequence number')
    for sequence_number, median, label in zip(sequence_numbers, medians, labels):
        axes[0].annotate(label, (sequence_number, median), textcoords='offset points', xytext=(0, 5), fontsize=6, ha='center')
    axes[0].set_title(f"{args.camera} {args.dayobs:%Y%m%d}: daily {calibration_type} frames compared with the master "
                      f"({len(comparison_statistics)} frames)")
    for axis in axes:
        axis.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(args.cachedir, f"{args.camera}_{args.dayobs:%Y%m%d}_{calibration_type}_daily_summary.pdf"))
    plt.close()


# --------------------------------------------------------------------------
# Per calibration type
# --------------------------------------------------------------------------
def process_calibration_type(calibration_type, args):
    """Masters of one calibration type: download and plot; daily frames of DAY-OBS: compare with their master."""
    masters = find_masters(args, calibration_type)
    for (readmode, filter_name), master in sorted(masters.items(), key=lambda item: str(item[0])):
        print(f" {calibration_type:<5s} {readmode:<15s} {'Filter: ' + filter_name if filter_name else '':<14s} "
              f"Filename: {master['FILENAME']}, FrameID: {master['frameid']}")
        download_to_cache(master, args.cachedir)
    for (readmode, _), master in masters.items():
        if readmode == 'full_frame':
            plot_calibration(master, calibration_type, args)

    # Daily flats only for the requested filters; daily bias and dark always
    if calibration_type == 'flat' and args.filter is None:
        return
    daily_frames = query_calibrations(args.dayobs.strftime('%Y%m%d'), args.camera, calibration_type, is_master=False,
                                      site=args.site, filters=args.filter)
    daily_frames = sorted(daily_frames, key=lambda frame: frame['FILENAME'])
    if args.max_daily_frames is not None:
        daily_frames = daily_frames[:args.max_daily_frames]

    comparison_statistics = []
    for frame in daily_frames:
        master = masters.get(master_key(frame, calibration_type))
        if master is None:
            _log.warning(f"no master {calibration_type} for {frame['FILENAME']} ({frame['CONFMODE']}, {frame['FILTER']}); skipped")
            continue
        download_to_cache(frame, args.cachedir)
        comparison_statistics.append(plot_calibration(frame, calibration_type, args, master=master))
    if comparison_statistics:
        plot_daily_summary(comparison_statistics, calibration_type, args)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Compare the daily calibration frames (flat, bias, dark) of a camera on a DAY-OBS with the master frames.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--camera", required=True, help="Instrument identifier, e.g. fa16")
    p.add_argument("--dayobs", type=parse_dayobs, default=default_dayobs(),
                   help="DAY-OBS of the daily calibration frames, YYYYMMDD (default: yesterday UTC)")
    p.add_argument("--ndays", type=int, default=1,
                   help="Number of days, ending on DAY-OBS, searched for the most recent master frames")
    p.add_argument("--caltypes", nargs="+", default=['flat', 'bias', 'dark'], choices=['flat', 'bias', 'dark'],
                   help="Calibration types to process")
    p.add_argument("--loglevel", default="DEBUG", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    p.add_argument("--site", default=None, help="Site identifier to query (default: None)")
    p.add_argument("--filter", default=None, nargs="*",
                   help="Flat filters to query. If none is set, master flats in all filters are downloaded and no daily "
                        "flats are compared. If set, only these filters are queried, and the daily flats of DAY-OBS are "
                        "compared with the master flat of the same filter.")
    p.add_argument("--max-daily-frames", type=int, default=None,
                   help="Maximum number of daily frames per calibration type to download and compare (default: all; "
                        "a day can have ~50 biases of ~20 MB each)")
    p.add_argument("--cachedir", default="flatcache", help="Directory to cache downloaded files and to write plots")
    args = p.parse_args()
    setup_logging(level=args.loglevel)
    _log.debug(f"Calibration comparison for DAY-OBS {args.dayobs} and camera {args.camera}")
    if not os.path.exists(args.cachedir):
        _log.info(f"Creating cache directory at {args.cachedir}")
        os.makedirs(args.cachedir)
    return args


def main():
    args = parse_args()
    for calibration_type in args.caltypes:
        process_calibration_type(calibration_type, args)


if __name__ == "__main__":
    main()
