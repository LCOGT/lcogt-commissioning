"""
Skeleton: query the LCO archive for flat field images of one camera for a given DAY-OBS.

Usage:
    export ARCHIVE_API_TOKEN=xxxxxxxx
    python flatfield_comaprison.py --camera fa16 --dayobs 20261001
"""
import argparse
import datetime as dt
import json
import logging
import astropy
import matplotlib.pyplot as plt
import os

import numpy as np

import lcocommissioning.common.lco_archive_utilities as lco_archive_utilities
from lcocommissioning.common.logging_config import setup_logging  # TODO: verify import path

_log = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# Authentication
# --------------------------------------------------------------------------
def get_archive_token():
    """Read the LCO archive API token from the ARCHIVE_API_TOKEN environment variable."""
    token = os.environ.get("ARCHIVE_API_TOKEN")
    if not token:
        raise SystemExit("Environment variable ARCHIVE_API_TOKEN is not set.")
    return token


def auth_headers(token):
    """HTTP headers used to authenticate requests against the LCO archive API."""
    return {"Authorization": f"Token {token}"}


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
    return (dt.datetime.utcnow() - dt.timedelta(days=1)).date()


def dayobs_to_window(dayobs):
    """Convert a DAY-OBS into a [start, end) UTC query window.
    TODO: adjust if the archive should be queried by the DAY-OBS keyword directly,
    or if the site's DAY-OBS boundary differs from 00:00 UTC."""
    start = dt.datetime.combine(dayobs, dt.time(0, 0))
    return start, start + dt.timedelta(days=1)


# --------------------------------------------------------------------------
# Archive access (adapter to lco_archive_utilities)
# --------------------------------------------------------------------------
def query_flats(dayobs,   camera,  obstype="SKYFLAT", rlevel=91, isMaster=False, site=None, filter = filter, readmode='full_frame'):
    """Return archive frame records for flats of `camera` taken on `dayobs`."""

    frames = lco_archive_utilities.get_frames_for_noisegainanalysis(dayobs=dayobs,
        camera=camera, obstype=['SKYFLAT'], rlevel=91,  site=site, isMaster=isMaster, filter=filter, readmode=readmode)
    _log.info(f"Date {dayobs} found {len(frames)} frames")
    return frames


def download_flat(token, frame):
    """Download one frame and return its image data (numpy array)."""
    # TODO: replace with the download function from lco_archive_utilities, e.g.:
    # return lco_archive_utilities.<download_function>(frame, headers=auth_headers(token))
    raise NotImplementedError


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Query LCO archive for flat fields of a camera on a DAY-OBS.")
    p.add_argument("--camera", required=True, help="Instrument identifier, e.g. fa16")
    p.add_argument("--dayobs", type=parse_dayobs, default=default_dayobs(),
                   help="DAY-OBS to query, YYYYMMDD (default: yesterday UTC)")
    p.add_argument("--ndays", type=int, default=1, help="Number of days to query starting from DAY-OBS (default: 1)")
    p.add_argument("--loglevel", default="DEBUG", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    p.add_argument("--site", default=None, help="Site identifier to query (default: None)")
    p.add_argument("--filter", default=None, nargs="*", help="Filter to query. If none is set, master flats in all filters will be downloaded. If set, only the specified filters will be queried, and daily individual flats will be compared against the master flat.")
    p.add_argument("--cachedir", default="flatcache", help="Directory to cache downloaded files (default: None)")
    args = p.parse_args()
    setup_logging(level=args.loglevel)
    _log.debug(f"Flat Field Comparison for DAY-OBS {args.dayobs} and camera {args.camera}")
    if not os.path.exists(args.cachedir):
        _log.info(f"Creating cache directory at {args.cachedir}")
        os.makedirs(args.cachedir)
    return args 


def plot_flats(entry, args, width=200, normalizeBy=None):

    fitsimage = astropy.io.fits.open(os.path.join(args.cachedir, entry['FILENAME']))
    normalizeImage = None
    normFileName = None
    imageFileName = entry['FILENAME']

    if normalizeBy is not None:
        normalizeImage = astropy.io.fits.open(os.path.join(args.cachedir, normalizeBy['FILENAME']))['SCI'].data
        fitsimage['SCI'].data /= normalizeImage
        normFileName = normalizeBy['FILENAME']

    plt.figure(figsize=(15//2,20//2))
    grid = plt.GridSpec(2, 1, wspace=0.1, hspace=0.1, height_ratios=[4, 1])
    plt.subplot(grid[0])

    plt.imshow(fitsimage['SCI'].data, origin='lower', cmap='gray', vmin=0.97, vmax=1.03)
    plt.colorbar()
    plt.title(f"{imageFileName}\n {'normalized by ' + normFileName if normFileName is not None else ''}\n{entry['INSTRUME']:<4s}   {entry['CONFMODE']:<15s}  Filter: {entry['FILTER']:<5s}   DateObs: {entry['DAY-OBS']}")
    plt.xlabel("Pixel X coordinate")
    plt.ylabel("Pixel Y coordinate")

    centery = fitsimage['SCI'].data.shape[0] // 2
    horizontalCut = fitsimage['SCI'].data[centery-width//2:centery+width//2, :]
 
    plt.axhline(y=centery-width//2, color='red', linestyle='--')
    plt.axhline(y=centery+width//2, color='red', linestyle='--')

    plt.subplot(grid[1])
    horizontalCut = np.median(horizontalCut, axis=0)
    horizontalCut /= np.median(horizontalCut)

    plt.axhline(y=1, color='red', linestyle='--')   
    plt.plot(horizontalCut)
    plt.xlabel("Pixel X coordinate")
    plt.ylabel("Median Pixel Value\nin cut")
    plt.ylim(0.97, 1.03)

    plt.tight_layout()
    figname=f"{entry['FILENAME'].replace(".fz","").replace(".fits","")}_{entry['CONFMODE']}_{entry['FILTER']}.pdf"
    plt.savefig(os.path.join(args.cachedir, figname), dpi=150)
    plt.close()
    fitsimage.close()


def main():
    args = parse_args()

    ### Get a libray of existing master flats
    superCalibrationDates = lco_archive_utilities.ArchiveDiskCrawler.get_last_n_days(args.ndays)

    masterflatlist = {} 
    


    ### Get the most up to datre superalibration flat and plot them
    for date in superCalibrationDates:
        frames = query_flats(date, args.camera,  site=args.site, filter=args.filter, isMaster=True, readmode=["central_2k_2x2", "full_frame"])
        for frame in frames:
            f = frame
            if f['CONFMODE'] not in masterflatlist:
                masterflatlist[f['CONFMODE']] = {}
            masterflatlist[f['CONFMODE']][f'{frame["FILTER"]}'] = f


    for confmode, flats in masterflatlist.items():
        for filtername, frame in flats.items():
            print(f" {confmode:<15s}  Filter: {filtername:>5s},  Filename: {frame['FILENAME']}, FrameID: {frame['frameid']   }")

            cachedImageName = os.path.join(args.cachedir, frame['FILENAME'])    
            if not os.path.exists(cachedImageName):
                fitsimage = lco_archive_utilities.download_from_archive(frame['frameid'])
                fitsimage.writeto(cachedImageName, overwrite=True)
                fitsimage.close()

    for filtername, frame in masterflatlist['full_frame'].items():
        plot_flats(frame, args)


    if args.filter is not None:
        readmode=[args.readmode if hasattr(args, 'readmode') else "full_frame",]
        dailyCalibrationDates = lco_archive_utilities.ArchiveDiskCrawler.get_last_n_days(2)    
        dailyFlatList = {}
        for date in dailyCalibrationDates:
            frames = query_flats(date, args.camera,  site=args.site, filter=args.filter, isMaster=False, readmode=readmode)
            for frame in frames:
                f = frame
                dailyFlatList[f'{frame["FILENAME"]}'] = f

    # Get the most recent daily falts and compare them against the super calibration frame to see what changed
        for filename,  frame in dailyFlatList.items():
                
            cachedImageName = os.path.join(args.cachedir, frame['FILENAME'])    
            if not os.path.exists(cachedImageName):
                fitsimage = lco_archive_utilities.download_from_archive(frame['frameid'])
                fitsimage.writeto(cachedImageName, overwrite=True)
                fitsimage.close()

        for filename, frame in dailyFlatList.items():
            plot_flats(frame, args, normalizeBy=masterflatlist['full_frame'][args.filter[0]])
        
if __name__ == "__main__":
    main()



