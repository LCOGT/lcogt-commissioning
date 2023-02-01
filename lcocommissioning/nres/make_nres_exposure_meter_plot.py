import datetime
from datetime import timezone

import astropy.io.fits as fits
import logging
import argparse
from astropy.time import Time

import astropy.table
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

from lcocommissioning.common.common import dateformat

log = logging.getLogger(__name__)
EXPOSURE_METER = 'EXPOSURE_METER'

def get_nresobject(filename):
    log.debug (f"Open fits object {filename}")
    hdul = fits.open (filename)
    if EXPOSURE_METER in hdul and len (hdul[EXPOSURE_METER].data)>1:
        return hdul

    return None

def process_fitsobj (hdul):

    imagename = hdul['SPECTRUM'].header['ORIGNAME']
    objects =  hdul['SPECTRUM'].header['OBJECTS']
    object = hdul['SPECTRUM'].header['OBJECT']
    exptime = hdul['SPECTRUM'].header['EXPTIME']

    dateobs = Time(hdul['SPECTRUM'].header['DATE-OBS'], format='fits').to_datetime()+datetime.timedelta(seconds=exptime/2.)

    fibers = objects.split ('&')
    actfiber = 'FIB0COUNTS' if fibers[0] == 'none' else 'FIB2COUNTS'

    data = hdul[EXPOSURE_METER].data
    data = astropy.table.Table(data)
    mjd_start = data['MJD_START']
    utctime  = [Time(mjd_start[x], format='mjd').to_datetime(None) for x in range (len(mjd_start))]

    calibcounts = data['FIB1COUNTS']
    sciencecounts = data[actfiber]

    weighted_mjd = np.sum ( sciencecounts * (mjd_start-mjd_start[0])) / np.sum (sciencecounts) + mjd_start[0]
    weighted_dateobs = Time(weighted_mjd, format='mjd').to_datetime(None)


    fig = plt.figure()
    ax = plt.subplot(111)
    plt.plot (utctime, calibcounts, label="Calibration Fiber")
    plt.plot (utctime, sciencecounts, label="Science Fiber")
    plt.title (f'{imagename}\nExposure meter {object}')
    plt.xlabel(f"Time [UTC], exptime={exptime} s")
    plt.ylabel ("Flux count  in fiber aperture [ADU]")

    plt.axvline (dateobs, color='black', label=f"DATE-OBS\n{dateobs}")
    plt.axvline (weighted_dateobs, color='red', label=f"flux-weighted DATE-OBS\n{weighted_dateobs}")

    plt.gcf().autofmt_xdate()
    plt.setp(plt.gca().xaxis.get_minorticklabels(), rotation=30)
    plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=30)
    plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(60) if exptime < 600 else range (0,60,5)))
    plt.gca().xaxis.set_minor_locator(mdates.SecondLocator(bysecond=[0,15,30, 45,] if exptime<600 else [0,60]))
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y%m%d\n%H:%M'))
    plt.gca().grid(which='major')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.5),
          fancybox=True, shadow=True, ncol=2)
    plt.savefig(f'{imagename}.exposuremeter.png',  dpi=150, bbox_inches='tight', )


    plt.close()
    hdul.close()
    return fig

def parse_args():
    """Function to harvest and parse the commandline arguments for noise
    analysis"""

    parser = argparse.ArgumentParser(description='Plot NRES Exposure Meter')


    parser.add_argument('fitsfiles', type=str, nargs='+',
                        help='Fits files for cross talk measurement')

    parser.add_argument('--loglevel', dest='loglevel', default='INFO', choices=['DEBUG', 'INFO'],
                        help='Set the debug level')
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.loglevel.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args




def main ():
    args = parse_args()
    for image in args.fitsfiles:
        hdul = get_nresobject(image)
        if hdul is not None:
            process_fitsobj (hdul)


main()