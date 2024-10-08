import datetime
from datetime import timezone
import os.path
import astropy.io.fits as fits
import logging
import argparse
from astropy.time import Time

import astropy.table
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


log = logging.getLogger(__name__)
EXPOSURE_METER = 'EXPOSURE_METER'

def get_nresobject(filename):
    log.debug (f"Open fits object {filename}")
    hdul = fits.open (filename)
    if EXPOSURE_METER in hdul:
        return hdul

    return None

def process_fitsobj (hdul):

    imagename = hdul['SPECTRUM'].header['ORIGNAME']
    objects =  hdul['SPECTRUM'].header['OBJECTS']
    object = hdul['SPECTRUM'].header['OBJECT']
    exptime = hdul['SPECTRUM'].header['EXPTIME']

    # get mid-time of exposure according to shutter
    dateobs = Time(hdul['SPECTRUM'].header['DATE-OBS'], format='fits').to_datetime()+datetime.timedelta(seconds=exptime/2.)

    # get the active science fiber
    fibers = objects.split ('&')
    actfiber = 'FIB0COUNTS' if fibers[0] == 'none' else 'FIB2COUNTS'

    data = hdul[EXPOSURE_METER].data

    fig_exposure_meter = plt.figure(figsize=(8.5, 11))
    ax = plt.subplot(111)
    if len(data) > 1:

        data = astropy.table.Table(data)
        mjd_start = data['MJD_START']

        # get mid-times of exposure meter measurements
        utctime  = [Time(mjd_start[x], format='mjd').to_datetime(None) + datetime.timedelta(seconds=data['EXP_TIME'][0]/2)
                for x in range (len(mjd_start))]

        calibcounts = data['FIB1COUNTS']
        sciencecounts = data[actfiber]

        # get flux-weighted mid time of exposure. These are large numbers, avoid some numerical overflow.
        total_flux = np.sum (sciencecounts)
        weighted_mjd = np.sum ( sciencecounts * (mjd_start-mjd_start[0])) / total_flux + mjd_start[0]
        weighted_dateobs = Time(weighted_mjd, format='mjd').to_datetime(None)

        plt.plot (utctime, calibcounts, label="Calibration Fiber")
        plt.plot (utctime, sciencecounts, label="Science Fiber")
        plt.title (f'{imagename}\nExposure meter {object}\nTotal Exp Meter Flux: {total_flux:.2E}')
        plt.xlabel(f"Time [UTC], exptime={exptime} s")
        plt.ylabel ("Flux count  in fiber aperture [ADU]")

        # mark open shutter and flux based mid-time of exposure
        plt.axvline (dateobs, color='black', label=f"mid-time\n{dateobs}")
        plt.axvline (weighted_dateobs, color='red', label=f"flux-weighted mid-time\n{weighted_dateobs}")

        # some more fancy axis formating.
        plt.gcf().autofmt_xdate()
        plt.setp(plt.gca().xaxis.get_minorticklabels(), rotation=30)
        plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=30)
        plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(60) if exptime < 600 else range (0,60,5)))
        plt.gca().xaxis.set_minor_locator(mdates.SecondLocator(bysecond=[0,15,30, 45,] if exptime<600 else [0,60]))
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y%m%d\n%H:%M'))
        plt.gca().grid(which='major')

        # I am legend.
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
          fancybox=True, shadow=True, ncol=2)

        log.debug (f'{imagename} Texp: {exptime} EM Texp: {data["EXP_TIME"][0]} Total Flux: {total_flux: .2E} '
                  f'time interval {(utctime[-1] - utctime[0]).seconds / (len (utctime)-1): 5.2f}')
    else:
        plt.title (" No Exposure meter data available")
        # TODO: Maybe make it a fancier note?
        log.debug (f"Not processing image {imagename} since not enough data from the exposure meter")

    hdul.close()
    return fig_exposure_meter

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
            fig = process_fitsobj (hdul)
            #fig.tight_layout()
            log.info ("Saving to pdf")
            pdffile = f'{os.path.basename(image).replace (".fits.fz", "exposure_meter.pdf")}'
            print (pdffile)
            pdf =  PdfPages(pdffile)
            pdf.savefig(fig)
            plt.close(fig)
            pdf.close()


main()