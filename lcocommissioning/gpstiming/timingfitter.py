import argparse
import functools

import logging

from astropy.io import ascii, fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy import signal
log = logging.getLogger(__name__)
plt.style.use("seaborn-v0_8")
from functools import partial
from multiprocessing import Pool
log = logging.getLogger(__name__)


def rampfunction (fractStartTime,  dt05, amplitude, bias):
    """

    :param fractStarttime: fraction of the exposure start time in seconds
    :param dt05: Time delay in seconds
    :return:
    """
    # triable is the absolute value of a sawtooth
    global period 
    deltaphase = 2* np.pi *dt05
    retval =  amplitude * np.abs(signal.sawtooth(2 * np.pi * fractStartTime / period + deltaphase)) + bias
    return retval

def maketestplots ():
    plt.figure()

    for dt in np.arange (0,1,0.25):
        fs = np.arange (0,1,0.01)
        y = rampfunction(fs,dt,1,1)
        plt.plot (fs,y,label=f"{dt}")
    plt.legend ()
    plt.savefig ('gpexamples.png')

#maketestplots()
#exit(0)

def do_gpsfitting (fractime, lightlevel, std, x,y, outpng, myperiod):
    global period
    period = myperiod
    bounds = [[-period/2, 0,-1000],[0,200000,200000]]
    try:
        (paramset, istat) = scipy.optimize.curve_fit(rampfunction, fractime, lightlevel, bounds=bounds, )
        delta = np.abs (lightlevel - rampfunction(fractime, paramset[0], amplitude = paramset[1], bias = paramset[2]))
        good = (delta < np.std (delta) * 4) & (lightlevel > 10)
        (paramset, istat) = scipy.optimize.curve_fit(rampfunction, fractime[good], lightlevel[good], bounds=bounds, )
        paramset[0] = paramset[0] / period # XXXXX Verify, but I think we need to correct the delay.
        perr = np.sqrt(np.diag(istat))
       
        if outpng is not None:
            plt.figure()
            _x = np.arange(0,period,0.01)
            bad = np.logical_not(good)
            plt.plot (fractime[bad], lightlevel[bad], 'x', c='grey', )
            plt.plot (fractime[good], lightlevel[good], '.',  c='black', label="data")
            plt.plot (_x,rampfunction(_x, paramset[0]*period, amplitude = paramset[1], bias = paramset[2]), '-', label=f"dt = {paramset[0]: 6.4f} +/- {perr[0]: 6.4f}s")
            plt.legend()
            plt.xlabel("Fractional UTSTART [s]")
            plt.ylabel ("Illumination Level [ADU]")
            plt.savefig (outpng,  bbox_inches='tight')
            plt.close()

        return (paramset,istat,x,y)
    except:
        log.exception ("fitting exception caught")
    return (None,None,x,y)


def getTestdata (testdeltaT=0, n=100, Amplitude=480*0.7, bias=500, ron=3.1, npixels=1):

    testdata_time = np.random.uniform(size=n)
    test_measurements = rampfunction(testdata_time,testdeltaT,Amplitude,0)
    perpixelnoise =  np.sqrt (test_measurements)* np.random.normal(scale=1, size=len (testdata_time))+ \
                     np.random.normal(scale=ron, size=len(testdata_time))
    test_measurements += bias  + perpixelnoise / np.sqrt(npixels)
    std = np.sqrt (test_measurements + ron**2) / np.sqrt(npixels)
    print (std)
    return (testdata_time, test_measurements, std)

def processfits(fitsname, makepng=False, title="", period=1, x=None, y=None):
    basename = fitsname[:-5]
    print (basename)
    f = fits.open (fitsname)
    binning = f[0].header['BLK']
    dimX = f[0].header['NAXIS1']
    dimY = f[0].header['NAXIS2']
    dimZ = f[0].header['NAXIS3']
    dt = np.zeros((dimY,dimX))
    dt_err = np.zeros((dimY,dimX))
    fracsec = f[2].data['fracsec']

    myargs = []
    log.info ("creating fitting tasks tasks")
    for xx in range (dimX):
        for yy in range (dimY):
            mean = f[0].data[:,yy,xx]
            std  = f[1].data[:,yy,xx]
            outpng=f"{basename}_gpsfit_{xx}_{yy}.pdf" if makepng else None
            if (x is not None) and (y is not None) and (xx == x)  and (yy == y):
                 outpng=f"{basename}_gpsfit_{xx}_{yy}.pdf" 
            myargs.append ( (fracsec,mean,std,xx,yy,outpng,period) )

            #paramset, pcov,x,y = do_gpsfitting(fracsec,mean, std, xx,yy, outpng)
            #perr = np.sqrt(np.diag(pcov))
            #dt[y,x] = paramset[0]
            #dt_err[y,x] = perr[0]
    log.info("Starting the fitting work")
    with Pool(processes=8) as pool:
        results = pool.starmap (do_gpsfitting, myargs)

        for result in results:
            paramset, pcov,x,y = result
            log.debug (f'Fitting result: {paramset}, {pcov}, {x}, {y}' )
            if paramset is not None:
                try:
                    perr = np.sqrt(np.diag(pcov))
                    dt[y,x] = paramset[0]
                    dt_err[y,x] = perr[0]
                except Exception as e:
                    log.error (f"Error while parsing fit {paramset} {e}")

    log.info ("Making nice graphs")

    mean_delay = np.median (dt)
    min = np.min (dt)
    max = np.max (dt)

    plt.figure()
    plt.imshow(dt, cmap='viridis', aspect='equal', vmin=min, vmax=max , origin='lower')
    plt.colorbar()
    plt.savefig (f'{basename}_gpsmap.png')

    plt.figure()
    meandt = np.mean (dt,axis=1)
    row = np.arange (dimY)
    row = row *binning + binning/2.

    fit =np.polyfit (row, meandt, 1)
    poly = np.poly1d (fit)
    for i in range(2):
        residual =  (meandt -poly (row))
        good = np.abs(residual) < 3 * np.std (residual)
        fit = np.polyfit (row[good], meandt[good], 1)
        poly = np.poly1d (fit)

    plt.plot (row, meandt,'.', label="Data")
    plt.plot (row, poly(row), label=str(fit))
    plt.legend()
    plt.xlabel ("Row number")
    plt.ylabel ("Delay of start of exposure [s]")
    plt.title (title)
    plt.savefig(f'{basename}_gps_perrow.png', bbox_inches='tight')


    plt.figure()
    meandt = np.mean (dt,axis=0)
    col = np.arange (dimX)
    col = col *binning + binning/2.

    fit =np.polyfit (col, meandt, 1)
    poly = np.poly1d (fit)
    for i in range(3):
        residual =  (meandt -poly (col))
        good = np.abs(residual) < 3 * np.std (residual)
        fit = np.polyfit (row[good], meandt[good], 1)
        poly = np.poly1d (fit)

    plt.plot (col, meandt,'.')
    plt.plot (col, poly(col), label=str(fit))
    plt.legend()
    plt.xlabel ("Col number")
    plt.ylabel ("Delay of start of exposure [s]")
    plt.title (title)
    plt.savefig(f'{basename}_gps_percolumn.png', bbox_inches='tight')





def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Calculate gps timing')
    parser.add_argument('inputfile', nargs=1)
    parser.add_argument('--title', type=str)
    parser.add_argument('--png', action='store_true')
    parser.add_argument('--x', type=int,default=None )
    parser.add_argument('--y', type=int,default=None )
    parser.add_argument('--period', default = 2, type=int, help="Signal period in seconds")

    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()
    args.inputfile = args.inputfile[0]
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args

def main():
    args=parseCommandLine()
    log.info (f'Reading in input file {args.inputfile}')
    processfits(args.inputfile, makepng=args.png, title=args.title, period=args.period, x=args.x, y=args.y)


if __name__ == '__main__':
    main()


