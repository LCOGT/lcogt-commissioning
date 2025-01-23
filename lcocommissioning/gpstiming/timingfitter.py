import argparse
import functools

import logging

from astropy.io import ascii, fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy import signal
log = logging.getLogger(__name__)
plt.style.use("seaborn")
from functools import partial
from multiprocessing import Pool


def rampfunction (fractStartTime, dt05, amplitude, bias):
    """

    :param fractStarttime: fraction of the exposure start time in seconds
    :param dt05: Time delay in seconds
    :return:
    """
    # triable is the absolute value of a sawtooth
    frequency = 1
    deltaphase = 2* np.pi *dt05
    retval =  amplitude * np.abs(signal.sawtooth(2 * np.pi * frequency * fractStartTime + deltaphase)) + bias
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

def do_gpsfitting (fractime, lightlevel, std, x,y, outpng):
    bounds = [[0, 0,-1000],[+1,10000,100000]]
    try:
        (paramset, istat) = scipy.optimize.curve_fit(rampfunction, fractime, lightlevel, bounds=bounds, )
        delta = np.abs (lightlevel - rampfunction(fractime, paramset[0], amplitude = paramset[1], bias = paramset[2]))
        good = (delta < np.std (delta) * 4) & (lightlevel > 10)
        (paramset, istat) = scipy.optimize.curve_fit(rampfunction, fractime[good], lightlevel[good], bounds=bounds, )
        perr = np.sqrt(np.diag(istat))
       
        if outpng is not None:
            plt.figure()
            _x = np.arange(0,1,0.01)
            bad = np.logical_not(good)
            plt.plot (fractime[bad], lightlevel[bad], 'x', c='grey', )
            plt.plot (fractime[good], lightlevel[good], '.',  c='black', label="data")
            plt.plot (_x,rampfunction(_x, paramset[0], amplitude = paramset[1], bias = paramset[2]), '-', label=f"dt = {paramset[0]: 6.4f} +/- {perr[0]: 6.4f}s")
            plt.legend()
            plt.xlabel("Fractional UTSTART [s]")
            plt.ylabel ("Illumination Level [ADU]")
            plt.savefig (outpng,  bbox_inches='tight')
            plt.close()

        return (paramset,istat,x,y)
    except:
        print ("fitting exception caught")
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

def readMeasurementData (fname):
    data = None
    with open (fname) as input:
        data = ascii.read(input.read(), delimiter=' ',guess=False, format='no_header')
    return data


def readsimplefile(filename):
    data = readMeasurementData(filename)

    fractime = data['col2'].data
    mean = data['col3'].data
    std = data['col4'].data

    #fractime,mean, std = getTestdata(testdeltaT=0.3,n=30)

    paramset, pcov = do_gpsfitting(fractime,mean, std=std)
    perr = np.sqrt(np.diag(pcov))

    print (f"paramters: {paramset}\nErrors: {perr}")
    plt.figure()

    x = np.arange(0,1,0.01)

    plt.errorbar (fractime, mean,yerr=std, fmt='.', label="data")
    plt.plot (x,rampfunction(x, paramset[0], amplitude = paramset[1], bias = paramset[2]), '-', label=f"dt = {paramset[0]: 6.4f} +/- {perr[0]: 6.4f}s")
    plt.legend()
    plt.savefig ("gpscorrelation.pdf")
    plt.close()



def processfits(fitsname, makepng=False, title=""):
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
    log.info ("creating tasks")
    for xx in range (dimX):
        for yy in range (dimY):
            mean = f[0].data[:,yy,xx]
            std  = f[1].data[:,yy,xx]
            outpng=f"gpsfit_{xx}_{yy}.png" if makepng else None

            myargs.append ( (fracsec,mean,std,xx,yy,outpng) )

            #paramset, pcov,x,y = do_gpsfitting(fracsec,mean, std, xx,yy, outpng)
            #perr = np.sqrt(np.diag(pcov))
            #dt[y,x] = paramset[0]
            #dt_err[y,x] = perr[0]
    log.info("Starting the fitting work")
    with Pool(processes=30) as pool:
        results = pool.starmap (do_gpsfitting, myargs)

        for result in results:
            paramset, pcov,x,y = result
            print (paramset, pcov, x, y)
            if paramset is not None:
                try:
                    perr = np.sqrt(np.diag(pcov))
                    dt[y,x] = paramset[0]
                    dt_err[y,x] = perr[0]
                except Exception as e:
                    print (f"Error while parsing fit {paramset} {e}")

    log.info ("Making nice graphs")

    plt.figure()
    plt.imshow(dt, cmap='viridis', aspect='equal', vmin=0.83, vmax=0.96, origin='lower')
    plt.colorbar()
    plt.savefig (f'{basename}_gpsmap.png')

    plt.figure()
    meandt = np.mean (dt,axis=1)

    row = np.arange (dimY)
    row = row *binning + binning /2.

    fit =np.polyfit (row, meandt, 1)
    fit = np.poly1d (fit)
    for i in range(3):
        residual =  (meandt -fit (row))
        good = np.abs(residual) < 3 * np.std (residual)
        fit = np.polyfit (row[good], meandt[good], 1)
        fit = np.poly1d (fit)

    plt.plot (row, meandt,'.')
    plt.plot (row, fit(row), label=fit)
    plt.legend()
    plt.xlabel ("Row number")
    plt.ylabel ("Delay of start of exposure [s]")
    plt.title (title)
    plt.savefig(f'{basename}_gps_perrow.png', bbox_inches='tight')


    plt.figure()
    plt.plot (row[:-2],residual[:-2], '.')
    plt.title (title)
    plt.xlabel ('Row number')
    plt.ylabel ('Timing residual from fir')
    plt.savefig(f'{basename}_gps_perrow_residual.png', bbox_inches='tight')



def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Calculate gps timing')
    parser.add_argument('inputfile', nargs=1)
    parser.add_argument('--title', type=str)
    parser.add_argument('--png', action='store_true')

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
    processfits(args.inputfile, makepng=args.png, title=args.title)


if __name__ == '__main__':
    main()


