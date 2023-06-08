from astropy.io import ascii, fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy import signal
import sys




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

def do_gpsfitting (fractime, lightlevel, std=None, outpng = None):
    bounds = [[0, 0,-1000],[+1,10000,100000]]
    (paramset, istat) = scipy.optimize.curve_fit(rampfunction, fractime, lightlevel, bounds=bounds, sigma=std)
    perr = np.sqrt(np.diag(istat))
    if outpng is not None:
        plt.figure()

        x = np.arange(0,1,0.01)

        plt.errorbar (fractime, lightlevel,yerr=std, fmt='.', label="data")
        plt.plot (x,rampfunction(x, paramset[0], amplitude = paramset[1], bias = paramset[2]), '-', label=f"dt = {paramset[0]: 6.4f} +/- {perr[0]: 6.4f}s")
        plt.legend()
        plt.savefig (outpng)
        plt.close()

    return (paramset,istat)


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



def processfits(fitsname):
    f = fits.open (fitsname)
    binning = f[0].header['BLK']

    dimX = f[0].header['NAXIS1']
    dimY = f[0].header['NAXIS2']
    dimZ = f[0].header['NAXIS3']
    BLK = f[0].header['BLK']
    dt = np.zeros((dimY,dimX))
    dt_err = np.zeros((dimY,dimX))

    fracsec = f[2].data['fracsec']
    for xx in range (dimX):
        for yy in range (dimY):
            mean = f[0].data[:,yy,xx]
            std  = f[1].data[:,yy,xx]
            outpng=f"gpsfit_{xx}_{yy}.png"
            paramset, pcov = do_gpsfitting(fracsec,mean, std=std)
            perr = np.sqrt(np.diag(pcov))

            dt[yy,xx] = paramset[0]
            dt_err[yy,xx] = perr[0]

    plt.figure()
    plt.imshow(dt)
    plt.colorbar()
    plt.savefig ('gpsmap.png')

    plt.figure()
    meandt = np.mean (dt,axis=1)

    row = np.arange (dimY)
    row = row *binning

    fit =np.polyfit (row, meandt, 1)
    fit = np.poly1d (fit)
    plt.plot (row, meandt,'.')
    plt.plot (row, fit(row), label=fit)
    plt.legend()
    plt.xlabel ("Row number")
    plt.ylabel ("Delayed of start of exposure")



    plt.savefig('gpstime-perrow')

    plt.figure()
    row = np.arange(dimY) * BLK
    time_offset = np.average(dt, axis=1)
    plt.plot (row, time_offset,'.')
    plt.xlabel('row number')
    plt.ylabel('dt [s] absolute - timestamp')
    c = np.polyfit(row,time_offset,1)
    plt.plot (row,(np.poly1d(c))(row),'-', label=f"{c[0]:5.4e}*row+{c[1]:5.4}")
    print (c)
    plt.legend()
    plt.savefig ("gps_row_time.png")



processfits(sys.argv[1])

#readsimplefile(sys.argv[1])


