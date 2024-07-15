import glob
import math
import os.path
import logging
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.optimize import curve_fit
import datetime
import re
import torch
from matplotlib import colormaps
plt.style.use("ggplot")
plt.rcParams['figure.figsize'] = [10, 5]

TIMINGRUNS=1
_logger = logging.getLogger(__name__)

_sqrt2pi = math.sqrt(2 * math.pi)

def np_gaussian(x, A, x0, sigma):
    """
     x vector of values
     A scalar overall gaussian peak scaling
     x0: center of gaussian distribution
     sigma: std deviation of gaussian distribution

     :returns: value of the gaussian distribution at that point.
    """
    return A * np.exp(-1 * (x - x0) ** 2 / (2 * sigma ** 2)) / sigma  / _sqrt2pi


def np_tripple_gaussian_function(x, A, x0,sigma, leftscale, rightscale,  delta):
    """Tri-modal distribution
    x: varible
    A: amplitude at center
    B: fraction of center amplitude in side lobes
    x0: center of central peak
    delta: offset left or right from center
    """
    return A * (np_gaussian(x, 1, x0, sigma) + np_gaussian(x, rightscale, x0 - delta, sigma) + np_gaussian(x, leftscale, x0 + delta,
                                                                                                  sigma))
def np_quintuple_gaussian_function(x, A, x0, sigma, leftscale1, rightscale1, delta1, leftscale2,rightscale2, delta2 ):
    return A * (np_gaussian(x, 1, x0, sigma) +
                np_gaussian(x, rightscale1, x0 - delta1, sigma) +
                np_gaussian(x, leftscale1,  x0 + delta1, sigma) +
                np_gaussian(x, rightscale2, x0 - delta2, sigma) +
                np_gaussian(x, leftscale2,  x0 + delta2, sigma)
                )



def summed_ln_likelihood_3(x, A, x0,sigma, leftscale, rightscale, delta):
    assert (sigma >0), "Negative sigma encountered"
    v = np_tripple_gaussian_function(x, A, x0,sigma, leftscale, rightscale, delta)
    return  np.sum(np.log(v+1e-12))


def summed_ln_likelihood_5(x, A, x0,sigma, leftscale1, rightscale1, delta1, leftscale2,rightscale2,delta2):
    assert (sigma >0), "Negative sigma encountered"
    v = np_quintuple_gaussian_function(x, A, x0, sigma, leftscale1, rightscale1, delta1, leftscale2,rightscale2,delta2)
    return  np.sum(np.log(v+1e-12))


def plotLHFunction (x, params, inputname):
    plt.figure ()
    plt.subplot(3,2,1)
    sigma = np.arange(0.01,30,0.01)
    lf = [-1. * summed_ln_likelihood_3(x, params[0], params[1], params[2], params[3], s, params[5]) for s in sigma]
    plt.plot (sigma, lf, label="\sigma")
    plt.title ("\sigma")

    plt.subplot(3,2,2)
    delta = np.arange(0,200,1)
    lf = [-1. * summed_ln_likelihood_3(x, params[0], params[1], params[2], params[3], params[4], s) for s in delta]
    plt.plot (delta, lf)
    plt.title ("\delta")


    plt.subplot(3,2,3)
    A = np.arange(0,params[0]*1.5,0.01)
    lf = [-1. * summed_ln_likelihood_3(x, s, params[1], params[2], params[3], params[4], params[5]) for s in A]
    plt.plot (A, lf,)
    plt.title ("Amplitude A")

    plt.subplot(3,2,5)
    Al = np.arange(0,1,0.01)
    lf = [-1. * summed_ln_likelihood_3(x, params[0], s, params[2], params[3], params[4], params[5]) for s in Al]
    plt.plot (Al, lf,)
    plt.title ("Amplitude Al")

    plt.subplot(3,2,6)
    Ar = np.arange(0,1,0.01)
    lf = [-1. * summed_ln_likelihood_3(x, params[0], params[1], s, params[3], params[4], params[5]) for s in Ar]
    plt.plot (Ar, lf,)
    plt.title ("Amplitude Ar")


    plt.savefig (f'temp/lhdetails_{inputname}.png')




def findditributionparamters(xvec, popt):
    """ derive the RTN coefficients for a series of readouts based on a maximum likelyhood fit.
        popt: 0-> overall scale
              1 -> x0
              2 -> sigma
              3-> left scale
              4-> right scale
              5 -> delta
    """
    guess = [popt[0], popt[1], popt[2], popt[3], popt[4],popt[5]]

    def lnprob(parameters, xvec):
        return -1. * summed_ln_likelihood_3(xvec, parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5])
    res = scipy.optimize.minimize(lnprob, guess, xvec,
                                  bounds=((popt[0], popt[0]),(np.min(xvec), np.max(xvec)), (1, 30),  (0, 0.95), (0, 0.95), (3,120) )
                                  )
    return res


def findx0(xvec, A, sigma, leftscale,rightscale, delta,leftscale2,rightscale2,delta2):
    """ Given a certain trimodal distribution, find the most likely value for series of values drawn from that distribution.
    For CMOS RTS: If we have the same exposure level on a given pixel, and N readouts, what is the most likely occurance?
    """

    def lnprob(paramters, mydata):
        # print (f"parameters: {parameters} data: {xvec}")
        x0 = paramters[0]
        return -1 * summed_ln_likelihood_5(mydata, A, x0, sigma, leftscale, rightscale, delta,leftscale2,rightscale2,delta2)

    guess = [np.mean(xvec), ]
    res = scipy.optimize.minimize(lnprob,guess, xvec)
    return (res)


def readdistribution(file):
    """ Rwead a serices of pixel values from a text file. One value per line. """
    return np.asarray(np.loadtxt(file))


def fitbinneddata (data, numbins=20):
    bins = np.linspace(np.min(data), np.max(data), numbins)
    histo1, bins1 = np.histogram(data, bins=bins, )
    histo1 = histo1 / np.max(histo1)
    binscenters = np.array([0.5 * (bins1[i] + bins1[i + 1]) for i in range(len(bins1) - 1)])
    popt = None
    try:
        ss_tot = np.sum((histo1-np.mean(histo1))**2)
        n = numbins
        # Do a three-peak fit
        popt_5, pcov_5 = curve_fit(np_quintuple_gaussian_function, xdata=binscenters, ydata=histo1,
                                   p0=[1, np.median(data), 5, 0.5, 0.5, (np.max(data)-np.median(data)), 0.5, 0.5, (np.max(data)-np.median(data)) ],
                                   bounds=[[0,np.median(data)-15, 2, 0, 0, 0, 0, 0, 0], [np.inf, np.median(data)+15, 30, 1, 1, np.max(data)-np.median(data),1, 1, np.max(data)-np.median(data)]],
                                   )
        err_5 = np.sum(np.sqrt(np.diag(pcov_5))) / (9+2) # 9 deg of freedom
        residuals = histo1- np_quintuple_gaussian_function(binscenters, *popt_5)
        ss_res = np.sum(residuals**2)
        r_squared_5 = 1 - ( (ss_res / (n-9-1)) / (ss_tot / (n-1)))
        # Do a three-peak fit
        popt_3, pcov_3 = curve_fit(np_tripple_gaussian_function, xdata=binscenters, ydata=histo1,
                                   p0=[1, np.median(data), 5, 0.5, 0.5, (np.max(data)-np.median(data))],
                                   bounds=[[0,np.median(data)-15, 2, 0, 0, 0], [np.inf, np.median(data)+15, 30, 1, 1, np.max(data)-np.median(data)]],
                                   )
        err_3 = np.sum(np.sqrt(np.diag(pcov_3))) / (7+2) # 7 deg of freedom
        residuals = histo1- np_tripple_gaussian_function(binscenters, *popt_3)
        ss_res = np.sum(residuals**2)
        r_squared_3 = 1 - ( (ss_res / (n-7-1)) / (ss_tot / (n-1)))
        # Do a single peak fit
        popt_1, pcov_1 = curve_fit(np_gaussian, xdata=binscenters, ydata=histo1, p0=[1, np.mean(data), 5],
                                   bounds=[[0, np.min(data), 2], [np.inf, np.max(data), np.inf]],
                                   )
        err_1 = np.sum(np.sqrt(np.diag(pcov_1))) / (5+2) # 5 deg of freedom
        residuals = histo1- np_gaussian(binscenters, *popt_1)
        ss_res = np.sum(residuals**2)
        r_squared_1 = 1 - ( (ss_res / (n-7-1)) / (ss_tot / (n-1)))
        errors = [err_1,err_3,err_5]

        rvalues = [r_squared_1,r_squared_3,r_squared_5]
        best = np.argmin(errors)
        bestr = np.argmax(rvalues)
        print ("Errors:", errors, " best ",best)
        print ("R^2:", rvalues, " best ", bestr)
        if bestr == 0: #single order
            print ("Single fit is best")
            popt_5[0] = popt_1[0]
            popt_5[1] = popt_1[1]
            popt_5[2] = popt_1[2]
            popt_5[3] = 0
            popt_5[4] = 0
            popt_5[5] = 0
            popt_5[6] = 0
            popt_5[7] = 0
            popt_5[8] = 0
        elif bestr == 1: # 3rd order
            print ("Triple fit is best")
            popt_5[0] = popt_3[0]
            popt_5[1] = popt_3[1]
            popt_5[2] = popt_3[2]
            popt_5[3] = popt_3[3]
            popt_5[4] = popt_3[4]
            popt_5[5] = popt_3[5]
            popt_5[6] = 0
            popt_5[7] = 0
            popt_5[8] = 0


    except Exception as e:
        _logger.exception("Fitting failed", e)
        popt_5 = None
    return popt_5

def fitandplot_binneddata(data, label, numbins=15, color = 'black',  fitcolor='black', legend=None):
    """ Derive RTS paramters for a given pixel based on a sqeuince of readouts. This is based on fitting gaussians to a binned data sample
        Returns a popt set of coefficients .
    """
    popt = fitbinneddata(data, numbins)
    counts, bins = np.histogram (data, bins=numbins)

    counts = counts / np.max (counts)
    plt.hist(bins[:-1], bins, weights=counts, color=color, label=legend)

    if popt is not None:
        xspace = np.linspace(np.min(data) - 20, np.max(data) + 20, 100)
        plt.plot(xspace, np_quintuple_gaussian_function(xspace, *popt), color=fitcolor, linewidth=2,
                 label=f"Fit: {popt[0]: 4.2f} {popt[1]: 4.2f} {popt[2]: 4.2f}  {popt[3]: 7.1f} "
                       f"{popt[4]:> 4.1f} {popt[5]:> 8.1f}")
    return popt



def plotlikelyhoods (function, popt, poptidx, min,max, fname):

    mypopt = np.popt
    xbase = np.arange(min,max,0.1)
    values = []
    for x in xbase:
        mypopt[poptidx] = x
        values.append(function ())


p = re.compile ("pixeldist_(\d*\.?\d*)_(\d*)_(\d*).txt")

def analyseandplotdata(minvec, probevec, name=None):

    pixel_noise, pixel_x, pixel_y = get_pixelxy(name)

    plt.figure()
    popt_probe = fitandplot_binneddata(probevec, "Worst case", numbins=30, color='lightgreen', fitcolor='green', legend=f"Pixel Noise {pixel_noise}")
    #p_best = fitandplot_binneddata(minvec, "Best case", numbins=10, color='yellow', fitcolor='darkblue', legend="best pixel")

    with np.printoptions(precision=3, suppress=True):
        print (f"Current case fit:\n\t{popt_probe}")

    if (popt_probe) is None:
        print ("Cannot fit", name)
        return

    plt.title (f"QHY600 pixel @ {pixel_x}/ {pixel_y} ")
    plt.xlabel ("Signal [ADU]")

    # Now do an unbinned fit
    res=findditributionparamters(probevec, popt_probe)

    xbase = np.arange(np.min(probevec)-20, np.max(probevec)+20,1)
    if (False):
        popt = res.x
        plt.plot(xbase, np_tripple_gaussian_function(xbase, *popt), color='aqua', linewidth=2,
             label=f"ML Fit: {popt[0]: 4.2f} {popt[1]: 4.2f} {popt[2]: 4.2f}  {popt[3]: 7.1f} "
                   f"{popt[4]:> 4.1f} {popt[5]:> 8.1f} {popt[6]: 7.1f} {popt[7]: 7.1f} {popt[8]: 7.1f}")

        with np.printoptions(precision=3, suppress=True):
            print(f"Likelihood fit:\n\tFit Success: {res.success}\n\t{res.x}")


    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
    plt.savefig(f"temp/distributionfit_{name}.png", bbox_inches="tight", dpi=150)
    plt.close()

    #plotLHFunction(probevec, popt_probe, name)

    x = range(len(probevec))
    means = [np.mean(probevec[0:x]) for x in range(len(probevec))]
    medians = [np.median(probevec[0:x]) for x in range(len(probevec))]
    mlhres = [findx0(probevec[0:x], popt_probe[0], popt_probe[2], popt_probe[3], popt_probe[4],popt_probe[5],popt_probe[6], popt_probe[7],popt_probe[8]) for x in range(len(probevec))]
    mlh = [res.x[0] for res in mlhres]
    plt.figure()
    plt.plot(x, probevec, '.')
    plt.plot(x, means, label="Mean of first N", color='blue')
    plt.plot(x, medians, label="Median of first N", color='red')

    plt.plot(x, mlh, color='aqua', linewidth=5, alpha=0.3)
    plt.plot(x, mlh, color='aqua', linewidth=3, alpha=0.5)
    plt.plot(x, mlh, label="maximum likelihood fit", color='aqua', linewidth=1, alpha = 0.95)

    plt.ylabel("Pixel value [ADU]")
    plt.xlabel ("Image number")
    plt.title (f"Bias stacking @ {pixel_x}/{pixel_y}")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)

    plt.savefig(f"temp/meanvsmle{name}.png", bbox_inches="tight", dpi=300)


def get_pixelxy(name):
    if name is not None:
        m = p.match(name)
        pixel_noise = m.group(1)
        pixel_x = m.group(2)
        pixel_y = m.group(3)
    else:
        pixel_noise = pixel_x = pixel_y = None
    return pixel_noise, pixel_x, pixel_y


if __name__ == '__main__':
    minvec = readdistribution('mindistr.txt')

    listofdata = glob.glob('temp/*.txt')
    for f in listofdata[0:20]:

        probevec = readdistribution(f)
        analyseandplotdata(minvec,probevec=probevec, name=os.path.basename(f) )
