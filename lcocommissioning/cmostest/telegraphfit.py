import glob
import math
import os.path

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.optimize import curve_fit
import datetime
import torch
from matplotlib import colormaps
plt.style.use("ggplot")

TIMINGRUNS=1

_sqrt2pi = math.sqrt(2 * math.pi)

def np_gaussian(x, A, x0, sigma):
    """
     x vector of values
     A scalar overall gaussian peak scaling
     x0: center of gaussian distribution
     sigma: std deviation of gaussian distribution

     :returns: value of the gaussian distribution at that point.
    """
    return A * np.exp(-1 * (x - x0) ** 2 / (2 * sigma ** 2)) / sigma /sigma /  _sqrt2pi


def np_tripple_gaussian_function(x, A, leftscale, rightscale, x0, sigma, delta):
    """Tri-modal distribution
    x: varible
    A: amplitude at center
    B: fraction of center amplitude in side lobes
    x0: center of central peak
    delta: offset left or right from center
    """
    return A * (np_gaussian(x, 1, x0, sigma) + np_gaussian(x, rightscale, x0 - delta, sigma) + np_gaussian(x, leftscale, x0 + delta,
                                                                                                  sigma))


def summed_ln_likelihood(x, A, leftscale, rightscale, x0, sigma, delta):
    sum = np.sum(np.log(np_tripple_gaussian_function(x, A, leftscale, rightscale, x0, sigma, delta)))
    return sum


def findditributionparamters(xvec, popt):
    """ derive the RTN coefficients for a series of readouts based on a maximum likelyhood fit.
        popt: 0-> overall scale
              1-> left scale
              2-> right scale
              3-> x0
              4 -> sigma
              5 -> delta
    """
    guess = [popt[0], popt[1], popt[2], popt[3], popt[4],popt[5]]

    def lnprob(parameters):
        return -1 * summed_ln_likelihood(xvec, parameters[0], parameters[1], parameters[2], parameters[3], parameters[4],parameters[5])

    res = scipy.optimize.minimize(lnprob, guess, bounds=(  (0,1),  (0,0.95), (0, 0.95), (np.min(xvec), np.max(xvec)), (1, 20), (0, 120)), )
    return res


def findx0(xvec, A, leftscale,rightscale, sigma, delta):
    """ Given a certain trimodal distribution, find the most likely value for series of values drawn from that distribution.
    For CMOS RTS: If we have the same exposure level on a given pixel, and N readouts, what is the most likely occurance?
    """

    def lnprob(paramters, xvec):
        # print (f"paramters: {paramters} data: {xvec}")
        x0 = paramters[0]
        return -1 * summed_ln_likelihood(xvec, A, leftscale,rightscale, x0, sigma, delta)

    start = datetime.datetime.utcnow()

    guess = [np.mean(xvec), ]
    res = scipy.optimize.minimize(lnprob, guess, xvec, )
    end = datetime.datetime.utcnow()
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
        # Do a three-peak fit
        popt_3, pcov_3 = curve_fit(np_tripple_gaussian_function, xdata=binscenters, ydata=histo1,
                                   p0=[1, 0.5, 0.5, np.median(data), 20, (np.max(data)-np.median(data))],
                                   bounds=[[0,0, 0, np.median(data)-15, 2., 0], [np.inf, 1, 1, np.median(data)+15, 200, np.max(data)-np.median(data)]],
                                   )
        err_3 = np.sum(np.sqrt(np.diag(pcov_3)))
        #print(f"Fit results 3 gauss: \n\t{popt_3}\n\t{err_3}")
        # Do a single peak fit
        popt_1, pcov_1 = curve_fit(np_gaussian, xdata=binscenters, ydata=histo1, p0=[1, np.mean(data), 5],
                                   bounds=[[0, np.min(data), 2], [np.inf, np.max(data), np.inf]],
                                   )
        err_1 = np.sum(np.sqrt(np.diag(pcov_1)))
        #print(f"Fit results 1 gauss: \n\t{popt_1}\n\t{err_1}")

        if 1.5*err_1 < err_3:
            print ("Single fit is better than tripple fit")
            popt_3[0] = popt_1[0]
            popt_3[1] = 0
            popt_3[2] = 0
            popt_3[3] = popt_1[1]
            popt_3[4] = popt_1[2]
            popt_3[5] = 0


        print (popt_3)
    except Exception as e:
        print("Fitting failed", e)
        popt_3 = None
    return popt_3

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
        plt.plot(xspace, np_tripple_gaussian_function(xspace, *popt), color=fitcolor, linewidth=2,
                 label=f"Fit: {popt[0]: 4.2f} {popt[1]: 4.2f} {popt[2]: 4.2f}  {popt[3]: 7.1f} "
                       f"{popt[4]:> 4.1f} {popt[5]:> 8.1f}")
    return popt


def analyseandplotdata(bestvec, probevec, name=None):

    plt.figure()

    popt_probe = fitandplot_binneddata(probevec, "Worst case", numbins=30, color='lightgreen', fitcolor='green', legend="worst pixel")
    p_best = fitandplot_binneddata(minvec, "Best case", numbins=10, color='yellow', fitcolor='darkblue', legend="best pixel")
    print (f"worst case fit {popt_probe}")

    if (popt_probe) is None:
        print ("Cannot fit", name)
        return
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15))
    plt.title (f"QHY600 best & {name}")
    plt.xlabel ("Signal [ADU]")
    plt.savefig(f"temp/distributionfit_{name}.png", bbox_inches="tight", dpi=150)
    plt.close()


    start = datetime.datetime.utcnow()
    for ii in range(TIMINGRUNS):
        res=findditributionparamters(probevec, popt_probe)
    end = datetime.datetime.utcnow()
    print(f"Likelihood fit:\n\t{res.success}\t{res.x}")
    print(f"Timing for likelyhood fit with knowing the answer : : {(end - start) / TIMINGRUNS}")


    x = range(len(probevec))
    means = [np.mean(probevec[0:x]) for x in range(len(probevec))]
    medians = [np.median(probevec[0:x]) for x in range(len(probevec))]
    mlhres = [findx0(probevec[0:x], popt_probe[0], popt_probe[1], popt_probe[3], popt_probe[4],popt_probe[5]) for x in range(len(probevec))]
    mlh = [res.x[0] for res in mlhres]
    plt.figure()
    plt.plot(x, probevec, '.')
    plt.plot(x, means, label="Mean of first N", color='blue')
    plt.plot(x, medians, label="Median of first N", color='red')
    plt.ylabel("Pixel value [ADU]")
    plt.xlabel ("Image number")
    plt.title ("Bias stacking")
    plt.plot(x, mlh,  color='aqua', linewidth=5, alpha=0.3)
    plt.plot(x, mlh,  color='aqua', linewidth=3, alpha=0.5)
    plt.plot(x, mlh, label="maximum likelihood fit", color='aqua', linewidth=1, alpha = 0.95)
    plt.legend(loc='upper center', )#bbox_to_anchor=(0.5, -0.15))
    plt.savefig(f"temp/meanvsmle{name}.png", bbox_inches="tight", dpi=150)



if __name__ == '__main__':
    minvec = readdistribution('mindistr.txt')

    listofdata = glob.glob('temp/*.txt')
    for f in listofdata:

        probevec = readdistribution(f)
        analyseandplotdata(minvec,probevec=probevec, name=os.path.basename(f) )
