import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.optimize import curve_fit
import datetime
import torch

"""

"""


_sqrt2pi = math.sqrt(2 * math.pi)


def np_gaussian(x, A, x0, sigma):
    """
     x vector of values
     A scalar overall gaussian peak scaling
     x0: center of gaussian distribution
     sigma: std deviation of gaussian distribution

     :returns: value of the gaussian distribution at that point.
    """
    return A * np.exp(-1 * (x - x0) ** 2 / (2 * sigma ** 2)) / sigma / sigma / _sqrt2pi


def np_tripple_gaussian_function(x, A, B, x0, sigma, delta):
    """Tri-modal distribution
    x: varible
    A: amplitude at center
    B: fraction of center amplitude in side lobes
    x0: center of central peak
    delta: offset left or right from center
    """
    return A * (np_gaussian(x, 1, x0, sigma) + np_gaussian(x, B, x0 - delta, sigma) + np_gaussian(x, B, x0 + delta,
                                                                                                  sigma))


def summed_ln_likelihood(x, A, B, x0, sigma, delta):
    sum = np.sum(np.log(np_tripple_gaussian_function(x, A, B, x0, sigma, delta)))
    return sum


def findditributionparamters(xvec, popt):
    """ derive the RTN coefficients for a series of readouts based on a maximum likelyhood fit.
    """
    guess = [popt[1], popt[2], popt[3], popt[4]]

    def lnprob(parameters):
        return -1 * summed_ln_likelihood(xvec, 1, parameters[0], parameters[1], parameters[2], parameters[3])

    res = scipy.optimize.minimize(lnprob, guess, bounds=((0, 1.0), (np.min(xvec), np.max(xvec)), (1, 100), (0, 100)), )
    print(f"Likelihood fit:\n\t{res}")
    return res


def findx0(xvec, A, B, sigma, delta):
    """ Given a certain trimodal distribution, find the most likely value for series of values drawn from that distribution.
    For CMOS RTS: If we have the same exposure level on a given pixel, and N readouts, what is the most likely occurance?
    """

    def lnprob(paramters, xvec):
        # print (f"paramters: {paramters} data: {xvec}")
        x0 = paramters[0]
        return -1 * summed_ln_likelihood(xvec, 1, B, x0, sigma, delta)

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
    pcov = None
    try:
        # Do a three-peak fit
        popt_3, pcov_3 = curve_fit(np_tripple_gaussian_function, xdata=binscenters, ydata=histo1,
                                   p0=[0.5, 0.25, np.mean(data), 5, 100],
                                   bounds=[[0, 0, np.min(data), 2.5, 0], [np.inf, 2, np.max(data), np.inf, np.inf]],
                                   )
        err_3 = np.sum(np.sqrt(np.diag(pcov_3)))
        # print(f"Fit results 3 gauss: \n\t{popt_3}\n\t{err_3}")
        # Do a single peak fit
        popt_1, pcov_1 = curve_fit(np_gaussian, xdata=binscenters, ydata=histo1, p0=[1, np.mean(data), 5],
                                   bounds=[[0, np.min(data), 2], [np.inf, np.max(data), np.inf]],
                                   )
        err_1 = np.sum(np.sqrt(np.diag(pcov_1)))
        # print(f"Fit results 1 gauss: \n\t{popt_1}\n\t{err_1}")

        if err_1 < err_3:
            popt_3[0] = popt_1[0]
            popt_3[1] = 0
            popt_3[2] = popt_1[1]
            popt_3[3] = popt_1[2]
            popt_3[4] = 0
        popt = popt_3

    except Exception as e:
        print("Fitting failed", e)
    return popt

def fitandplot_binneddata(data, label, numbins=20):
    """ Derive RTS paramters for a given pixel based on a sqeuince of readouts. This is based on fitting gaussians to a binned data sample
        Returns a popt set of coefficients .
    """
    start = datetime.datetime.utcnow()
    for ii in range(100):
        popt = fitbinneddata(data, numbins)
    end = datetime.datetime.utcnow()
    print(f"Timing for parameter determination: : {(end - start) / 100}")

    bins = np.linspace(np.min(data), np.max(data), numbins)
    histo1, bins1 = np.histogram(data, bins=bins, )
    histo1 = histo1 / np.max(histo1)
    binscenters = np.array([0.5 * (bins1[i] + bins1[i + 1]) for i in range(len(bins1) - 1)])

    plt.hist(binscenters, len(binscenters), weights=histo1, label=f"data {label}")

    if popt is not None:
        xspace = np.linspace(np.min(data) - 20, np.max(data) + 20, 100)
        plt.plot(xspace, np_tripple_gaussian_function(xspace, *popt),
                 label=f"Fit: {popt[0]: 4.2f} {popt[1]: 4.2f} {popt[2]: 7.1f} "
                       f"{popt[3]:> 4.1f} {popt[4]:> 8.1f}")
    return popt


def main():
    plt.figure()

    data = readdistribution('exampledata/maxdistr.txt')

    popt = fitandplot_binneddata(data, "Worst case")
    data = readdistribution('exampledata/mindistr.txt')
    p_best = fitandplot_binneddata(data, "Best case")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05))
    plt.savefig("distributionfit.png", bbox_inches="tight")
    plt.close()

    # data = readdistribution('exampledata/mindistr.txt')

    plt.figure()
    # popt=p_best
    x = np.linspace(0, 10, num=100)
    y = [summed_ln_likelihood(data[0:], 1, popt[1], popt[2], _x, popt[4]) for _x in x]
    plt.plot(x, y, label="ln likelihood")
    plt.savefig("lhe.png")

    start = datetime.datetime.utcnow()
    for ii in range(100):
        findditributionparamters(data, popt)
    end = datetime.datetime.utcnow()
    print(f"Timing for likelyhodd fit: : : {(end - start) / 100}")

    x = range(len(data))
    means = [np.mean(data[0:x]) for x in range(len(data))]
    mlhres = [findx0(data[0:x], popt[0], popt[1], popt[3], popt[4]) for x in range(len(data))]
    mlh = [res.x[0] for res in mlhres]
    plt.figure()
    plt.plot(x, data, '.')
    plt.plot(x, means, label="Mean of first N")
    plt.hlines(np.mean(data), xmin=0, xmax=len(data), color='black')
    plt.plot(x, mlh, label="maximum likelihood fit")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05))
    plt.savefig("meanvsmle.png", bbox_inches="tight", dpi=300)

    plt.figure()
    plt.plot(x, data, '.')
    plt.plot(x, means, label="Mean of first N")
    plt.hlines(np.mean(data), xmin=0, xmax=len(data), color='black')
    plt.plot(x, mlh, label="maximum likelihood fit")
    plt.ylim([np.mean(data) - 5, np.mean(data) + 5])
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05))
    plt.savefig("meanvsmle_zoom.png", bbox_inches="tight", dpi=300)
    plt.close()

    _data = data[0:]
    start = datetime.datetime.utcnow()
    for ii in range(10000):
        findx0(_data, popt[0], popt[1], popt[3], popt[4])
    end = datetime.datetime.utcnow()
    print(f"Timing for fit: : {(end - start) / 10000}")


if __name__ == '__main__':
    main()
