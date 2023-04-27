import json
import logging

import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import pickle

from lcocommissioning.common.SourceCatalogProvider import getImageFWHM

L1FWHM = "L1FWHM"
FOCDMD = "FOCDMD"

_LIMIT_EXPONENT_U = 1.0
_MIN_NUMBER_OF_POINTS = 5

# This describes our model for a focus curve: Seeing and defocus add in quadrature.
sqrtfit = lambda x, seeing, bestfocus, slope, tweak: (seeing ** 2 + (slope * (x - bestfocus)) ** 2) ** tweak
polyfit = lambda x, seeing, bestfocus, slope: slope * (x - bestfocus) ** 2 + seeing


def focus_curve_fit(xdata, ydata, func=sqrtfit):
    """
    Generic iterative fit with sigma rejection.

    :param xdata:
    :param ydata:
    :param func:
    :param plot:
    :return:
    """

    # TODO: Verification that we have enough points.

    # TODO: Boundaries and intial guess externally driven by context.
    initial_guess = [3, 0, 60, 0.5] if func == sqrtfit else None
    bounds = [[2, -0.1, 0, 0.4], [6, 0.1, 100, _LIMIT_EXPONENT_U]] if func == sqrtfit else [-math.inf, math.inf]

    for iter in range(2):
        try:
            (paramset, istat) = optimize.curve_fit(func, xdata, ydata, p0=initial_guess, bounds=bounds)
        except:
            paramset = None
            istat = None

        if paramset is None:
            continue
        # sigma rejection and rms estimate:
        fit = func(xdata, *paramset)
        delta = ydata - fit
        s = np.std(delta)
        good = (delta < 3 * s)
        xdata = xdata[good]
        ydata = ydata[good]

    paramerrors = np.sqrt(np.diag(istat)) if istat is not None else None

    return paramset, paramerrors


def overplot_fit(func, paramset):
    if paramset is None:
        return
    base = np.arange(-0.1, 0.1, 0.01)
    y = func(base, *paramset)
    plt.plot(base, y, "--", color='orange' if func == sqrtfit else 'grey',
             label="sqrt {:5.2f}".format(paramset[3]) if func == sqrtfit else "parabola")


def anaylse_deltarho_tilt(bestfits):

    # Get rid of focus zeropint, not relevant for us
    for ii in bestfits.keys():
        bestfits[ii] = bestfits[ii] - bestfits[4]

    throwx = 2/3 * 36 # mm
    throwy = 2/3 * 24 # mm

    delta_focus_x = (bestfits[0] + bestfits[3] + bestfits[6]) / 3 - (bestfits[2]+bestfits[5]+bestfits[8]) / 3.
    delta_focus_y = (bestfits[0] + bestfits[1] + bestfits[2]) / 3 - (bestfits[6]+bestfits[7]+bestfits[8]) / 3.

    angle_x = math.atan(delta_focus_x/throwx) # in radians
    angle_y = math.atan(delta_focus_y/throwy)

    screwthrowx_from_center = 60 #mm
    screwthrowy_from_center = 60 #mm

    correction_x = math.tan(angle_x) * screwthrowx_from_center
    correction_y = math.tan(angle_y) * screwthrowy_from_center
    screwpitch = 0.01 # (mm/rev)



    print (f"Focal plane offsets X: {delta_focus_x:7.5f} mm Y: {delta_focus_y:7.5f} mm")
    print (f"Focal plane tilts are along x axis: {angle_x:7.5f} rad, along y axis: {angle_y:7.5f} rad")
    print (f"Shim throw x {screwthrowx_from_center:5.2f} mm shim delta X: {correction_x:7.5f} mm ")
    print (f"Shim throw y {screwthrowy_from_center:5.2f} mm shim delta Y: {correction_y:7.5f} mm")
    # plot focal plane
    xx=[-throwx/2,throwx/2]
    xy = [-delta_focus_x/2,delta_focus_x/2]
    yx=[-throwy/2,throwy/2]
    yy = [-delta_focus_y/2,delta_focus_y/2]
    plt.plot (xx,xy,label="x-direction,")
    plt.plot (yx,yy,label="y-direction,")

    plt.plot ([0,screwthrowx_from_center], [0,correction_x], label=f'cor_x={correction_x:5.4f}mm {correction_x/screwpitch:5.4} revs')
    plt.plot ([0,screwthrowy_from_center], [0,correction_y], label=f'cor_y={correction_y:5.4f}mm {correction_y/screwpitch:5.4} revs')


    plt.xlabel ("distance focal plane [mm]")
    plt.ylabel ("defocus [mm]")
    plt.legend()
# Full frame readout is 9576*6388, 3.75um pixels
    # |  |  |  |  throw in x from center is (9576/3)*3.76um = 12.00192 mm



def main():

    logging.basicConfig(level=getattr(logging, 'INFO'),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    error_string = None
    focusdict = {}
    fwhmdict = {}

    if '.pickle' in sys.argv[1]:
        with open(sys.argv[1], 'rb') as f:
            bestfits =  pickle.load(f)
        plt.figure()
        anaylse_deltarho_tilt (bestfits)
        plt.savefig ('deltarhofocus.png')
        exit(0)
    for image in sys.argv[1:]:
        fwhmdict[image] = {}
        focus, fwhm = getImageFWHM(image, minarea=5, sections=True)
        focusdict[image] = focus
        fwhmdict[image] = fwhm
        print (fwhmdict[image])

    bestfits = {}
    plt.figure()
    for section in range (0,9):
        focuslist = []
        fwhmlist = []
        for image in focusdict.keys():
            focuslist.append (focusdict[image])
            fwhmlist.append (fwhmdict[image][section])
        focuslist = np.asarray(focuslist)
        fwhmlist = np.asarray(fwhmlist)
        print ("Focus input: {}\nFWHM: {}".format (np.round(focuslist,3), np.round (fwhmlist,3)))

        print (f"Fitting FWHM in Section {section}")
        exponential_p, exponential_rms = focus_curve_fit(focuslist, fwhmlist, sqrtfit)

        # we will need this a few times - meaningful references here
        if exponential_p is not None:
            bestfocus = exponential_p[1]
            bestfocus_error = exponential_rms[1]
            bestfits[section] = bestfocus

            if not math.isfinite(bestfocus_error):
                error_string = "fit did not converge"
            if bestfocus_error > 0.25:
                error_string = "focus fit is too noisy"
            if abs(exponential_p[1]) > 2.5:
                error_string = "Focus offset too large to be credible."

            return_package = {'fitok': True if error_string is None else False,
                          'fit_seeing': round(exponential_p[0], 2),
                          'fit_focus': round(bestfocus, 2),
                          'fit_slope': round(exponential_p[2], 2),
                          'fit_exponent': round(exponential_p[3], 2),
                          'fit_rms': round(bestfocus_error, 2),
                          'errormsg': error_string}
        else:
            return_package = None

        plt.subplot (3,3,section+1)
        if (return_package is not None) and math.isfinite(bestfocus_error):
            plt.axvline(x=bestfocus, color='orange', label="best focus sqrt")
            #plt.axes().axvspan(bestfocus - bestfocus_error, bestfocus + bestfocus_error, alpha=0.1, color='grey')
            plt.xlabel("FOCUS  [mm foc plane]")
            plt.ylabel("FWHM (pix)")
            overplot_fit(sqrtfit, exponential_p)
            plt.plot(focuslist, fwhmlist, 'o')
            plt.xlim([-0.1, 0.1])
            plt.ylim([0, 6])
            plt.title(f"{section} {bestfocus:5.3f}" if math.isfinite(
                      bestfocus_error) else "Fit failed")
    plt.savefig("{}".format("focus_0.png"), bbox_inches='tight')
    with open('deltarho_focus' + '.pickle', 'wb') as f:
        pickle.dump(bestfits, f, pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    main()
