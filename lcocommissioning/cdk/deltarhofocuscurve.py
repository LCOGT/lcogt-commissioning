import json
import logging
import os

import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import pickle
import matplotlib.patches as patches


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


def overplot_fit(func, paramset, ax):
    if paramset is None:
        return
    base = np.arange(-0.1, 0.1, 0.01)
    y = func(base, *paramset)
    ax.plot(base, y, "--", color='orange' if func == sqrtfit else 'grey',
             label="sqrt {:5.2f}".format(paramset[3]) if func == sqrtfit else "parabola")


def anaylse_deltarho_tilt(bestfits):

    # Get rid of focus zeropint, not relevant for us
    focuszeropoint = bestfits[4]
    for ii in bestfits.keys():
        bestfits[ii] = bestfits[ii] - focuszeropoint

    detectorsizeX = 36
    detectorsizeY = 24

    throwx = 2/3 * detectorsizeX # mm
    throwy = 2/3 * detectorsizeY # mm

    # Sector ids: x increases right, y increases up. Not as plotted!!!
    # 6 7 8
    # 3 4 5
    # 0 1 2

    delta_focus_x = (bestfits[0] + bestfits[3] + bestfits[6]) / 3. - (bestfits[2]+bestfits[5]+bestfits[8]) / 3.
    delta_focus_y = (bestfits[0] + bestfits[1] + bestfits[2]) / 3. - (bestfits[6]+bestfits[7]+bestfits[8]) / 3.

    angle_x = math.atan(delta_focus_x/throwx) # in radians
    angle_y = math.atan(delta_focus_y/throwy)

    shimtrow_x = 137.8 #mm
    shimthrow_y = 220.5 #mm

    correction_x = math.tan(angle_x) * shimtrow_x
    correction_y = math.tan(angle_y) * shimthrow_y
    screwpitch = 0.01 # (mm/rev)



    Narrative ='\n'.join ( (
        f"Focal plane offsets\n X: {delta_focus_x:7.4f} mm left - right\n Y: {delta_focus_y:7.4f} mm bottom - top\n",
        f"Focal plane tilts are along\n  x axis: {angle_x: 7.5f} rad\n  y axis: {angle_y: 7.5f} rad\n",
        f"Shim throw x {shimtrow_x: 5.2f} mm\n -> delta X: {correction_x: 6.4f} mm ",
        f' ??? -> add to {"left of image" if delta_focus_x < 0 else "right of image"}\n',
        f"Shim throw y {shimthrow_y: 5.2f} mm\n -> delta Y: {correction_y: 6.4f} mm",
        f' ??? -> add to {"bottom of image" if delta_focus_y < 0 else "top of image"}',
    ))
    # plot focal plane
    print (Narrative)


    xx=[-throwx/2,throwx/2]
    xy = [-delta_focus_x/2,delta_focus_x/2]
    yx=[-throwy/2,throwy/2]
    yy = [-delta_focus_y/2,delta_focus_y/2]


    fig, ax = plt.subplots()

    DetectorImprint = patches.Rectangle ((-detectorsizeX/2, -detectorsizeY/2), detectorsizeX, detectorsizeY, )

    ul = (-shimtrow_x/2, shimthrow_y/2)
    ur = (+shimtrow_x/2, shimthrow_y/2)
    ll = (-shimtrow_x/2, -shimthrow_y/2)
    lr = (+shimtrow_x/2, -shimthrow_y/2)

    ax.add_patch (DetectorImprint)
    ax.add_patch (patches.Circle (ul))
    ax.add_patch (patches.Circle (ur))
    ax.add_patch (patches.Circle (ll))
    ax.add_patch (patches.Circle (lr))

    ax.set_xlim((-shimtrow_x/2*1.2, shimtrow_x/2*1.2))
    ax.set_ylim((-shimthrow_y/2*1.2, shimthrow_y/2*1.2))

    plt.xlabel ("distance focal plane CCD X [mm]")
    plt.xlabel ("distance focal plane CCD Y [mm]")

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# place a text box in upper left in axes coords
    ax.text(1.05, 0.99, Narrative, transform=ax.transAxes, fontsize=8,
        horizontalalignment='left', verticalalignment='top', bbox=props)

    plt.gca().set_aspect('equal')



def main():

    logging.basicConfig(level=getattr(logging, 'INFO'),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    error_string = None
    focusdict = {}
    fwhmdict = {}

    if '.pickle' in sys.argv[1]:
        analyse_tilts(sys.argv[1])
        exit(0)

    for image in sys.argv[1:]:
        fwhmdict[image] = {}
        focus, fwhm = getImageFWHM(image, minarea=5, sections=True)
        focusdict[image] = focus
        fwhmdict[image] = fwhm
        print (fwhmdict[image])

    bestfits = {}
    plt.figure()
    fig, axes = plt.subplots(3, 3, figsize=(7 ,7))
    plt.subplots_adjust (hspace=0.4)

    for section in range (0,9):
        ax = axes[2 - section // 3, section % 3]
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

        if (return_package is not None) and math.isfinite(bestfocus_error):

            ax.axvline(x=bestfocus, color='orange', label="best focus sqrt")
            ax.set_xlabel("FOCUS  [mm]")
            ax.set_ylabel("FWHM (pix)")
            overplot_fit(sqrtfit, exponential_p, ax=ax)
            ax.plot(focuslist, fwhmlist, 'o')
            ax.set_xlim([-0.1, 0.1])
            ax.set_ylim([0, 8])
            ax.grid(True)
            ax.set_title(f"{section} {bestfocus:5.3f}mm" if math.isfinite(
                      bestfocus_error) else "Fit failed")

    plt.suptitle (f"Focus Gradient {os.path.basename(sys.argv[1])}")
    plt.tight_layout()
    plt.savefig("{}".format("deltarho_focusgradient.pdf"), bbox_inches='tight', dpi=150)

    with open('deltarho_focus' + '.pickle', 'wb') as f:
        pickle.dump(bestfits, f, pickle.HIGHEST_PROTOCOL)

    analyse_tilts('deltarho_focus' + '.pickle')


def analyse_tilts(filename):
    with open(filename, 'rb') as f:
        bestfits = pickle.load(f)
    plt.figure()
    anaylse_deltarho_tilt(bestfits)
    plt.tight_layout()
    plt.savefig('deltarhofocus.pdf')


if __name__ == '__main__':
    main()
