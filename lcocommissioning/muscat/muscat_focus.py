import argparse
import datetime as dt
import logging
import os
import re

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.io.fits import ImageHDU, CompImageHDU
from lcocommissioning.common.SourceCatalogProvider import SEPSourceCatalogProvider
from lcocommissioning.common.lco_archive_utilities import get_frames_from_request, download_from_archive, \
    get_muscat_focus_requesetids
from lcocommissioning.common.muscatfocusdb_orm import muscatfocusdb, MuscatFocusMeasurement
from scipy import optimize

_log = logging.getLogger(__name__)
logging.getLogger('matplotlib.font_manager').disabled = True
L1FWHM = "L1FWHM"
FOCDMD = "FOCDMD"

_LIMIT_EXPONENT_U = 0.7
_MIN_NUMBER_OF_POINTS = 5

# This describes our model for a focus curve: Seeing and defocus add in quadrature.
sqrtfit = lambda x, seeing, bestfocus, slope, tweak: (seeing ** 2 + (slope * (x - bestfocus)) ** 2) ** tweak
fermifit = lambda x, thetazero, x0, T, a: a / (1 + np.e ** ((x - x0) / T)) + thetazero


def focus_curve_fit(xdata, ydata, func=sqrtfit):
    """
    Generic iterative fit with sigma rejection.

    :param xdata:
    :param ydata:
    :param func:
    :param plot:
    :return:
    """

    # TODO: Boundaries and initial guess externally driven by context.
    initial_guess = None
    bound = None
    label = None

    if func == sqrtfit:
        initial_guess = [2, 0, 1, 0.6]
        bounds = [[0, -3, 0, 0.5], [5, 3, 5, _LIMIT_EXPONENT_U]]

    if func == fermifit:
        initial_guess = [0, 0, 100, 0]
        ydata = ydata * np.pi / 180.
        bounds = [[-np.pi, -np.inf, 0, -np.pi], [+np.pi, np.inf, np.inf, np.pi]]

    for iter in range(2):
        try:
            (paramset, istat) = optimize.curve_fit(func, xdata, ydata, p0=initial_guess, bounds=bounds)
        except:
            paramset = None
            istat = None
            _log.debug("fit error")

        if paramset is None:
            continue
        # sigma rejection and rms estimate:
        fit = func(xdata, *paramset)
        delta = ydata - fit
        s = np.std(delta)
        good = (delta < 2 * s)
        xdata = xdata[good]
        ydata = ydata[good]

    paramerrors = np.sqrt(np.diag(istat)) if istat is not None else None
    if (func == fermifit) & (paramset is not None):
        paramset[3] *= 180./np.pi
        paramset[0] *= 180./np.pi

    return paramset, paramerrors


def overplot_fit(func, paramset, color=None):
    if paramset is None:
        _log.info ("Not plotting function since plot parameters are None")
        return
    base = np.arange(-3.6, 3.6, 0.1)
    y = func(base, *paramset)
    if color is None:
        color = 'blue' if func == sqrtfit else 'grey'
    plt.plot(base, y, "--", color=color,
             label="sqrt {:5.2f}".format(paramset[3]) if func == sqrtfit else "parabola")


def getImageData(imagename, minarea=20, deblend=0.5, archive=False):
    """ Measure the FWHM of an image, tuned to get a reasonable FWHM also for defocussed images.
    """
    _log.info (f"Getting image datga for image {imagename}")
    if archive:
        hdul = download_from_archive (imagename)
    else:
        hdul = fits.open(imagename, 'readonly', ignore_missing_end=True)

    deltaFocus = None
    pixelscale = None
    foctemp = None
    dateobs = None
    for ii in range(len(hdul)):
        if FOCDMD in hdul[ii].header:
            deltaFocus = hdul[ii].header[FOCDMD]
        if 'PIXSCALE' in hdul[ii].header:
            pixelscale = hdul[ii].header['PIXSCALE']
        if 'FOCTEMP' in hdul[ii].header:
            foctemp = hdul[ii].header['FOCTEMP']
        if 'DATE-OBS' in hdul[ii].header:
            dateobs = hdul[ii].header['DATE-OBS']

    _log.info (f"Delta Focus:  {deltaFocus}")
    catalog = SEPSourceCatalogProvider(refineWCSViaLCO=False)
    fwhmcat = np.asarray([])
    thetacat = np.asarray([])
    ellcat = np.asarray([])
    for ii in range(len(hdul)):
        if isinstance(hdul[ii], ImageHDU) or isinstance(hdul[ii], CompImageHDU):
            cat, wcs = catalog.get_source_catalog_from_fitsobject(hdul, ext=ii, minarea=minarea, deblend=deblend)
            fwhmcat = np.append(fwhmcat, cat['fwhm'])
            thetacat = np.append(thetacat, cat['theta'])
            ellcat = np.append(thetacat, cat['ellipticity'])
    hdul.close()

    # comprehension of the object catalog....
    good = fwhmcat > 0
    goodtheta = thetacat > -720
    goodell = (ellcat > 0) & (ellcat <= 1)

    for iter in range(3):
        medianfwhm = np.median(fwhmcat[good])
        mediantheta = np.median(thetacat[goodtheta])
        medianell = np.median(ellcat[goodell])
        fwhmstd = np.std(fwhmcat[good])
        thetastd = np.std(thetacat[goodtheta])
        ellstd = np.std(ellcat[goodell])

        good = abs(fwhmcat - medianfwhm) < 2 * fwhmstd
        goodtheta = abs(thetacat - mediantheta) < 2 * thetastd
        goodell = abs(ellcat - medianell) < 2 * ellstd

        if np.sum(good) > 5:
            medianfwhm = np.median(fwhmcat[good])
        if np.sum(goodtheta) > 5:
            mediantheta = np.median(thetacat[good])
        if np.sum(goodell) > 5:
            medianell = np.median(ellcat[goodell])

    if pixelscale is None:
        pixelscale = 1
    medianfwhm *= pixelscale

    _log.info ("{}  FOCCMD {: 5.3f} FWHM (\" : pix) ({: 5.2f} : {: 5.2f}) \pm {: 5.2f} pixel  {: 5.2f} {: 6.4f}".format(imagename, deltaFocus,
                                                                                                                        medianfwhm, medianfwhm / pixelscale,
                                                                                                                        fwhmstd, mediantheta, pixelscale))



    return deltaFocus, medianfwhm, mediantheta, medianell, foctemp, dateobs


def sort_input_images (inputlist):
    ''' Sort the input image list by camera. Return a dictionary {cameraname -> [listof iamges]}'''
    cameraregex = '.*[012]m.*-([feas][apfkq]\d\d)-.*'

    returndict = {}

    for image in inputlist:
        m = re.search (cameraregex, image)
        if (m is not None):
            camera = m.group(1)
            if camera not in returndict:
                returndict[camera] = []
            returndict[camera].append(image)

        else:
            _log.debug (f"{image}: No camera identifier found, skipping!")
    return returndict

def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Do a simultaneous focus fit of Muscat ep cameras f')


    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument ('--files', type=str, nargs='+')
    group.add_argument ('--requestid', type=int, nargs='?')
    group.add_argument ('--crawl-after', type=dt.datetime.fromisoformat, help="Consider data only after the given date, in ISO format (e.g., 2022-09-01 00:00:00")

    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')
    parser.add_argument('--noplot', action='store_true',)
    parser.add_argument('--before',  type=dt.datetime.fromisoformat, default=dt.datetime.utcnow().isoformat())
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args

def get_auto_focus_frames(requestid):
    candidates = get_frames_from_request(requestid)
    # get the camera ids
    cameraregex = '.*[012]m.*-([sfea][aqpfk]\d\d)-.*'
    cameras = []
    for imageinfo in candidates['results']:
        m = re.search (cameraregex, imageinfo['basename'])
        if (m is not None):
            cameras.append(m.group(1))
    camaras = sorted (set(cameras))
    _log.info (f"camera set detected: {camaras}")
    focusimagedict = {}
    for camera in cameras:
        focusimagedict[camera]=[]


    for imageinfo in candidates['results']:
        if 'x00' in imageinfo['basename']:
            m = re.search (cameraregex, imageinfo['basename'])
            camera =   m.group(1)  if (m is not None) else None
            if camera is not None:
                focusimagedict[camera].append(imageinfo['id'])

    return focusimagedict

def get_focus_measurements(args, inputimagedict):
    measurementlist = {}
    for camera in inputimagedict:
        measurementlist[camera] = {'focuslist': [], 'fwhmlist': []}

        for image in inputimagedict[camera]:
            focus, fwhm, theta, ell, foctemp, dateobs  = getImageData(image, minarea=5, deblend=0.5, archive=(args.requestid is not None) or (args.crawl_after is not None))
            measurementlist[camera]['focuslist'].append(focus)
            measurementlist[camera]['fwhmlist'].append(fwhm)

        measurementlist[camera]['focuslist'] = np.asarray(measurementlist[camera]['focuslist'])
        measurementlist[camera]['fwhmlist'] = np.asarray(measurementlist[camera]['fwhmlist'])

        exponential_p, exponential_rms = focus_curve_fit(measurementlist[camera]['focuslist'],
                                                         measurementlist[camera]['fwhmlist'], sqrtfit)

        measurementlist[camera]['exponential_p'] = exponential_p
        measurementlist[camera]['exponential_rms'] = exponential_p

    return measurementlist, foctemp, dateobs

def get_focus_measurements_requestid (requestid, args):
    inputimagedict = get_auto_focus_frames(requestid)
    measurementlist, foctemp, dateobs  = get_focus_measurements(args, inputimagedict)
    return measurementlist, foctemp, dateobs

def get_focusmeasurements_filelist(filelist, args):
    inputimagedict = sort_input_images(filelist)
    measurementlist, foctemp, dateobs  = get_focus_measurements(args, inputimagedict)
    return measurementlist, foctemp, dateobs


def  process_single_requestid (requestid, args):
    measurementlist, foctemp, dateobs  = get_focus_measurements_requestid(requestid, args)
    camera_g = 'ep06' if args.muscat=='mc04' else 'ep04'
    camera_r = 'ep07'  if args.muscat=='mc04' else 'ep02'
    camera_i = 'ep08'  if args.muscat=='mc04' else 'ep03'
    camera_z = 'ep09'  if args.muscat=='mc04' else 'ep05'

    _log.info (f"Measureemnt results: {measurementlist}")
    goodresult = measurementlist[camera_g]['exponential_p'] is not None
    goodresult = goodresult and (measurementlist[camera_r]['exponential_p'] is not None)
    goodresult = goodresult and (measurementlist[camera_i]['exponential_p'] is not None)
    goodresult = goodresult and (measurementlist[camera_z]['exponential_p'] is not None)

    if goodresult:
        newitem = MuscatFocusMeasurement(requestid=int (args.requestid),
                                     muscat = str (args.muscat),
                                     dateobs=str(dateobs),
                                     temperature = float(foctemp),
                                     camera_g = camera_g,
                                     focus_g  = float(measurementlist[camera_g]['exponential_p'][1]),
                                     seeing_g = float(measurementlist[camera_g]['exponential_p'][0]),
                                     error_g  = float(measurementlist[camera_g]['exponential_rms'][1]),

                                     camera_r = camera_r,
                                     focus_r  = float(measurementlist[camera_r]['exponential_p'][1]),
                                     seeing_r = float(measurementlist[camera_r]['exponential_p'][0]),
                                     error_r  = float(measurementlist[camera_r]['exponential_rms'][1]),

                                     camera_i = camera_i,
                                     focus_i  = float(measurementlist[camera_i]['exponential_p'][1]),
                                     seeing_i = float(measurementlist[camera_i]['exponential_p'][0]),
                                     error_i  = float(measurementlist[camera_i]['exponential_rms'][1]),

                                     camera_z = camera_z,
                                     focus_z  = float(measurementlist[camera_z]['exponential_p'][1]),
                                     seeing_z = float(measurementlist[camera_z]['exponential_p'][0]),
                                     error_z  = float(measurementlist[camera_z]['exponential_rms'][1]),
                                     )
    else:
        newitem = MuscatFocusMeasurement(requestid=int (args.requestid))
    _database.addMeasurement(newitem)
    if not args.noplot:
        plot_focuscurve(measurementlist, args)


def main():
    args = parseCommandLine()
    args.muscat = 'mc04'
    if (args.files):
        measurementlist, foctemp, dateobs = get_focusmeasurements_filelist(args.files, args)
        if not args.noplot:
            plot_focuscurve(measurementlist, args)
    if args.requestid:
        process_single_requestid(args.requestid, args)


    if args.crawl_after:
        requestsids = get_muscat_focus_requesetids(args.muscat)
        _log.info(f"Got {len(requestsids)} auto focus entries")
        for requestid in requestsids:
            _log.info (f" PROCESSING {requestid}")
            args.requestid = requestid
            if  _database.exists(requestid) is None:
                process_single_requestid(requestid, args)
            else:
                _log.info (f"Ommitting request file {requestid} since already processed.")




def plot_focuscurve(measurementlist, args):
    plt.style.use('ggplot')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.xlim([-3.6, 3.6])
    plt.ylim([0, 4])
    plt.xlabel("FOCUS Demand [mm foc plane]")
    plt.ylabel("FWHM ['']")
    plothandles = []
    bestfocus_yu = 0.1
    colors = ['g', 'r', 'b', 'orange']
    coloridx = 0
    for camera in measurementlist:

        try:

            color = colors[coloridx % len(colors)]
            coloridx = coloridx + 1
            bestfocus = measurementlist[camera]['exponential_p'][1]
            bestfocus_error = measurementlist[camera]['exponential_rms'][1]
            overplot_fit(sqrtfit, measurementlist[camera]['exponential_p'], color=color)
            plt1, = plt.plot(measurementlist[camera]['focuslist'], measurementlist[camera]['fwhmlist'], 'o',
                             color=color, label=f"{camera} {bestfocus:6.2f}")

            plt.axvline(x=bestfocus, color=color)
            plothandles.append(plt1)
            print(f"Best focus for camera {camera}: {bestfocus:5.2f}")
        except:
            _log.error("Something bad")
    ax1.legend(handles=plothandles, loc="lower right", bbox_to_anchor=(1, -0.1),
               bbox_transform=fig.transFigure, ncol=4)
    plt.title(f'Muscat Focus ID {args.requestid if args.requestid is not None else ""}')
    plt.savefig("{}".format(
        f'muscat_focus_{args.requestid if args.requestid is not None else os.path.basename(args.files[0])}.png'),
                bbox_inches="tight")

_database = muscatfocusdb('sqlite:///muscatfocus.sqlite')
if __name__ == '__main__':


    main()
    _database.close()
