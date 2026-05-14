"""
Program to calculate the noise and gain for each extension of a given mef file
"""

import logging
import sys
import os
import os.path
import numpy as np
import argparse
import astropy.time as astt
import astropy.stats
from astropy.stats import sigma_clip
from astropy.table import Table
from astropy.io import ascii
from astropy.modeling import models, fitting
import math
import json

from lcocommissioning.common import lco_archive_utilities
from lcocommissioning.common.ccd_noisegain import dosingleLevelGain
from lcocommissioning.common.noisegaindb_orm import NoiseGainMeasurement, noisegaindb
from lcocommissioning.common.Image import Image
from lcocommissioning.common.logging_config import setup_logging
import matplotlib.pyplot as plt
from astropy.io import fits

_logger = logging.getLogger(__name__)
mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)


def findkeywordinhdul(hdulist, keyword):
    # f&*& fpack!
    for ext in hdulist:
        val = ext.header.get(keyword)
        if val is not None:
            return val

    return None


def sortinputfitsfiles(
    listoffiles,
    sortby="exptime",
    selectedreadmode="full_frame",
    ignoretemp=False,
    useaws=False,
    useoverscan=True,
):
    """go through the list of input and sort files by bias and flat type.
    Find pairs of flats that are of the same exposure time / closest in illumination level.
    """

    sortedlistofFiles = {}
    filemetrics = {}
    listoffiles.sort(['FILENAME'])
    # random.shuffle(listoffiles)
    for filecandidate in listoffiles:
        filename = str(filecandidate["FILENAME"])
        # First stage: go through the images and derive the metrics from them to pair.
        # TODO: avoid opening all biases, it is pointless, need only 2!

        if ("bias" in sortedlistofFiles):
            if (len(sortedlistofFiles["bias"]) == 2) and ("b00" in filename):
                _logger.debug(f'Already have two bias frames, skipping further bias candidate {filename}.')
                continue              
                
        if useaws:
            hdu = lco_archive_utilities.download_from_archive(filecandidate["frameid"])
        else:
            _logger.info (f"Candidate:  {filename}")
            filecandidate = {"FILENAME": filename}
            fitsfilepath = str(filecandidate["FILENAME"])
            hdu = fits.open(fitsfilepath)

        # Note: The values below could come from the elasticsearch already. But not if working on local files.
        ccdstemp = findkeywordinhdul(hdu, "CCDSTEMP")
        ccdatemp = findkeywordinhdul(hdu, "CCDATEMP")
        readoutmode = findkeywordinhdul(hdu, "CONFMODE")
        if readoutmode is None:
            readoutmode = str(findkeywordinhdul(hdu, "READMODE"))

        if readoutmode not in selectedreadmode:
            _logger.info(
                "Rejecting file as it is not in the correct readout mode ({} != {})".format(
                    readoutmode, selectedreadmode
                )
            )
            hdu.close()
            continue

        tempdiff = 0
        if (ccdstemp is not None) & (ccdatemp is not None) & (not ignoretemp):
            tempdiff = float(ccdatemp) - float(ccdstemp)

        if abs(tempdiff) > 2:
            hdu.close()
            _logger.info(
                "rejecting file {}: CCD temp is not near set point, delta = {:5.2f}".format(
                    filecandidate, tempdiff
                )
            )
            continue

        if "b00" in filecandidate["FILENAME"]:  # it is a bias
            if "bias" not in sortedlistofFiles:
                sortedlistofFiles["bias"] = []
            if len(sortedlistofFiles["bias"]) < 2:
                sortedlistofFiles["bias"].append(filename)

        else:  # it is (interpreted as) a flat
            if sortby == "exptime":
                exptime = findkeywordinhdul(hdu, "EXPTIME")
                if exptime is not None:
                    filemetrics[str(filecandidate["FILENAME"])] = str(exptime)

            if sortby == "filterlevel":
                filter = findkeywordinhdul(hdu, "FILTER")
               

                if (filter is not None) and ("b00" not in filename):
                    image = Image(
                        hdu, overscancorrect=useoverscan, alreadyopenedhdu=True
                    )
                    if image.data is None:
                        level = -1
                    else:
                        naxis1 = image.data[0].shape[1]
                        naxis2 = image.data[0].shape[0]
                        level = np.nanmean(
                            image.data[0][
                                naxis2 // 4 : naxis2 * 3 // 4,
                                naxis1 // 4 : naxis1 * 3 // 4,
                            ]
                        )
                        if (level == 0):
                            level = -1
                        
                    _logger.info(
                        f'Input file metrics {filename} filter:{filter} light level: {level: 8.1f} naxis 1/2: {naxis1} {naxis2} overscan corrected: {useoverscan}'
                    )
                    filemetrics[filename] = (filter, level)

        hdu.close()

    if "bias" not in sortedlistofFiles:
        _logger.fatal("No suitable bias frames found in list!")
        return sortedlistofFiles
    # pair the flat fields
    if sortby == "exptime":
        # task is simple: just find flats with the same exposure time
        unique = set(filemetrics.values())

        for et in sorted(unique, key=float):
            if float(et) > 0.001:
                sortedlistofFiles[str(et)] = []

        for filename in filemetrics.keys():
            if ("x00" in filename) or ("f00" in filename):
                # identified a bias exposure
                # print ("%s, %s" % (filename, filemetrics[filename]))
                if len(sortedlistofFiles[filemetrics[filename]]) < 2:
                    sortedlistofFiles[filemetrics[filename]].append(filename)

    if sortby == "filterlevel":
        tempsortedListofFiles = {}

        for filename in filemetrics.keys():
            (filter, level) = filemetrics[filename]
            if level < 5:
                _logger.info("rejecting image {}, level is to low".format(filename))

            if filter not in tempsortedListofFiles.keys():
                tempsortedListofFiles[filter] = {}
                ## a new level detected
                tempsortedListofFiles[filter][level] = [
                    filename,
                ]
                _logger.debug(
                    "Starting first entry for new filter: %s %s %s"
                    % (filename, filter, level)
                )
                continue

            matchfound = False
            for knownfilter in tempsortedListofFiles.keys():
                if filter == knownfilter:
                    for knownlevel in tempsortedListofFiles[knownfilter].keys():
                        if (0.95 < knownlevel / level) and (knownlevel / level < 1.05):

                            if len(tempsortedListofFiles[knownfilter][knownlevel]) < 2:
                                _logger.debug(
                                    "adding to set: %s level of %f is within 5 percent of level %f"
                                    % (filename, level, knownlevel)
                                )
                                tempsortedListofFiles[knownfilter][knownlevel].append(
                                    filename
                                )
                                matchfound = True
                                continue

            if not matchfound:
                tempsortedListofFiles[filter][level] = [
                    filename,
                ]
                _logger.debug("Starting new level pair with file %s " % (filename))

        for knownfilter in tempsortedListofFiles.keys():
            for knownlevel in tempsortedListofFiles[knownfilter].keys():
                sortedlistofFiles["%s% 8f" % (knownfilter, knownlevel)] = (
                    tempsortedListofFiles[knownfilter][knownlevel]
                )

    _logger.debug(sortedlistofFiles)
    return sortedlistofFiles

def find_nearest(array, value):
    delta = np.abs(array - value)
    idx = delta.argmin()
    print ("Index " ,  idx, delta)
    return idx

def graphresults(
    alllevels,
    allgains,
    allnoises,
    allshotnoises,
    allexptimes,
    alldateobs,
    args,
    maxlinearity=160000,
):

    adurange = 1 << args.adubits

    # plot flux vs dateobs to check for any temporal trends in the data
    plt.figure()
    for ext in alllevels:
        myexptimes = np.asarray(allexptimes[ext])
        flux = np.asarray(alllevels[ext]) / myexptimes
        mydateobs = np.asarray(alldateobs[ext])
        plt.plot(mydateobs, flux, "o", label="extension %s data" % (ext))
        plt.plot(
            mydateobs[myexptimes == 3],
            flux[myexptimes == 3],
            "o",
            label="extension %s data[Texp=3s]" % (ext),
        )
    plt.legend()
    plt.ylabel(("Flux [ADU/s]"))
    plt.xlabel("DATE-OBS")
    plt.gcf().autofmt_xdate()
    plt.title(args.readmode)
    plt.savefig(f"ptc_{args.readmode}_dateobs_flux.png", bbox_inches="tight")
    plt.close()

    # plot flux vs level to check for any trends in the data 
    plt.figure()
    for ext in alllevels:
        flux = np.asarray(alllevels[ext]) / np.asarray(allexptimes[ext])
        level = np.asarray(alllevels[ext])
        print(level, flux)
        plt.plot(level, flux, "o", label="extension %s data" % (ext))
    plt.legend()
    plt.ylabel(("Flux [ADU/s]"))
    plt.xlabel("Level  [ADU]")
    plt.title(f"Flux vs Level - {args.readmode}")
    plt.savefig(f"ptc_{args.readmode}_level_flux.png", bbox_inches="tight")
    plt.close()

    
    # plot gain vs level to check for any trends in the data and identify the linearity regime
    plt.figure()
    f, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 1]})
    for ext in alllevels:
        gains = np.asarray(allgains[ext])
        levels = np.asarray(alllevels[ext])
        sortedlevels = np.sort(levels)
        noises = np.asarray(allnoises[ext])
        good = (levels >100) & (levels < maxlinearity) & (gains < 20)
        gaindata = gains[(levels > 100) & (levels < maxlinearity) & (gains < 20)]
        noisedata = noises[(levels > 100) & (levels < maxlinearity) & (gains < 20)] 
        bestgain, med, gainstddev = astropy.stats.sigma_clipped_stats(gaindata, sigma=3, maxiters=3)
        bestnoise, med, noisestddev = astropy.stats.sigma_clipped_stats(noisedata, sigma=3, maxiters=3)
      
        fit = fitting.LinearLSQFitter()
        or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip, niter=3, sigma=3.0)
        line_init = models.Linear1D()
        fitted_gain_line, mask = or_fit(line_init, levels[good], gains[good], )
        filtered_levels = np.ma.masked_array(levels[good], mask=mask)
        filtered_gains = np.ma.masked_array(gains[good], mask=mask)
        bestgain = fitted_gain_line.intercept.value

        ax1.plot(levels, gains, "ko", fillstyle = "none", )
        ax1.plot(filtered_levels, filtered_gains, "ko", label=f"extension {ext} Gain" )
        ax1.plot(sortedlevels, fitted_gain_line(sortedlevels), "-", label=f"extension {ext} fit: gain = {bestgain:.2f} e-/ADU")

        line_init = models.Linear1D()
        fitted_noise_line, mask = or_fit(line_init, levels[good], noises[good], )
        filtered_levels = np.ma.masked_array(levels[good], mask=mask)
        filtered_noises = np.ma.masked_array(noises[good], mask=mask)
        bestnoise = fitted_noise_line.intercept.value
                
        ax2.plot(levels, noises, "ko", )
        ax2.plot(filtered_levels, filtered_noises, "ko", label=f"extension {ext} Noise" )
        ax2.plot(sortedlevels, fitted_noise_line(sortedlevels), "-", label=f"extension {ext} fit: noise = {bestnoise:.2f} e-")

        print(f"Best gain for ext {ext}: {bestgain:.2f} e-/ADU, best noise: {bestnoise:.2f} e-")

    ax1.set_ylim([0, 7])
    ax1.set_xlim([1, 1.1 * adurange])
    ax1.legend()
    
    ax1.set_ylabel("Gain [e-/ADU]")
    ax1.set_title(f"Gain & Readnoise {args.readmode}")
    
    ax2.set_ylim([0, 20])
    ax2.set_xlim([1, 1.1 * adurange])
    ax2.legend()
    ax2.set_ylabel("Readnoise [e-]")
    ax2.set_xlabel(("Exposure level [ADU]"))

    plt.savefig(f"ptc_{args.readmode}_levelgain.png", bbox_inches="tight")
    plt.close()

    print(adurange)
    _logger.debug("Plotting ptc")
    plt.figure()
    for ext in alllevels:
        plt.loglog(alllevels[ext], allshotnoises[ext], ".", label="extension %s" % ext)
    plt.legend()
    plt.xlim([1, 1.1 * adurange])
    plt.ylim([1, 3 * math.sqrt(adurange)])
    plt.xlabel("Exposure Level [ADU]")
    plt.ylabel("Measured Noise [ADU]")
    plt.title(args.readmode)
    plt.savefig(f"ptc_{args.readmode}_ptc.png", bbox_inches="tight")
    plt.close()


    _logger.info("Plotting level vs exptime")
    plt.figure()
    f, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={"height_ratios": [2, 1]})
    for ext in alllevels:

        exptimes = np.asarray(allexptimes[ext])
        levels = np.asarray(alllevels[ext])
        texp_sorted = np.sort(exptimes)     
        good = (exptimes >= 2) &(levels > 0)

        fit = fitting.LinearLSQFitter()
        or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip, niter=3, sigma=3.0)
        line_init = models.Linear1D()
        fitted_line, mask = or_fit(line_init, exptimes[good], levels[good], )
        filtered_texp = np.ma.masked_array(exptimes[good], mask=mask)
        filtered_levels = np.ma.masked_array(levels[good], mask=mask)

        ax1.plot(exptimes, levels, "ko", fillstyle = "none",)
        ax1.plot(filtered_texp, filtered_levels, "ko", label=f"extension {ext}")
        ax1.plot(texp_sorted, fitted_line(texp_sorted), "-", label=f"fit: {fitted_line.slope.value:.2f} ADU/s *texp + {fitted_line.intercept.value:.2f} ADU")

        midlevel = np.max (filtered_levels) / 2.
        mp_idx = max (0, find_nearest(levels, midlevel)-1)
        print ("Midlevel",  mp_idx,levels[mp_idx]) 
        mp_level = levels[mp_idx]
        mp_texp = exptimes[mp_idx]
        mp_flux = mp_level / mp_texp
        #linearity residual
        LR = 100 * (1 - mp_flux / ( levels / exptimes ) )

        ax2.plot(
            levels[good], LR[good], ".", label=f"ext {ext}"
        )

    ax1.legend()
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position("bottom")
    ax1.set_xlabel("Exposure time [s]")
    ax1.set_ylabel("Exposure level [ADU]")
    ax2.set_xlabel("Exposure Level [ADU]")
    ax2.set_ylabel("Residual (%)")

    plt.suptitle(f"{args.readmode}\n")
    # ax1.set_ylim([0, 65000])
    # ax2.set_xlim([0, 65000])
    ax2.set_ylim([-2,2])
    plt.savefig(f"ptc_{args.readmode}_texplevel.png", bbox_inches="tight")
    plt.close()


def frameidfromname(fname, filelist):
    """Tool to look up a lco archive frame id for a filename."""
    return filelist[filelist["FILENAME"] == fname]["frameid"][0]


def do_noisegain_for_fileset(
    inputlist, database: noisegaindb, args, frameidtranslationtable=None
):
    """Go through a list of files (bias and flats) to measure noise and gain
    on them. Optionally store results into a database backend make a nice graph.

    :param inputlist: List full path to flat, bias images.
    :param database: Database storage backend
    :param args:
    :param frameidtranslationtable: a table containing the 'FILENAME' and 'frameid' columns to translate
                                    image path into a lco archive request. Used if the args.useaws flag is set,
                                    in which case the images are downloaded via the archive API and not read from disk.
    :return:
    """
    alllevels = {}
    allgains = {}
    allnoises = {}
    allshotnoises = {}
    alllevel1s = {}
    alllevel2s = {}
    allexptimes = {}
    alldateobs = {}

    _logger.debug(
        f"Sifting through the {len(inputlist)} input files of readmode {args.readmode} and finding viable flat pair candidates"
    )
    sortedinputlist = sortinputfitsfiles(
        inputlist,
        sortby=args.sortby,
        selectedreadmode=args.readmode,
        ignoretemp=args.ignoretemp,
        useaws=args.useaws,
        useoverscan=not args.ignoreov,
    )
    _logger.info(
        f"Found {len(sortedinputlist)} viable sets for readmode {args.readmode}.")
    _logger.debug(f"Sorted input list: {json.dumps(sortedinputlist, indent=4)}")
    

    if 'bias' in sortedinputlist and len (sortedinputlist["bias"]) >1:
        bias1_fname = sortedinputlist["bias"][0]
        bias2_fname = sortedinputlist["bias"][1]
    else:
        _logger.warning(f'Not enough bias frames found to perform noise/gain measurement, skipping')
        return
    
    if frameidtranslationtable is not None:
        bias1_frameid = frameidfromname(
            sortedinputlist["bias"][0], frameidtranslationtable
        )
        bias2_frameid = frameidfromname(
            sortedinputlist["bias"][1], frameidtranslationtable
        )
        _logger.info(f"Bias1 id {bias1_fname} -> {bias1_frameid}")
        _logger.info(f"Bias2 id {bias2_fname} -> {bias2_frameid}")

    bias1 = (
        fits.open(bias1_fname)
        if not args.useaws
        else lco_archive_utilities.download_from_archive(bias1_frameid)
    )
    bias2 = (
        fits.open(bias2_fname)
        if not args.useaws
        else lco_archive_utilities.download_from_archive(bias2_frameid)
    )

    for pair_ii in sortedinputlist:

        if "bias" not in pair_ii:
            if len(sortedinputlist[pair_ii]) == 2:
                flat_1_fname = sortedinputlist[pair_ii][0]
                flat_2_fname = sortedinputlist[pair_ii][1]
                _logger.debug(f"\nNoise / Gain measurement based on metric {pair_ii}")
                _logger.debug(f" Flat names {os.path.basename(flat_1_fname)} {os.path.basename(flat_2_fname)}"
                )
                _logger.debug(
                    f" Bias names {os.path.basename(bias1_fname)} {os.path.basename(bias2_fname)}"
                )

                flat1 = (
                    fits.open(flat_1_fname)
                    if not args.useaws
                    else lco_archive_utilities.download_from_archive(
                        frameidfromname(flat_1_fname, frameidtranslationtable)
                    )
                )
                flat2 = (
                    fits.open(flat_2_fname)
                    if not args.useaws
                    else lco_archive_utilities.download_from_archive(
                        frameidfromname(flat_2_fname, frameidtranslationtable)
                    )
                )

                gains, levels, noises, shotnoises, level1s, level2s, exptimes = (
                    dosingleLevelGain(
                        bias1,
                        bias2,
                        flat1,
                        flat2,
                        args,
                        overscancorrect=not args.ignoreov,
                    )
                )

                # grabbing some meta data while we can
                dateobs = findkeywordinhdul(flat1, "DATE-OBS")
                camera = findkeywordinhdul(flat1, "INSTRUME")
                filter = findkeywordinhdul(flat1, "FILTER")
                readmode = findkeywordinhdul(flat1, "CONFMODE")
                dateobst = astt.Time(dateobs, scale="utc", format=None).to_datetime()
                flat1.close()
                flat2.close()

                for extension in range(len(levels)):
                    if extension not in alllevels:
                        alllevels[extension] = []
                        allgains[extension] = []
                        allnoises[extension] = []
                        allshotnoises[extension] = []
                        alllevel1s[extension] = []
                        alllevel2s[extension] = []
                        allexptimes[extension] = []
                        alldateobs[extension] = []

                    if database is not None:
                        identifier = "%s-%s-%s" % (
                            os.path.basename(sortedinputlist[pair_ii][0]),
                            os.path.basename(sortedinputlist[pair_ii][1]),
                            extension,
                        )
                        m = NoiseGainMeasurement(
                            name=identifier,
                            dateobs=dateobs,
                            camera=camera,
                            filter=filter,
                            extension=extension,
                            gain=float(gains[extension]),
                            readnoise=float(noises[extension]),
                            level=float(levels[extension]),
                            differencenoise=float(shotnoises[extension]),
                            level1=float(level1s[extension]),
                            level2=float(level2s[extension]),
                            readmode=readmode,
                        )
                        database.addMeasurement(m)
                        _logger.info(f"Added to database: {m}")

                    alllevels[extension].append(levels[extension])
                    allgains[extension].append(gains[extension])
                    allnoises[extension].append(noises[extension])
                    allshotnoises[extension].append(shotnoises[extension])
                    alllevel1s[extension].append(level1s[extension])
                    alllevel2s[extension].append(level2s[extension])
                    allexptimes[extension].append(exptimes[extension])
                    alldateobs[extension].append(dateobst)

    bias1.close()
    bias2.close()

    for ext in alllevels:
        t = Table()
        t["avglevel"] = alllevels[ext]
        t["level1"] = alllevel1s[ext]
        t["level2"] = alllevel1s[ext]
        t["gain"] = allgains[ext]
        t["noise"] = allnoises[ext]
        t["flatnoise"] = allshotnoises[ext]
        t["exptime"] = allexptimes[ext]
        t["dateobs"] = alldateobs[ext]
        try:
            ascii.write(t, f"ptc_data_{ext}.dat", overwrite=True)
        except Exception as e:
            _logger.debug(f"Error writing PTC data for extension {ext}: {e}")

    if args.makepng:
        graphresults(
            alllevels, allgains, allnoises, allshotnoises, allexptimes, alldateobs, args
        )


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description="General purpose CCD noise and gain measurement from pairs of flat fields and biases.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "fitsfile",
        type=str,
        nargs="+",
        help="Input fits files, must include at least two bias and two flat field frames.",
    )

    group = parser.add_argument_group(
        "Optionally, specify the location of statistics window. All units in pixels with an FITS image extension"
    )
    group.add_argument("--minx", type=int, default=None, help="minimum x.")
    group.add_argument("--maxx", type=int, default=None, help="maximum x.")
    group.add_argument("--miny", type=int, default=None, help="miniumm y.")
    group.add_argument("--maxy", type=int, default=None, help="maximum y.")

    parser.add_argument(
        "--imagepath",
        dest="opt_imagepath",
        type=str,
        default=None,
        help="pathname to prepend to fits file names.",
    )
    parser.add_argument(
        "--database",
        default="sqlite:///noisegain.sqlite",
        help="sqlite database where to store results.",
    )
    parser.add_argument(
        "--sortby",
        type=str,
        default="exptime",
        choices=["exptime", "filterlevel"],
        help="Automatically group flat fiel;ds by exposure time (great if using dome flas, or lab flats)."
        ", or by measured light level (great when using sky flats, but more computing intensive",
    )
    parser.add_argument("--readmode", default="full_frame")
    parser.add_argument("--adubits", default=18)
    parser.add_argument(
        "--noreprocessing",
        action="store_true",
        help="Do not reprocess if data are already in database",
    )
    parser.add_argument(
        "--ignoretemp",
        action="store_true",
        help="ignore if actual temperature differs from set point temperature. Reject by default.",
    )
    parser.add_argument(
        "--ignoreov", action="store_true", help="ignore overscan definition"
    )

    parser.add_argument(
        "--loglevel",
        dest="log_level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARN"],
        help="Set the debug level",
    )
    parser.add_argument(
        "--showimages",
        action="store_true",
        help="Interactively show difference flat and bias images.",
    )
    parser.add_argument(
        "--makepng",
        action="store_true",
        help="Create a png output image of noise, gain, and ptc.",
    )

    parser.add_argument(
        "--texpdelta",
        default = 0.0,
        type=float,
        help="Fixed value to add to Texp to compensatate for shutter mismatch'"
    )

    args = parser.parse_args()
    args.useaws = False

    setup_logging(args.log_level)

    for ii in range(len(args.fitsfile)):
        if args.opt_imagepath is not None:
            args.fitsfile[ii] = "%s/%s" % (args.opt_imagepath, args.fitsfile[ii])
        if not os.path.isfile(args.fitsfile[ii]):
            _logger.error("File %s does not exist. Giving up." % (args.fitsfile[ii]))
            sys.exit(0)

    args.fitsfile = Table(data=[np.asarray(args.fitsfile)], names=['FILENAME', ])

    return args


def main():
    args = parseCommandLine()
    database = noisegaindb(args.database) if args.database is not None else None

    if (database is not None) and args.noreprocessing:
        for inputname in args.fitsfile:
            if database.checkifalreadyused(os.path.basename(inputname)):
                _logger.info(
                    "File %s was already used in a noise measurement. Skipping this entire batch."
                    % inputname
                )
                exit(0)

    do_noisegain_for_fileset(args.fitsfile, database, args)

    database.close()


if __name__ == "__main__":
    main()
