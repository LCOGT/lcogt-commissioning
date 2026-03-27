import abc
import logging
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import sep
from astropy.io import fits
from astropy.io.fits import ImageHDU, CompImageHDU
from astropy.table import Table
from astropy.wcs import WCS

from lcocommissioning.gaiaastrometryservicetools import (
    astrometryServiceRefineWCSFromCatalog,
)

__author__ = "drharbeck@gmail.com"

log = logging.getLogger(__name__)
MIN_AREA = 9
THRESHOLD = 10.0


class SourceCatalogProvider(metaclass=abc.ABCMeta):
    """Interface to get an source catalog in pixel x/y coordinates out of a FITS image"""

    @abc.abstractmethod
    def get_source_catalog(self, imagename) -> (Table, WCS):
        """

        :param imagename:
        :return: (sourcecatalog, WCS)
        """
        pass


class e91SourceCatalogProvider(SourceCatalogProvider):
    """Read a source catalog and WCS from a LCO level 91 reduced image frame."""

    def get_source_catalog(self, imagename):
        e91image = fits.open(imagename)
        # ra = e91image['SCI'].header['CRVAL1']
        # dec = e91image['SCI'].header['CRVAL2']
        # log.debug("Image is at RA / Dec %f %f " % (ra, dec))
        try:
            sourceCatalog = e91image["CAT"].data
            sourceCatalog["x"] = sourceCatalog["xwin"]
            sourceCatalog["y"] = sourceCatalog["ywin"]
            log.debug("Source Catalog has %d entries" % len(sourceCatalog))

        except:
            log.warning("%s - No extension 'CAT' available, skipping." % (e91image))
            e91image.close()
            return (None, None)

        # instantiate the initial guess WCS from the image header
        image_wcs = WCS(e91image["SCI"].header)
        e91image.close()
        return (sourceCatalog, image_wcs)


class SEPSourceCatalogProvider(SourceCatalogProvider):
    """Submit fits image to LCO GAIA astrometry service for WCS refinenement and run SEP source finder on image data."""

    def __init__(self, refineWCSViaLCO=True):
        self.refineWCSViaLCO = refineWCSViaLCO

    def get_source_catalog(self, imagename, ext=1, minarea=20, deblend=0.5):
        fitsobject = fits.open(imagename)
        result = self.get_source_catalog_from_fitsobject(
            fitsobject, ext=ext, minarea=minarea, deblend=deblend
        )
        fitsobject.close()
        return result

    def get_source_catalog_from_fitsobject(
        self, fitsimage, ext=1, minarea=20, deblend=0.5, det_thres = 10, noprocess=False
    ):

        image_wcs = None

        for hdu in fitsimage:
            # search the first valid WCS. Search since we might have a fz compressed file which complicates stuff.
            try:
                original_wcs = WCS(hdu.header)
                continue
            except:
                log.warning("NO RA/DEC found yet, trying next extension")
                original_wcs = None

        log.debug("FITS file provided WCS is:\n{}".format(original_wcs))

        # Create a source catalog
        # TODO:    Better job of identifying the correct fits extension
        gain = fitsimage[ext].header.get("GAIN", 3.2)
        try:
            if "TRIMSEC" in fitsimage[ext].header:
                datasec = fitsimage[ext].header["TRIMSEC"]
            else:
                datasec = fitsimage[ext].header["DATASEC"]
            cs = [int(n) for n in re.split(",|:", datasec[1:-1])]
            bs = [
                int(n) for n in re.split(",|:", fitsimage[ext].header["BIASSEC"][1:-1])
            ]
            ovpixels = fitsimage[ext].data[bs[2] + 1 : bs[3] - 1, bs[0] + 1 : bs[1] - 1]
            overscan = np.median_nan(ovpixels)
            std = np.std_nan(ovpixels)
            overscan = np.mean_nan(ovpixels[np.abs(ovpixels - overscan) < 2 * std])
            if not noprocess:
                log.debug(
                    "Overscan correction: %f +/- %f" % (overscan, std)
                )
                image_data = (
                    fitsimage[ext].data[cs[2] - 1 : cs[3], cs[0] - 1 : cs[1]] - overscan
            )
        except:
            log.debug("No overscan specified")
            image_data = fitsimage[ext].data

        image_data = image_data[40:-40, :]  # Sinsitro: cut away the thinning artifacts.

        image_data = image_data.astype(float)
        error = ((np.abs(image_data) * gain) + 8**2.0) ** 0.5 / gain
        backGround = sep.Background(image_data, bw=32, bh=32, fw=3, fh=3)
        image_data = image_data - backGround

        # find sources
        
        objects, segmap = sep.extract(
            image_data,
            det_thres,
            err=backGround.rms(),
            deblend_cont=deblend,
            minarea=minarea,
            segmentation_map=True,
        )
        objects = Table(objects)
        # cleanup
        objects = objects[objects["flag"] == 0]
        objects = prune_nans_from_table(objects)
        # fwhm = 2.0 * (np.log(2) * (objects['a'] ** 2.0 + objects['b'] ** 2.0)) ** 0.5
        fwhm = np.sqrt((objects["x2"] + objects["y2"]) / 2) * 2.3548
        objects["theta"][objects["theta"] > (np.pi / 2.0)] -= np.pi
        objects["theta"][objects["theta"] < (-np.pi / 2.0)] += np.pi
        objects["theta"] = objects["theta"] * 180 / np.pi
        objects["ellipticity"] = 1.0 - (objects["b"] / objects["a"])
        objects = objects[fwhm > 1.0]

        flux_radii, flag = sep.flux_radius(
            image_data,
            objects["x"],
            objects["y"],
            6.0 * objects["a"],
            [0.25, 0.5, 0.75],
            normflux=objects["flux"],
            subpix=5,
        )
        sig = 2.0 / 2.35 * flux_radii[:, 1]
        xwin, ywin, flag = sep.winpos(image_data, objects["x"], objects["y"], sig)
        # python to FITS zero point convention. lower left pixel in image is 1/1, not 0/0
        sourcecatalog = Table(
            [
                xwin + 1,
                ywin + 1,
                objects["flux"],
                objects["theta"],
                objects["ellipticity"],
                fwhm,
            ],
            names=["x", "y", "flux", "theta", "ellipticity", "fwhm"],
        )

        log.debug("Sep found {} sources in image".format(len(sourcecatalog["x"])))

        # Lets refine the WCS solution.
        # TODO: Define condition when we want to refine the WCS
        submitImageInsteadofCatalog = False
        if self.refineWCSViaLCO:

            log.info("Sending raw source catalog to astrometry.net service")
            image_wcs = astrometryServiceRefineWCSFromCatalog(
                sourcecatalog, original_wcs
            )
            if image_wcs is None:
                image_wcs = original_wcs
        else:
            image_wcs = original_wcs

        return sourcecatalog, image_wcs


def prune_nans_from_table(table):
    nan_in_row = np.zeros(len(table), dtype=bool)
    for col in table.colnames:
        nan_in_row |= np.isnan(table[col])
    return table[~nan_in_row]


L1FWHM = "L1FWHM"
FOCDMD = "FOCDMD"


def getboundariesForSection(sectionid, dimX, dimY):
    """Get bounddaries for a 3x3 grid on 2d coordiantes"""
    ### Quadrants:
    ###   6 7 8
    ###   3 4 5
    ###   0 1 2

    bx = np.linspace(0, dimX - 1, 4)
    by = np.linspace(0, dimY - 1, 4)

    minx = bx[sectionid % 3]
    maxx = bx[sectionid % 3 + 1]
    miny = by[sectionid // 3]
    maxy = by[sectionid // 3 + 1]

    return minx, maxx, miny, maxy


def getImageFWHM(imagename, minarea=20, deblend=0.5, sections=False):
    """Measure the FWHM of an image, tuned to get a reasonable FWHM also for defocussed images."""

    hdul = fits.open(imagename, "readonly", ignore_missing_end=True)

    deltaFocus = None
    for ii in range(len(hdul)):
        if FOCDMD in hdul[ii].header:
            deltaFocus = hdul[ii].header[FOCDMD]
            continue
    catalog = SEPSourceCatalogProvider(refineWCSViaLCO=False)
    fwhmcat = np.asarray([])
    thetacat = np.asarray([])
    ellipticitycat = np.asarray([])

    for ii in range(len(hdul)):
        if isinstance(hdul[ii], ImageHDU) or isinstance(hdul[ii], CompImageHDU):
            cat, wcs = catalog.get_source_catalog(
                imagename, ext=ii, minarea=minarea, deblend=deblend
            )
            fwhms = cat["fwhm"]
            thetas = cat["theta"]
            ellipticities = cat["ellipticity"]
            fwhmcat = np.append(fwhmcat, fwhms)
            thetacat = np.append(thetacat, thetas)
            ellipticitycat = np.append(ellipticitycat, ellipticities)

    if not sections:
        medianfwhm = fit_fwhm(fwhmcat)
        mediantheta = np.median(thetacat)
        medianellipticity = np.median(ellipticitycat)
        log.debug(
            f"Median FWHM: {medianfwhm} Median Theta: {mediantheta} Median Ellipticity: {medianellipticity}"
        )

    else:
        medianfwhm = {}
        mediantheta = {}
        medianellipticity = {}
        for sectionID in range(0, 9):
            minx, maxx, miny, maxy = getboundariesForSection(
                sectionID, hdul[ii].data.shape[1], hdul[ii].data.shape[0]
            )
            log.debug(f"Section {sectionID} {minx} < X < {maxx} ;  {miny} < Y < {maxy}")
            selectarea = cat["x"] > minx
            selectarea = selectarea & (cat["x"] <= maxx)
            selectarea = selectarea & (cat["y"] > miny)
            selectarea = selectarea & (cat["y"] <= maxy)
            fwhms = cat["fwhm"][selectarea]
            thetas = cat["theta"][selectarea]
            ellipticities = cat["ellipticity"][selectarea]

            medianfwhm[sectionID] = fit_fwhm(fwhms)
            mediantheta[sectionID] = np.median(thetas)
            medianellipticity[sectionID] = np.median(ellipticities)
            log.debug(
                f"Median FWHM: {medianfwhm[sectionID]} Median Theta: {mediantheta[sectionID]} Median Ellipticity: {medianellipticity[sectionID]}"
            )

    hdul.close()
    return deltaFocus, medianfwhm, mediantheta, medianellipticity


def fit_fwhm(fwhmcat):
    good = fwhmcat > 2
    for iter in range(4):
        medianfwhm = np.median(fwhmcat[good])

        fwhmstd = np.std(fwhmcat[good])
        good = (abs(fwhmcat - medianfwhm) < 2 * fwhmstd) & (fwhmcat > 2)
        if np.sum(good) > 10:
            medianfwhm = np.median(fwhmcat[good])
        else:
            log.debug("Not enough sources left in array. aborting")
            continue
    return medianfwhm


import os


def doimagegrid(image):
    minareas = [9, 15, 20, 25, 30]
    deblends = [0.005, 0.05, 0.1, 0.5, 1.0]
    titlename = os.path.basename(image)
    print(titlename)

    X, Y = np.meshgrid(minareas, deblends)
    Z = X * 0 + Y * 0
    for i in range(len(minareas)):
        for j in range(len(deblends)):
            focus, fwhm = getImageFWHM(image, minarea=X[i, j], deblend=Y[i, j])
            Z[i, j] = fwhm
            print("{} {} {}".format(X[i, j], Y[i, j], Z[i, j]))
    plt.pcolormesh(X, Y, Z)

    plt.xlabel("MINAREA")
    plt.ylabel("DEBLEND")
    plt.title("Estimated FWHM [pixels]\n{}".format(titlename))
    plt.colorbar()
    # plt.show()
    plt.savefig("fwhm_deblendminarea_grid_{}.png".format(titlename))
    plt.close()


if __name__ == "__main__":
    # TODO: Make this a test code
    logging.basicConfig(
        level=getattr(logging, "INFO"),
        format="%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s",
    )
    sourcecatalogProvider = SEPSourceCatalogProvider(refineWCSViaLCO=False)

    for image in sys.argv[1:]:
        doimagegrid(image)
        # getImageFWHM(image)
