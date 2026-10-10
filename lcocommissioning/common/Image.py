import logging
import re

import numpy as np
from astropy.io import fits

from lcocommissioning.common.common import quadrantboundaries

_log = logging.getLogger(__name__)


class Image(object):
    """ Generic class to read in all SCI extensions from a fits file, be it fz compressed or not.

    Code is taken from LCO Banzai pipeline.

    """

    filename = None

    def __init__(self, filename, alreadyopenedhdu=False, overscancorrect=False, gaincorrect=False, skycorrect=False,
                 trim=True, minx=None, maxx=None, miny=None, maxy=None):
        """
        Load an image from a FITS file

        copypaste from banzai

        Parameters
        ----------
        filename: str
              Full path of the file to open

        Returns
        -------
        data: numpy array
            image data; will have 3 dimensions if the file was either multi-extension or
             a datacube
        header: astropy.io.fits.Header
              Header from the primary extension
        bpm: numpy array
            Array of bad pixel mask values if the BPM extension exists. None otherwise.
            extension_headers: list of astropy.io.fits.Header
                           List of headers from other SCI extensions that are not the
                           primary extension

        Notes
        -----
        The file can be either compressed or not. If there are multiple extensions,
        e.g. Sinistros, the extensions should be (SCI, 1), (SCI, 2), ...
        Sinsitro frames that were taken as datacubes will be munged later so that the
        output images are consistent
        """

        self.skycorrect = skycorrect
        self.gaincorrect = gaincorrect
        self.overscancorrect = overscancorrect

        if alreadyopenedhdu:
            hdulist = filename
        else:
            hdulist = fits.open(filename, 'readonly')

        # Get the main header
        self.primaryheader = hdulist[0].header
        if (alreadyopenedhdu or filename.endswith(".fz")) and (len(hdulist) == 2):
            for card in hdulist[1].header:
                if len(card) > 0:
                    self.primaryheader.append((card, hdulist[1].header[card]))
                    pass
        # Check for multi-extension fits
        self.extension_headers = []
        self.ccdsec = []
        self.ccdsum = []
        self.extver = []

        sci_extensions = self.get_extensions_by_name(hdulist, ['SCI', 'SPECTRUM', 'COMPRESSED_IMAGE'])
        _log.debug(f"SCI extensions found: {sci_extensions}")
        if len(sci_extensions) == 0:
            _log.warning("No SCI extenstion found in image %s. Forcing primary ." % filename)
            sci_extensions = [hdulist[0]]

        # Find out whow big dat are. Warning: assumption is that all extensions have same dimensions.
        datasec = sci_extensions[0].header.get('DATASEC')
        _log.debug("DATASEC: {}".format(datasec))

        if ((datasec is None) or not trim):
            _log.debug("Not trimming")
            cs = [1, sci_extensions[0].header['NAXIS1'], 1, sci_extensions[0].header['NAXIS2']]
        else:
            cs = self.fitssection_to_slice(datasec)
        _log.debug(f"Parsed datasec is: {cs}")

        if len(sci_extensions) > 0:
            # Generate internal storage array for pre-processed data
            self.data = np.zeros((len(sci_extensions), cs[3] - cs[2] + 1, cs[1] - cs[0] + 1), dtype=np.float32)

            for i, hdu in enumerate(sci_extensions):

                gain = 1.
                overscan = 0.
                extver = hdu.header.get('EXTVER', str(i + 1))
                self.extver.append(extver)

                self.ccdsec.append(cs)

                if overscancorrect:
                    overscan = self.get_overscan_from_hdu(hdu)

                if gaincorrect:
                    gain = float(hdu.header.get('GAIN', "1.0"))
                    hdu.header['GAIN'] = "1.0"

                hdu.header['OVLEVEL'] = overscan

                self.data[i, :, :] = (hdu.data[cs[2] - 1:cs[3], cs[0] - 1:cs[1]] - overscan) * gain

                skylevel = 0
                if skycorrect:
                    imagepixels = self.data[
                                  i, 20:-20, 20:-20]
                    skylevel = np.median(imagepixels)
                    std = np.std(imagepixels - skylevel)
                    skylevel = np.median(imagepixels[np.abs(imagepixels - skylevel) < 3 * std])
                    hdu.header['SKYLEVEL'] = skylevel
                    self.data[i, :, :] = self.data[i, :, :] - skylevel

                _log.debug(
                    "Correcting image extension #%d with gain / overscan / sky: % 5.3f % 8.1f  % 8.2f" % (
                        i, gain, overscan, skylevel))

                self.extension_headers.append(hdu.header)

        else:

            self.data = hdulist[0].data.astype(np.float32)

        try:
            self.bpm = hdulist['BPM'].data.astype(np.uint8)
        except KeyError:
            self.bpm = None

        if not alreadyopenedhdu:
            hdulist.close()

    def getccddata(self, extension, simulatext=True):
        """
        Returns the data as define by DATASEC in ehader
        :param extension:
        :return:
        """
        retdata = None
        if self.data.shape[0] > 1:
            retdata = self.data[extension]
        elif simulatext:
            extensions = quadrantboundaries(self.data[0])
            
            label, boundaries = extensions[extension]
            _log.debug (f"Simulating MEF extension {extension} on flat image {self.data.shape}: {label}, {boundaries}")
            y0, y1, x0, x1 = boundaries
            retdata = self.data[0,y0:y1, x0:x1]

            if self.skycorrect:
                imagepixels = retdata[20:-20, 20:-20]
                skylevel = np.median(imagepixels)
                std = np.std(imagepixels - skylevel)
                skylevel = np.median(imagepixels[np.abs(imagepixels - skylevel) < 3 * std])
                retdata = retdata - skylevel

        else:
            _log.warning("Image has only one extension, returning full image data")
            retdata = self.data
        return retdata

    def get_overscan_from_hdu(self, hdu, sig_rej=2, biassecheader='BIASSEC'):
        """ Calculate the overscan level of an image extension.
            Calculation is based on slice defined by header keyword.
        """
        overkeyword = hdu.header.get(biassecheader)
        if overkeyword is None:
            return 0

        if  hdu.header.get(biassecheader) in ('UNKNOWN', 'N/A'):
            _log.debug(f"Bias Section is undefined for extension {hdu}")
            return 0

        biassecslice = self.fitssection_to_slice(overkeyword)
        ovpixels = hdu.data[
                   biassecslice[2] + 1:biassecslice[3] - 1, biassecslice[0] + 1: biassecslice[1] - 1]
        overscanlevel = np.median(ovpixels)
        std = np.std(ovpixels)
        overscanlevel = np.mean(ovpixels[np.abs(ovpixels - overscanlevel) < sig_rej * std])
        return overscanlevel

    @staticmethod
    def fitssection_to_slice(keyword):
        integers = [int(n) for n in re.split(',|:', keyword[1:-1])]
        return integers

    def get_extensions_by_name(self, fits_hdulist, name):
        """
        Get a list of the science extensions from a multi-extension fits file (HDU list)

        Parameters
        ----------
        fits_hdulist: HDUList
                  input fits HDUList to search for SCI extensions

        name: str
          Extension name to collect, e.g. SCI

        Returns
        -------
        HDUList: an HDUList object with only the SCI extensions
        """
        # The following of using False is just an awful convention and will probably be
        # deprecated at some point
        extension_info = fits_hdulist.info(False)
        return fits.HDUList([fits_hdulist[ext[0]] for ext in extension_info if
                             ((ext[1] in name) and (fits_hdulist[ext[0]].data is not None))])


class ImageRegionReader(object):
    """ Read rectangular regions of the trimmed and overscan corrected science extensions of a fits file.

    reader.region(ext, y0, y1, x0, x1) returns the same pixels as Image(...).data[ext, y0:y1, x0:x1], but only the
    requested region is read via the HDU's .section, i.e., for compressed images only the overlapping tiles are
    decompressed, and the full image data are never loaded into (and cached by) the HDU. Use this instead of Image
    when only parts of large images are needed.
    """

    def __init__(self, hdulist, overscancorrect=False, trim=True):
        self.overscancorrect = overscancorrect

        # Same header merging as Image, but on a copy so the input hdulist is not modified.
        self.primaryheader = hdulist[0].header.copy()
        if len(hdulist) == 2:
            for card in hdulist[1].header:
                if len(card) > 0:
                    self.primaryheader.append((card, hdulist[1].header[card]))

        # Unlike Image.get_extensions_by_name, select by header only; accessing hdu.data would load the full image.
        self.sci_extensions = [hdu for hdu in hdulist if hdu.name in ('SCI', 'SPECTRUM', 'COMPRESSED_IMAGE')
                               and hdu.header.get('NAXIS', 0) > 0]
        if len(self.sci_extensions) == 0:
            _log.warning("No SCI extenstion found in image. Forcing primary .")
            self.sci_extensions = [hdulist[0]]

        # Assumption is that all extensions have same dimensions.
        header = self.sci_extensions[0].header
        datasec = header.get('DATASEC')
        if (datasec is None) or not trim:
            self.cs = [1, header['NAXIS1'], 1, header['NAXIS2']]
        else:
            self.cs = Image.fitssection_to_slice(datasec)
        self.shape = (len(self.sci_extensions), self.cs[3] - self.cs[2] + 1, self.cs[1] - self.cs[0] + 1)
        self._overscan = {}

    def overscan(self, extension):
        if not self.overscancorrect:
            return 0.
        if extension not in self._overscan:
            self._overscan[extension] = self._get_overscan_from_hdu(self.sci_extensions[extension])
        return self._overscan[extension]

    def region(self, extension, y0, y1, x0, x1):
        """ Return trimmed, overscan corrected pixels [y0:y1, x0:x1] (in trimmed coordinates) of an extension."""
        y0, y1 = max(0, y0), min(self.shape[1], y1)
        x0, x1 = max(0, x0), min(self.shape[2], x1)
        oy, ox = self.cs[2] - 1, self.cs[0] - 1
        pixels = self.sci_extensions[extension].section[oy + y0:oy + y1, ox + x0:ox + x1]
        return (pixels - self.overscan(extension)).astype(np.float32)

    @staticmethod
    def _get_overscan_from_hdu(hdu, sig_rej=2, biassecheader='BIASSEC'):
        """ Same as Image.get_overscan_from_hdu, but reads only the overscan pixels."""
        overkeyword = hdu.header.get(biassecheader)
        if overkeyword is None or overkeyword in ('UNKNOWN', 'N/A'):
            return 0
        biassecslice = Image.fitssection_to_slice(overkeyword)
        ovpixels = hdu.section[biassecslice[2] + 1:biassecslice[3] - 1, biassecslice[0] + 1: biassecslice[1] - 1]
        overscanlevel = np.median(ovpixels)
        std = np.std(ovpixels)
        overscanlevel = np.mean(ovpixels[np.abs(ovpixels - overscanlevel) < sig_rej * std])
        return overscanlevel
