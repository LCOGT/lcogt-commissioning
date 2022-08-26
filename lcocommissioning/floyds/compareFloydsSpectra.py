import logging
import sys
import matplotlib.pyplot as plt

_logger = logging.getLogger(__name__)
from astropy.io import fits
import numpy as np

plt.style.use('seaborn-darkgrid')


def plotSpectrum(fitsobject):
    xmin = float(fitsobject[0].header['XMIN'])
    xmax = float(fitsobject[0].header['XMAX'])
    npix = int(fitsobject[0].header['NAXIS1'])
    obj = fitsobject[0].header['OBJECT']
    dateobs = fitsobject[0].header['DATE-OBS']
    print(xmin, xmax, npix)
    x = np.arange(0, npix, 1)
    x = x * (xmax - xmin) / npix + xmin
    factor = 1
    if npix > 2000:
        factor = 0.27

    plt.plot(x, factor * fitsobject[0].data[1, 0, :], label=f"{obj} {dateobs}")


def main():
    filelist = sys.argv[1:]
    print(filelist)

    plt.figure(figsize=(20, 8), dpi=300)

    for file in filelist:
        spectrumobject = fits.open(file)
        plotSpectrum(spectrumobject)

    plt.ylabel("10**-22 ergs/s/cm2/A")
    plt.xlabel("lambda [A]")
    plt.title("Floyds Spectrum")
    plt.legend()

    plt.savefig("spectrum.png", bbox_inches="tight")


if __name__ == '__main__':
    main()
