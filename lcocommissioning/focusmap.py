
import json
import os

import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import optimize

from lcocommissioning.common.SourceCatalogProvider import SEPSourceCatalogProvider

''' make  1D and 2D plots of focus / ellipticity / orientaion distriobitions

'''


def plotdistribution(imagecatalog, name=None):


    medianfwhm = np.median (imagecatalog['fwhm'])
    std = np.std (imagecatalog['fwhm'])
    medianfwhm = np.median (imagecatalog['fwhm'][ np.abs (imagecatalog['fwhm'] - medianfwhm) < std])
    print (f"Median fwhm is {medianfwhm} +/- {std}")
    good =  np.abs (imagecatalog['fwhm'] - medianfwhm) < std*2


    fig, ax = plt.subplots()
    #ax.set_box_aspect(1)
    plt.scatter (x=imagecatalog['x'][good], y=imagecatalog['y'][good], c=imagecatalog['ellipticity'][good], vmin=0,vmax=0.15)
    plt.colorbar()
    plt.xlabel ("x center")
    plt.ylabel ("y center")
    plt.savefig ('ellipticity2d.png')



    fig, ax = plt.subplots()
    ax.set_box_aspect(1)
    plt.scatter (x=imagecatalog['x'][good], y=imagecatalog['y'][good], c=imagecatalog['fwhm'][good], vmin=medianfwhm-std,vmax=medianfwhm+std)
    plt.colorbar()
    plt.title (f"FWHM {name}")

    plt.xlabel ("x center")
    plt.ylabel ("y center")
    plt.savefig ('fwhm2d.png')


    fig = plt.figure()

    plt.subplot (2,1,2)
    plt.plot (imagecatalog['x'][good], imagecatalog['ellipticity'][good],'.')
    plt.xlabel ("x position")
    plt.ylabel ("ellipticity")



    plt.ylim([0,0.2])
    plt.subplot (2,1,1)
    plt.title (f"Ellipticity {name}")
    plt.plot (imagecatalog['y'][good], imagecatalog['ellipticity'][good], '.')
    plt.xlabel ("y position")
    plt.ylabel ("ellipticity")
    plt.ylim([0,0.2])
    plt.savefig ('ellipticity1d.png')


    fig = plt.figure()

    plt.subplot (2,1,2)
    plt.plot (imagecatalog['x'][good], imagecatalog['fwhm'][good],'.')
    plt.xlabel ("x position")
    plt.ylabel ("fwhm [pixel]")



    plt.ylim([medianfwhm-2*std,medianfwhm+2*std])
    plt.subplot (2,1,1)
    plt.plot (imagecatalog['y'][good], imagecatalog['fwhm'][good], '.')
    plt.xlabel ("y position")
    plt.ylabel ("fwhm [pixel]")
    plt.ylim([medianfwhm-std,medianfwhm+std])
    plt.title (f"FWHM {name}")
    plt.savefig ('fwhm.png')


def main():
    catalogmaker = SEPSourceCatalogProvider(refineWCSViaLCO=False)

    for image in sys.argv[1:]:

        catalog, wcs = catalogmaker.get_source_catalog(image, )


        plotdistribution(catalog, name=os.path.basename(image))



if __name__ == '__main__':
    main()
