import logging
import os
import astropy
import numpy as np
from astropy.io import fits
import argparse
import matplotlib.pyplot as plt
import lcocommissioning.common.SourceCatalogProvider  as cp
from astropy.stats import sigma_clipped_stats   



mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

# Example function to create an astigmatism map
def create_astigmatism_map(fitsobj):
    """
    Create a simple astigmatism map using matplotlib.

    Parameters:
        data (numpy.ndarray): 2D array representing the astigmatism data.
        title (str): Title of the plot.
    """
    # Extract data from the FITS file

    myCatalogProvider = cp.SEPSourceCatalogProvider()
    object_catalog, wcs  = myCatalogProvider.get_source_catalog_from_fitsobject(fitsobj, ext='SCI', minarea=20, deblend=0.5, noprocess=True, det_thres=20)
    imagename = os.path.basename (fitsobj.filename())

    imagedata = fitsobj['SCI'].data

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Display the data as an image
    mean, median, stddev = astropy.stats.sigma_clipped_stats(imagedata, sigma=3.0)
    cax = ax.imshow(imagedata, cmap='gist_gray', origin='lower',  vmin=median-2*stddev, vmax=median+5*stddev)

   
    u = object_catalog['ellipticity']
    v = object_catalog['ellipticity']
    u = u * np.cos(object_catalog['theta']*np.pi/180) 
    v = v * np.sin(object_catalog['theta']*np.pi/180) 
    ax.quiver(object_catalog['x'], object_catalog['y'],u,v, pivot='tail', angles='xy', color='cyan')

    # Set title and labels
    ax.set_title("Astigmatism Map")
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")


    plt.title (imagename)
    # Show the plot
    plt.savefig (f'astigmatism-{imagename}.png', dpi=150)
    # Close the FITS file


def parse_arguments():
    """
    Parse command-line arguments.
     Returns:
        argparse.Namespace: Parsed arguments containing the file path.
    """
    parser = argparse.ArgumentParser(description="Generate an astigmatism map from a FITS file.")
    parser.add_argument("file_path", type=str, nargs='+', help="Path to the FITS file containing astigmatism data.")

    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO'],
                        help='Set the debug level')
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args

# Example usage
if __name__ == "__main__":
    args = parse_arguments()
    
    for image in args.file_path:
        print ("opening fits file " + image)
        fitsobj = fits.open(image)
        create_astigmatism_map(fitsobj)
        fitsobj.close()