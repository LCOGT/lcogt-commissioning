import astropy.stats
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import argparse
import logging
import matplotlib.pyplot as plt
import os
import json 

log=logging.getLogger(__name__)
logging.getLogger('matplotlib.font_manager').disabled = True

def parse_arguments():
    parser = argparse.ArgumentParser(description='Annotate LCO guide frame.')
    parser.add_argument('filenames', type=str, nargs="+", help='Path to the FITS file')
    parser.add_argument('--loglevel', dest='log_level', default='DEBUG', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    return args

import xmltodict
import xml.etree.ElementTree as ET

def read_xml_file(xml_filename):
    tree = ET.parse(xml_filename)
    xml_data = tree.getroot()
    #here you can change the encoding type to be able to set it to the one you need
    xmlstr = ET.tostring(xml_data, encoding='utf-8', method='xml')
    tree =xmltodict.parse(xmlstr, process_namespaces=False)
    #print (json.dumps(tree, indent=3))
    xce = []
    yce = []
    if int(tree['ns0:processedImage']['numberOfValidCentroids']) > 1:
        for centroid in tree['ns0:processedImage']['centroids']:
            xce.append(float(centroid['pixel']['x']))
            yce.append(float(centroid['pixel']['y']))
    if int(tree['ns0:processedImage']['numberOfValidCentroids']) == 1:
            centroid = tree['ns0:processedImage']['centroids']
            xce.append(float(centroid['pixel']['x']))
            yce.append(float(centroid['pixel']['y']))
        
    xce = np.asarray(xce)
    yce = yce = np.asarray(yce)
    log.debug(f"Read {len(xce)} centroids from {xml_filename}")
    return xce, yce


def replace_raw_with_cat(file_path):

    retval =  file_path.replace("raw", "cat").replace( "fits", "fits.guide.xml").replace (".fz", "").replace ("g00", "g01").replace ("flash", "cat")
    log.debug(f"Replacing {file_path} with {retval}")
    return retval


def main():
    args = parse_arguments()
    # Load the FITS file
    for filename in args.filenames:
        log.info (filename)
        basename = os.path.basename(filename)
        catlogname = replace_raw_with_cat(filename)
        xce, yce = read_xml_file(catlogname)

        hdu_list = fits.open(filename)
        hdu = hdu_list[0]
        data = hdu.data

        # Plot the image with WCS projection
        plt.figure()
        ax = plt.subplot()
        mean, median, std = astropy.stats.sigma_clipped_stats(data, sigma=3.0)
        ax.imshow(data, origin='lower', cmap='gray', vmin=median-std, vmax=median+3* std)
        plt.scatter (xce, yce, s=25, marker="o", color='none', linewidths=1, edgecolors='red')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        plt.grid(color='white', ls='solid')
        plt.savefig (f"{basename}_soirces.png", dpi=200)



if __name__ == "__main__":
    main()
   

  