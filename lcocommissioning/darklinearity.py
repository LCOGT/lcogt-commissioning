import logging
import sys
import os
import os.path
import numpy as np
import argparse
import astropy.time as astt
from astropy.table import QTable, Table, Column
from astropy.io import ascii
import math
from astropy.stats import sigma_clipped_stats


from lcocommissioning.common import lco_archive_utilities
from lcocommissioning.common.ccd_noisegain import dosingleLevelGain
from lcocommissioning.common.common import dateformat
from lcocommissioning.common.noisegaindb_orm import NoiseGainMeasurement, noisegaindb
from lcocommissioning.common.Image import Image
import matplotlib.pyplot as plt
from astropy.io import fits

_logger = logging.getLogger(__name__)
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

CONFMODE = None

def findkeywordinhdul(hdulist, keyword):
    # f&*& fpack!
    for ext in hdulist:
        val = ext.header.get(keyword)
        if val is not None:
            return val

    return None

def sort_by_darktime (filelist):
    global CONFMODE
    darkdict = {}
    for filename in filelist:
        fitsobj = fits.open (filename)
        exptime = findkeywordinhdul(fitsobj, 'EXPTIME')
        confmode = findkeywordinhdul(fitsobj, 'CONFMODE')
        if CONFMODE is None:
            CONFMODE = confmode
        assert (CONFMODE == confmode), f"CONFMODE mismatch {CONFMODE} {confmode}"
        fitsobj.close()
        if exptime not in darkdict:
            darkdict[exptime] = [filename]
        else:
            darkdict[exptime].append (filename)
    return darkdict

def parseargs():
    parser = argparse.ArgumentParser(
    description='General purpose single extension dark measurements',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fitsfiles', type=str, nargs='+',
                        help='Input fits files, must include at least two bias and two flat field frames.')

    args = parser.parse_args()
    return args

def combine_stack(list_of_files):
    assert (len(list_of_files)>0)
    print (f"Stacking {len(list_of_files)} images")
    targetbuffer = None
    for ii in range (len(list_of_files)):
        filename = list_of_files[ii]
        fitsobj=fits.open (filename)
        if targetbuffer is None:
            shape = fitsobj['SCI'].data.shape
            targetbuffer = np.zeros( (len(list_of_files),shape[0],shape[1]))
        targetbuffer[ii,::] = fitsobj['SCI'].data * 1.0
        fitsobj.close(fitsobj)
    average = np.average(targetbuffer,axis=0)
    return average


def main():
    args = parseargs()
    dark_dict = sort_by_darktime(args.fitsfiles)
    print (f"dark times set: {dark_dict.keys()}")
    assert (0 in dark_dict)

    maxexptime = np.max(list(dark_dict.keys()))
    print (f"Maximum Tdark={maxexptime}")
    stacked_images = {}
    
    levelstats = {}

    for exptime in dark_dict:
        stacked_images[exptime] = combine_stack(dark_dict[exptime])
        

    biasimage = stacked_images[0]
    stacked_images.pop(0,None)
    print ("Bias-correcting images")
    for exptime in stacked_images:
        stacked_images[exptime] = stacked_images[exptime] - biasimage
        mean = np.average(stacked_images[exptime][1000:2000,1000:2000])
        std = np.std(stacked_images[exptime][1000:2000,1000:2000])
        plt.imshow(stacked_images[exptime], clim=(mean-2*std,mean+3*std))
        plt.title (f"{CONFMODE}\nDark {exptime} sec\nDark Current {mean/exptime} ADU/sec")
        plt.colorbar()
        plt.savefig (f"dark_{CONFMODE}_{exptime}.png", bbox_inches="tight")
        plt.close()

    print ("Making ratio images")
   
    for exptime in stacked_images:
        if exptime == maxexptime:
            continue
        expected_ratio = maxexptime / exptime
        print (f" {maxexptime} / {exptime} = {expected_ratio}")
        
        correctiondelta = stacked_images[maxexptime] - expected_ratio *  stacked_images[exptime]
        image_mean = np.average(correctiondelta[1000:2000,3000:4000])
        std = np.std(correctiondelta[1000:2000,3000:4000])

        plt.imshow (correctiondelta, clim=(-8, 8))
        plt.title (f"{CONFMODE}\nDark residual Dark{maxexptime} - {expected_ratio} * {exptime} sec\nmean:{image_mean}")
        plt.colorbar()
        plt.savefig (f"dark_ratio_{CONFMODE}_{maxexptime}_{exptime}.png", bbox_inches="tight")
        plt.close()



if __name__ == '__main__':
    main()









