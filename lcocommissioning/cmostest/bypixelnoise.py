import sys

import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

''' Derive per pixel staistics: 
   output image of mean value per pixel
   output image of std deviation per pixel
   output histrogram of noise per pixel'''



def create_memmap (fitsfile):
    fimage = fits.open(fitsfile)
    data = np.copy(np.asarray(fimage[1].data[100:350,100:350]))
    data = data - np.mean (data)
    del fimage[0].data
    fimage.close()
    return data



def main():
    inputfiles = sys.argv[1:]
    print ("Reading in all the images")
    imagedata =  np.asarray([create_memmap(inputfile) for inputfile in inputfiles])
    print (imagedata.shape)


    # look at temporal values: is the illumination constant?

    perimagemedian = [np.mean (image) for image in imagedata]
    print (perimagemedian)
    plt.figure()
    plt.plot (perimagemedian, "o")
    plt.title ("Iamge mean value in series")
    plt.savefig("timeline.png", dpi=150)
    plt.close()





    # look at the mean image
    print ("Making a median stacked image")
    mean_image = np.mean (imagedata, axis=0)
    print (mean_image.shape)
    plt.figure()
    median = np.median (mean_image)
    std = scipy.stats.median_abs_deviation(mean_image, axis=None)
    plt.imshow (mean_image, clim=[median - 2 * std, median + 2*std])
    plt.colorbar()
    plt.title ("Mean stacked image")
    plt.savefig("mean_image.png", dpi=150)
    plt.close()

    print ("Making an image of of trhe noise distribution")
    stdimage = np.std (imagedata, axis=0)

    plt.figure()
    median = np.median (stdimage)
    std = scipy.stats.median_abs_deviation(stdimage, axis=None)
    plt.imshow (stdimage, clim=[median - 3*std, median+3*std])
    plt.colorbar()
    plt.title ("By pixel standard deviation")
    plt.savefig("std_image.png", dpi=150)
    plt.close()






    # look at the per pxiel distribution for a dedicated pixel

    maxx, maxy = np.unravel_index(np.argmax (stdimage), stdimage.shape)
    minx, miny = np.unravel_index(np.argmin (stdimage), stdimage.shape)

    print ("Looking at a single pixel distribution" )
    plt.figure()
    minsinglepixel = imagedata[:,minx,miny]
    _ = plt.hist(minsinglepixel.flatten(),  bins=50, density = True, label=f"minimum noise {np.min(stdimage)} @ {minx}/{miny}")

    np.savetxt ("mindistr.txt",  minsinglepixel)


    maxsinglepixel = imagedata[:,maxx,maxy]
    _ = plt.hist(maxsinglepixel.flatten(),  bins=50, density = True, label=f"maximum noise {np.max(stdimage)} @ {maxx}/{maxy}")
    np.savetxt ("maxdistr.txt" , maxsinglepixel)

    plt.title ("Distribution of Single pixel values")
    plt.xlabel("Per pixel value [ADU]")
    plt.ylabel("Density")
    plt.legend()
    plt.savefig("pixelhistogram.png", dpi=150)
    plt.close()


    print ("Looking at noise histogram distribution in image.")
    plt.figure()
    _ = plt.hist(stdimage.flatten(), density = True, bins=100, range=[median-7*std, median+ 7* std])
    plt.title ("Distribution of per pixel noise")
    plt.xlabel("Per pixel noise [ADU]")
    plt.ylabel("Density")
    plt.vlines(median,0,np.max(_[0]), color='red')
    plt.savefig("noisehistogram.png", dpi=150)
    plt.close()


if __name__ == '__main__':
    main()
