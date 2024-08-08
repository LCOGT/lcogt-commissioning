import sys

import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
from  scipy import stats

''' Derive per pixel staistics: 
   output image of mean value per pixel
   output image of std deviation per pixel
   output histrogram of noise per pixel'''
plt.style.use("ggplot")



def create_memmap (fitsfile):
    fimage = fits.open(fitsfile)
    data = np.copy(np.asarray(fimage['SCI'].data[-4096:,-4096:]))
    data = data #- np.median (data)
    del fimage[0].data
    fimage.close()
    return data



def main():
    inputfiles = sys.argv[1:]
    print ("Reading in all the images")
    imagedata =  np.asarray([create_memmap(inputfile) for inputfile in inputfiles])
    print (f"input data shape: {imagedata.shape}")


    # look at temporal values: is the illumination constant?

    perimagemedian = [np.mean (image) for image in imagedata]
    print (perimagemedian)
    plt.figure()
    plt.plot (perimagemedian, "o")
    plt.title ("Image mean value in series")
    plt.savefig("timeline.png", dpi=150)
    plt.close()

    offset = np.mean (perimagemedian)
    for ii in range (len (inputfiles)):
        imagedata[ii] = imagedata[ii] - perimagemedian[ii] + offset

    # look at the mean image
    print ("Making a median stacked image")
    mean_image = np.mean (imagedata, axis=0)
    print (mean_image.shape)
    plt.figure()
    median = np.median (mean_image)
    std = stats.median_abs_deviation(mean_image, axis=None)
    plt.imshow (mean_image, clim=[median - 2 * std, median + 2*std])
    plt.colorbar()
    plt.title ("Mean stacked image")
    plt.savefig("mean_image.png", dpi=150)
    del mean_image
    plt.close()

    print ("Making an image of of the noise distribution")
    stdimage = np.std (imagedata, axis=0)

    plt.figure()
    median = np.median (stdimage)
    std = stats.median_abs_deviation(stdimage, axis=None)
    plt.imshow (stdimage, clim=[median - 3*std, median+3*std])
    plt.colorbar()
    plt.title ("By pixel standard deviation")
    plt.savefig("std_image.png", dpi=150)
    plt.close()


    stdimage_flattened = stdimage.flatten();
    idx_1d = stdimage_flattened.argsort()[-500:]

    # convert the idx_1d back into indices arrays for each dimension
    x_idx, y_idx = np.unravel_index(idx_1d, stdimage.shape)

    # Check that we got the largest values.
    for x, y, in zip(x_idx, y_idx):
        print(x,y,stdimage[x][y])
        plt.figure()
        minsinglepixel = imagedata[:,x,y]
        noise = stdimage[x][y]
        _ = plt.hist(minsinglepixel.flatten(),  bins=50, density = False, label=f"noise {noise} @ {x}/{y}")
        plt.title ("Distribution of Single pixel values")
        plt.xlabel("Per pixel value [ADU]")
        plt.ylabel("Density")
        #plt.xlim([350,610])
        plt.legend()
        plt.savefig(f"temp/pixelhistogram_{noise:06.3f}_{x}_{y}.png", dpi=150, bbox_inches="tight")
        np.savetxt (f"temp/pixeldist_{noise:06.3f}_{x}_{y}.txt" , minsinglepixel.flatten())
        plt.close()



    # look at the per pxiel distribution for a dedicated pixel

    maxy, maxx = np.unravel_index(np.argmax (stdimage), stdimage.shape)
    miny, minx = np.unravel_index(np.argmin (stdimage), stdimage.shape)

    print ("Looking at a single pixel distribution" )
    plt.figure()
    minsinglepixel = imagedata[:,miny,minx]
    _ = plt.hist(minsinglepixel.flatten(),  bins=50, density = True, label=f"minimum noise {np.min(stdimage)} @ {minx}/{miny}")
    np.savetxt ("mindistr.txt",  minsinglepixel.flatten())

    maxsinglepixel = imagedata[:,maxy,maxx]
    _ = plt.hist(maxsinglepixel.flatten(),  bins=50, density = True, label=f"maximum noise {np.max(stdimage)} @ {maxx}/{maxy}")
    np.savetxt ("maxdistr.txt" , maxsinglepixel.flatten())

    dualpixl = imagedata[:,155,229]
    np.savetxt ("bimodal.txt" , dualpixl.flatten())

    maxsinglepixel = imagedata[:,maxy+1,maxx]
    np.savetxt ("maxdistrp1.txt" , maxsinglepixel)
    maxsinglepixel = imagedata[:,maxy-1,maxx]
    np.savetxt ("maxdistrm1.txt" , maxsinglepixel)

    plt.title ("Distribution of Single pixel values")
    plt.xlabel("Per pixel value [ADU]")
    plt.ylabel("Density")
    plt.legend()
    plt.savefig("pixelhistogram.png", dpi=150, bbox_inches="tight")
    plt.close()

    print ("Looking at noise histogram distribution in image.")
    plt.figure()
    flattenedimage = stdimage.flatten()
    typicalnoise = stats.mode (flattenedimage)[0][0]

    _ = plt.hist(flattenedimage, density = True, bins=80, range=[0,80])
    plt.title ("Distribution of per pixel noise")
    plt.xlabel("Per pixel rms noise [ADU]")
    plt.ylabel("Density")
    plt.vlines(typicalnoise,0,np.max(_[0]), color='yellow', label=f"Mode: {typicalnoise:.2f} ADU")
    num = np.sum (flattenedimage[flattenedimage > 3*typicalnoise])
    all =  len (flattenedimage)
    print (f"Pixels >  {3*typicalnoise:.2f} ADU readnoise: {num} out of {all}. This is {num*100./all} %")
    plt.yscale('log')
    plt.legend()
    plt.savefig("noisehistogram.png", dpi=150, bbox_inches="tight")

    plt.close()


if __name__ == '__main__':
    main()
