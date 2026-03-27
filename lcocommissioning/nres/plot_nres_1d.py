from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import sys


def main ():
    file_path = sys.argv[1]
    hdu = fits.open(file_path)
    data = hdu["SPECTRUM"].data
    wavelength = data['wavelength']

    order = 71

    fig, ax = plt.subplots()
    for i, w in enumerate(wavelength):
        if (data['fiber'][i] != 1) & (data['order'][i] == order):
            pixel_range = np.where(w > 0)
            ax.plot(w, data["normflux"][i],)
            print(f"|{data['order'][i]} | {min(w[pixel_range])} | {max(w[pixel_range])}|")
            ax.set_xlim(min(w[pixel_range]),max(w[pixel_range]))

    #ax.set_ylim(-0.5,2.0)
    ax.set_ylabel("normflux")
    ax.set_xlabel("wavelength")
    plt.show()
    hdu.close()


if __name__ == '__main__':
    main()