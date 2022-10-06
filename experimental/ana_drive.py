''' Analyse a time series, generated with astroimagej,
    to diagnose issues with telescope axis drive oscillations. '''
import os
from astropy.io import ascii as pyascii
from astropy.timeseries import LombScargle
import numpy as np
import matplotlib.pyplot as plt


def read_table(file):
    ''' Read in a astroimageJ data file'''
    data = pyascii.read(file)
    return data


d = read_table(os.sys.argv[1])
x = d['X(FITS)_C1']
x = x - np.median(x)
y = d['Y(FITS)_C1']
y = y - np.median(y)
t = d['JD_UTC']
t = (t - t[0]) * 24 * 60 * 60

plt.plot(t, x, label='x')
plt.plot(t, y, label='y')
plt.ylim([-10, 10])
plt.legend()
plt.title (f"Sampling {(t[-1] / len(t)):5.2f}s")
plt.xlabel("t[s] since start")
plt.ylabel ("Relative drift from median position [pixel]")but the exi
plt.savefig("t-xy.png")

plt.clf()
frequency_x, power_x = LombScargle(t, x).autopower()
frequency_y, power_y = LombScargle(t, y).autopower()
plt.plot(1 / frequency_x, power_x, label='power x')
plt.plot(1 / frequency_y, power_y, label='power y')
plt.xlim([3, 300])
plt.yscale('log')
plt.legend()
plt.xlabel('period[s]')
plt.savefig('freq-xy.png')
