import astropy
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle




def readTable (file):
    data = ascii.read (file)
    return (data)
    print (data['X(FITS)_T1'])




d = readTable ("/Users/dharbeck/lco/deltarho/lco_data-20220805-93/Measurements.xls")
x = d['X(FITS)_T1']
x = x-np.median (x)
y = d['Y(FITS)_T1']
y = y-np.median(y)
t = d['JD_UTC']
t = (t-t[0] ) * 24*60*60

plt.plot (t,x,label='x')
plt.plot (t,y,label='y')
plt.ylim([-10,10])
plt.legend()
plt.xlabel ("t[s]")
plt.savefig ("t-xy.png")

plt.clf()
frequency_x, power_x = LombScargle(t, x).autopower()
frequency_y, power_y = LombScargle(t, y).autopower()
plt.plot (frequency_x,power_x, label='power x')
plt.plot (frequency_y,power_y, label='power y')
#plt.xlim([0,100])
plt.legend ()
plt.xlabel ('period[s]')
plt.savefig ('freq-xy.png')
