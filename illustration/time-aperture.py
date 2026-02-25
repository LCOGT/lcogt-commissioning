import math
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates


def dateformat(starttime=None, endtime=None):
    """ Utility to prettify a plot with dates.
    """

    plt.xlim([starttime, endtime])
    plt.gcf().autofmt_xdate()
    years = mdates.YearLocator(5)  # every year
    months = mdates.MonthLocator(bymonth=[4, 7, 10])  # every month
    yearsFmt = mdates.DateFormatter('%Y')
    monthformat = mdates.DateFormatter('%Y')
    plt.gca().xaxis.set_major_locator(years)
    plt.gca().xaxis.set_major_formatter(yearsFmt)
    #plt.gca().xaxis.set_minor_locator(mdates.YearLocator(1))
    #plt.gca().xaxis.set_minor_formatter(monthformat)
    #plt.setp(plt.gca().xaxis.get_minorticklabels(), rotation=45)
    plt.yscale('log')
    plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=45)
    plt.gca().grid(which='minor')




year =     [ 2004, 2005, 2010,2017,2020,2023,2027,2030]
year = [datetime(year=x, month=1, day=1) for x in year]
aperture = np.asarray([   11,    5,   3.5, 1 , 2, 0.35, 1,  0.15 ])
multi = np.asarray(   [    1.,   1,     1, 13, 2,   10, 13,   30])
label = ["Salt", "Palomar", "WIYN", "LCO-1m", "LCO-2m", 'LCO Delta Rho', 'LCO-1m', 'Cube Sats?']

plt.style.use('ggplot')


plt.figure()
y = np.log10 (aperture)
y=aperture
collecting = (math.pi * (y*multi)**2)
plt.plot (year,y,'o', label="single aperture")
plt.plot (year,collecting, 'o', label="system collecting area" )
for ii in range(len (year)):
    plt.text (year[ii], y[ii], label[ii], rotation = 20 )

plt.legend()
plt.xlabel("Year")
plt.ylabel ("Aperture [m]\nSystem Collecting Area [m^2]")

dateformat( datetime(2000, 1, 1), datetime(2035, 12, 31))
#plt.ylim([-1,2])
plt.savefig ("yearap.png", bbox_inches='tight')
