import re
import sys
import numpy as np
from astropy.table import Table
import astropy.time as astt
from datetime import datetime
import matplotlib.pyplot as plt


pid_to_name = {'01550819': 'Sophia',
               '04573119': 'ep02',
               '04574619': 'ep03',
               '04573819': 'ep04',
}

def parselogfile (fname):
    regexstring = '.*picam: (\d+) Expected.*times: (\d\d:\d\d:\d\d.\d\d\d) (\d\d:\d\d:\d\d.\d\d\d) (\d\d:\d\d:\d\d.\d\d\d) (\d\d:\d\d:\d\d.\d\d\d)'

    f = open (fname)
    lines = re.findall (regexstring, f.read())
    lines = np.asarray (lines).T
    f.close()

    header = ['cameraid','expstart','measstart','expend','measend']
    expstart = [ datetime.strptime(v,'%H:%M:%S.%f') for v in lines[1]  ]
    measstart = [ datetime.strptime(v,'%H:%M:%S.%f') for v in lines[2]  ]

    expend = [ datetime.strptime(v,'%H:%M:%S.%f') for v in lines[3]  ]
    measend = [ datetime.strptime(v,'%H:%M:%S.%f') for v in lines[4]  ]

    T=Table ([lines[0],expstart,measstart, expend, measend], names=header)
    return (T)


def analyselogfile(T):
    cameras = np.unique (T['cameraid'])
    timedelta = T['measstart'] - T['expstart']
    exptimemeas = T['measend'] - T['measstart']
    exptimereq = T['expend'] - T['expstart']
    starttime = T['measstart']
    # Exp start vs measer start
    plt.figure()
    for camera in cameras:
        ctimedelta = timedelta[camera==T['cameraid']]
        ctimedelta = [ v.total_seconds() for v in ctimedelta]
        avg = np.average(ctimedelta)
        std = np.std (ctimedelta)
        plt.hist(ctimedelta, bins=600,  histtype='step', label=f"{pid_to_name[camera]:7s} {avg: 4.2f} +/- {std: 4.2f}")
    plt.legend()
    plt.xlabel ('Measured - Predicted Start of exposure [s]')
    plt.ylabel ('frequency')
    plt.title ('Muscat3 DATE-OBS  differences')
    plt.savefig ('princeton_timestampbdelta.pdf')

    # measured exposure time vs expected exposure time
    plt.figure()
    for camera in cameras:
        exptime_m  = exptimemeas[camera==T['cameraid']]
        exptime_m = np.asarray([ v.total_seconds() for v in exptime_m])
        ctimedelta = timedelta[camera==T['cameraid']]
        ctimedelta = [ v.total_seconds() for v in ctimedelta]
        stime = starttime[camera==T['cameraid']]
        print (ctimedelta)
        avg = np.average(ctimedelta)
        std = np.std (ctimedelta)
        plt.plot (stime, ctimedelta, ',', label=f"{pid_to_name[camera]:7s} {avg: 4.2f} +/- {std: 4.2f}")
    plt.legend()
    plt.xlabel ('Start time (hhmm)')
    plt.ylabel ('Delta DATE-OBS')

    plt.savefig ('princeton_starttimedelta.pdf')





t = parselogfile (sys.argv[1])
analyselogfile(t)