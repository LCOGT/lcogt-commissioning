from lcocommissioning.common.common import simpledateformat
from lcocommissioning.common.muscatfocusdb_orm import muscatfocusdb
import matplotlib.pyplot as plt
import logging
import datetime as dt

_logger = logging.getLogger(__name__)
logging.getLogger('matplotlib').setLevel(logging.FATAL)

_database = muscatfocusdb('sqlite:///muscatfocus.sqlite')


for muscat in ['mc03','mc04']:

    bestafter = dt.datetime.fromisoformat("2023-10-18") if muscat == 'mc04' else dt.datetime.fromisoformat("2023-12-10")

    alldata =_database.getMeasurementsForMuscat(muscat)
    alldata = alldata[ alldata ['seeing_r'] < 1.5]
    gooddata = alldata [ alldata['dateobs'] > bestafter]


    df_g = gooddata['focus_g'] - gooddata['focus_r']
    df_i = gooddata['focus_i'] - gooddata['focus_r']
    df_z = gooddata['focus_z'] - gooddata['focus_r']
    df_ga = alldata['focus_g'] - alldata['focus_r']
    df_ia = alldata['focus_i'] - alldata['focus_r']
    df_za = alldata['focus_z'] - alldata['focus_r']

    plt.figure()
    plt.plot (alldata['dateobs'], alldata['temperature'], '.', color='blue')
    plt.xlabel("DATE-OBS")
    plt.ylabel("Temp [deg C]")
    plt.title(f"{muscat}")
    simpledateformat()
    plt.savefig (f'{muscat}_focus_temperature.png')

    plt.figure()
    plt.plot (gooddata['temperature'], df_g, '.', label="g")
    plt.plot (gooddata['temperature'], df_i, '.', label="i")
    plt.plot (gooddata['temperature'], df_z, '.', label="z")
    plt.xlabel ("Focus temperature [deg C]")
    plt.ylabel ("Focus offset from rp channel")
    plt.legend()
    plt.title(f"{muscat}")
    plt.savefig (f"{muscat}_deltafocus_t.png")

    plt.figure()
    plt.plot (alldata['dateobs'], df_ga, '.',  label="g")
    plt.plot (alldata['dateobs'], df_ia, '.',  label="i")
    plt.plot (alldata['dateobs'], df_za, '.', label="z")
    plt.xlabel ("Date-Obs")
    plt.ylabel ("Focus offset from rp channel")
    plt.legend()
    plt.title(f"{muscat}")
    simpledateformat()
    plt.savefig (f"{muscat}_deltafocus_dateobs.png")

    plt.figure()
    plt.plot (alldata['dateobs'], alldata['focus_r'], '.', label="r focus correction")
    plt.xlabel ("Date-Obs")
    plt.ylabel ("Focus correction in rp")
    plt.title(f"{muscat}")
    plt.legend()
    simpledateformat()
    plt.savefig (f"{muscat}_rpfocus_dateobs.png")



_database.close()