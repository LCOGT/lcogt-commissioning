from lcocommissioning.common.common import simpledateformat
from lcocommissioning.common.muscatfocusdb_orm import muscatfocusdb
import matplotlib.pyplot as plt
import logging
_logger = logging.getLogger(__name__)
logging.getLogger('matplotlib').setLevel(logging.FATAL)

_database = muscatfocusdb('sqlite:///muscatfocus.sqlite')


alldata =_database.getMeasurementsForMuscat('mc04')

gooddata = alldata[ alldata ['seeing_r'] < 2]

plt.figure()
plt.plot (gooddata['dateobs'], gooddata['temperature'],'.')
plt.xlabel("DATE-OBS")
plt.ylabel("Temp [deg C]")
simpledateformat()
plt.savefig ('muscatfocus_temperature.png')


plt.figure()
df_g = gooddata['focus_g'] - gooddata['focus_r']
df_i = gooddata['focus_i'] - gooddata['focus_r']
df_z = gooddata['focus_z'] - gooddata['focus_r']
plt.plot (gooddata['temperature'], df_g, '.', label="g")
plt.plot (gooddata['temperature'], df_i, '.', label="g")
plt.plot (gooddata['temperature'], df_z, '.', label="g")

plt.xlabel ("Focus temeprature [deg C]")
plt.ylabel ("Focus offset from rp channel")
plt.legend()
plt.savefig ("muscat_deltafocus_t.png")

plt.figure()
df_g = gooddata['focus_g'] - gooddata['focus_r']
df_i = gooddata['focus_i'] - gooddata['focus_r']
df_z = gooddata['focus_z'] - gooddata['focus_r']
plt.plot (gooddata['dateobs'], df_g, '.', label="g")
plt.plot (gooddata['dateobs'], df_i, '.', label="g")
plt.plot (gooddata['dateobs'], df_z, '.', label="g")

plt.xlabel ("Danie-Obs")
plt.ylabel ("Focus offset from rp channel")
plt.legend()
simpledateformat()

plt.savefig ("muscat_deltafocus_dateobs.png")

_database.close()