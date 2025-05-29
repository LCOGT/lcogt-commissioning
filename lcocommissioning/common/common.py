import datetime as dt
import json
import logging
import os

import ephem
import math
import requests
from astropy.coordinates import SkyCoord

_log = logging.getLogger(__name__)
# LCO Request submission definitions
VALHALLA_URL = os.getenv('VALHALLA_URL', 'http://internal-observation-portal.lco.gtn')
VALHALLA_API_TOKEN = os.getenv('VALHALLA_API_TOKEN', '')

# LCO sites
lco_site_lonlat = {'bpl': (-119.863103, 34.433161),
                   'coj': (149.0708466, -31.2728196),
                   'cpt': (20.8124, -32.3826),
                   'elp': (-104.015173, 30.679833),
                   'lsc': (-70.8049, -30.1673666667),
                   'ogg': (-156.2589, 34.433161),
                   'sqa': (-120.04222167, 34.691453333),
                   'tfn': (-16.511544, 28.300433),
                   'tlv': (30.595833, 30.595833),
                   }

lco_site_alt = {'ogg': 3055,
                'coj':  1116,
                'cpt': 1760,
                'lsc': 2198,
                'elp': 2070,
                'tlv': 875,
                'ngq': 5100,

                }

# Dictionary of NRES instances
nres_instruments = {'lsc': 'nres01',
                    'elp': 'nres02',
                    'cpt': 'nres03',
                    'tlv': 'nres04',
                    }

lco_1meter_sites = ['lsc', 'cpt', 'coj', 'elp', 'tfn', 'bpl']
lco_2meter_sites = ['ogg', 'coj']
lco_MuSCAT_sites = ['ogg']
lco_nres_sites = nres_instruments.keys()

lco_sinistro1m_cameras = ['fa02', 'fa03', 'fa04', 'fa05', 'fa06', 'fa07', 'fa08', 'fa11', 'fa12', 'fa14', 'fa15',
                          'fa16', 'fa19', ]

archon_readout_modes = ["full_frame", "central_2k_2x2"]

lco_muscat_instruments = ['mc03']
lco_muscat_readmodes = ['MUSCAT_FAST', 'MUSCAT_SLOW']

goodXTalkTargets = ['auto', '30Psc', 'HD30562', '15 Sex',   '51Hya', 'Zet Boo', '9Peg', '91 Aqr',   ]
goodNRESFluxTargets = ['auto', 'HR9087', 'HR1544', 'HR4468', 'HD93521', 'HR3454', 'HR5501', 'HR7596',  'HD60753' ]
goodFloydsFluxStandards = ['auto', 'HZ 43', 'GD 71', 'BD+284211', 'HZ 44', 'L745-46A', 'Feige 110', 'EGGR274', 'EG21']

# list of proper motions for some stars:
listofpm = {'HR9087'  : [18.844,    -9.700],
            'HR1544'  : [  1.41,   -29.91],
            'HR3454'  : [-19.39,    -1.08],
            'HR5501'  : [-40.419,   -8.096],
            'HR7596'  : [ 39.126,  -13.931],
            'HR4468'  : [-59.38,     2.55 ],
            'HD93521' : [  0.220,    1.717],
            'L745-46A':[1138.690, -542.556],
            'HZ 43': [ 	-157.969, -107.168],
            'GD 71':   [  76.728, -172.960 ],
            'BD+284211': [ 	-34.809, -56.937 ],
            'HZ 44':     [ 	-66.004, -4.408],
            'L745-46A' : [ 	1138.690, -542.556],
            'Feige 110': [ 	-8.533, -0.592 ],
            'EG21': [ 	39.668, -103.237],
            'EGGR274' : [ 	77.398, 0.386 ],
            'HD60753' : [-2.897,5.617 ],

          }

# have a list of cached coordiantes so we do not alwasys need to look it up
listofchachedcoordiantes = {
    'HR9087': [0.4560297510252, -03.0275060205485],
    'HR1544': [72.6530124271, +8.9001803703],
    'HR3454': [130.8061458029, +03.3986629753],
    'HR5501': [221.3758584061525, +0.7172718096401],
    'HR7596': [298.6866476452963, +0.2736259408895],
    'HD93521': [162.0979650581188, +37.5703027604956],
    'HR4468': [174.1704723071, -9.8022475661]
}

default_constraints = {"max_airmass": 2,
                       "min_lunar_distance": 30.0, }


def get_ephem_obj_for_site(sitecode, dateobs):
    site = ephem.Observer()
    lon, lat = lco_site_lonlat[sitecode]
    site.lat = lat * math.pi / 180
    site.lon = lon * math.pi / 180
    site.date = ephem.Date(dateobs)
    return site


def is_valid_lco_site(sitecode):
    return sitecode in lco_site_lonlat


def get_auto_target(targetlist, sitecode, starttime, moonseparation=29, minalt=45):
    """ Go through a list of Simbad-resolvable objects and return the first visible object at the given site and time.

    :param targetlist: List of possible target names, as strings. Must resolve via simbad
    :param site:  LSC three letter site code
    :param starttime: datetime.datetime object for the time to check
    :param moonseparation:  minimum moon separation in degrees to consider object viable.
    :param minalt:  minimum altitude in degrees to consider an object viable.
    :return: Name of viable target or None
    """

    site = get_ephem_obj_for_site(sitecode, starttime + dt.timedelta(minutes=1))
    moon = ephem.Moon()
    moon.compute(site)
    _log.debug(f"Finding suitable star for site {sitecode} at LST {site.sidereal_time()}. Moon phase is  {moon.moon_phase*100} %%" )

    for starcandidate in targetlist:
        if 'auto' in starcandidate:
            continue
        if starcandidate in listofchachedcoordiantes:
            cradec = listofchachedcoordiantes[starcandidate]
            radec = SkyCoord(cradec[0], cradec[1], unit='deg', frame='icrs', )
        else:
            radec = SkyCoord.from_name(starcandidate)

        s = ephem.FixedBody()
        s._ra = radec.ra.degree * math.pi / 180
        s._dec = radec.dec.degree * math.pi / 180
        s.compute(site)

        separation = math.fabs (ephem.separation((moon.ra, moon.dec), (s.ra, s.dec)))

        alt = s.alt * 180 / math.pi
        separation = separation * 180 / math.pi

        altok = alt >= minalt
        sepok = separation >= moonseparation

        if (altok and sepok):
            _log.debug("\nViable star found: %s altitude % 4f moon separation % 4f" % (starcandidate, alt, separation))
            return starcandidate
        else:
            _log.debug(f"rejecting star {starcandidate:10s}  {radec.ra:9.5f} {radec.dec:9.5f}- altitude ok: {altok!s:^5} {alt:5.1f}    moon separation ok: {sepok!s:^5} {separation: 6.1f} < {moonseparation: 4.1f}")

    _log.debug("No viable star was found! full moon? returning None!")
    return None


def send_request_to_portal(requestgroup, dosubmit=False, url="https://observe.lco.global"):
    _log.debug(json.dumps(requestgroup, indent=4))
    if VALHALLA_API_TOKEN=='':
        _log.error("Environment Variable VALHALLA_API_TOKEN is not set; please set with your credentials. ")
        exit(1)
    response = requests.post(
        f'{url}/api/requestgroups/validate/',
        headers={'Authorization': 'Token {}'.format(VALHALLA_API_TOKEN)},
        json=requestgroup  # Make sure you use json!
    )
    print('API call return: {}'.format(response.content))

    if not dosubmit:
        print("Not submitting as per user request")
        return

    response = requests.post(
        f'{url}/api/requestgroups/',
        headers={'Authorization': 'Token {}'.format(VALHALLA_API_TOKEN)},
        json=requestgroup  # Make sure you use json!
    )
    # Make sure the API call was successful
    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as exc:
        print('API call failed: {}'.format(response.content))
        return

    requestgroup_dict = response.json()  # The API will return the newly submitted requestgroup as json

    # Print out the url on the portal where we can view the submitted request
    print(f'View the observing request: {url}/requestgroups/{requestgroup_dict["id"]}/')


def send_to_scheduler(user_request, dosubmit=False):
    """Submit a user request to LCO Scheduler via Valhalla interface
    """
    _log.fatal("Not supported any more")
    exit(1)


def submit_request_group(observation, dosubmit=False):
    """ Submit to LCO via an API call, direct submission"""

    if dosubmit:
        if VALHALLA_API_TOKEN=='':
            _log.error("Environment Variable VALHALLA_API_TOKEN is not set; please set with your credentials. ")
            exit(1)
        headers = {'Authorization': 'Token {token}'.format(token=VALHALLA_API_TOKEN)}
        response = requests.post(VALHALLA_URL + '/api/schedule/', json=observation, headers=headers)
        try:
            response.raise_for_status()
            _log.info(
                'Submitted request group with id: {0}. Check it at {1}/observations/{0}'.format(
                    response.json()['id'], VALHALLA_URL
                )
            )
        except Exception:
            _log.error(
                'Failed to submit request group: error code {}: {}'.format(response.status_code, response.content))
    else:
        _log.info("Not submitting block since CONFIRM  is not set")


import matplotlib.pyplot as plt
import matplotlib.dates as mdates



def simpledateformat ():
    plt.gcf().autofmt_xdate()
    plt.setp(plt.gca().xaxis.get_minorticklabels(), rotation=45)
    plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=45)
    plt.gca().grid(which='minor')
    
def dateformat(starttime=None, endtime=None):
    """ Utility to prettify a plot with dates.
    """

    plt.xlim([starttime, endtime])
    plt.gcf().autofmt_xdate()
    years = mdates.YearLocator()  # every year
    months = mdates.MonthLocator(bymonth=[4, 7, 10])  # every month
    yearsFmt = mdates.DateFormatter('%Y %b')
    monthformat = mdates.DateFormatter('%b')
    plt.gca().xaxis.set_major_locator(years)
    plt.gca().xaxis.set_major_formatter(yearsFmt)
    plt.gca().xaxis.set_minor_locator(months)
    plt.gca().xaxis.set_minor_formatter(monthformat)
    plt.setp(plt.gca().xaxis.get_minorticklabels(), rotation=45)
    plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=45)
    plt.gca().grid(which='minor')
