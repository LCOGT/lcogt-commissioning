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
                   'tlv': (0, 0), }

# Dictionary of NRES instances
nres_instruments = {'lsc': 'nres01',
                    'elp': 'nres02',
                    'cpt': 'nres03',
                    'tlv': 'nres04',
                    }

lco_1meter_sites = ['lsc', 'cpt', 'coj', 'elp', 'bpl']
lco_2meter_sites = ['ogg', 'coj']
lco_MuSCAT_sites = ['ogg']
lco_nres_sites = nres_instruments.keys()

lco_sinistro1m_cameras = ['fa02', 'fa03', 'fa04', 'fa05', 'fa06', 'fa07', 'fa08', 'fa11', 'fa12', 'fa14', 'fa15',
                          'fa16', 'fa19', ]

archon_readout_modes = ["full_frame", "central_2k_2x2"]

lco_muscat_instruments = ['mc03']
lco_muscat_readmodes = ['MUSCAT_FAST','MUSCAT_SLOW']

goodXTalkTargets = ['auto', '91 Aqr', 'HD30562', '15 Sex', '30Psc', '51Hya', 'Zet Boo']

default_constraints = {"max_airmass": 3,
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


def get_auto_target(targetlist, site, starttime, moonseparation=30, minalt=35):
    """ Go through a list of Simbad-resolvable objects and return the first visible object at the given site and time.

    :param targetlist: List of possible target names, as strings. Must resolve via simbad
    :param site:  LSC three letter site code
    :param starttime: datetime.datetime object for the time to check
    :param moonseparation:  minimum moon separation in degrees to consider object viable.
    :param minalt:  minimum altitude in degrees to consider an object viable.
    :return: Name of viable target or None
    """

    site = get_ephem_obj_for_site(site, starttime + dt.timedelta(minutes=30))
    moon = ephem.Moon()
    moon.compute(site)
    _log.debug("Finding suitable star for site %s. Moon phase is  %i %%" % (site, moon.moon_phase * 100))

    for starcandidate in targetlist:
        if 'auto' in starcandidate:
            continue
        radec = SkyCoord.from_name(starcandidate)
        s = ephem.FixedBody()
        s._ra = radec.ra.degree * math.pi / 180
        s._dec = radec.dec.degree * math.pi / 180
        s.compute(site)

        separation = (ephem.separation((moon.ra, moon.dec), (s.ra, s.dec)))

        alt = s.alt * 180 / math.pi
        separation = separation * 180 / math.pi

        altok = alt > minalt
        sepok = separation > moonseparation

        if (altok and sepok):
            _log.debug("\nViable star found: %s altitude % 4f moon separation % 4f" % (starcandidate, alt, separation))
            return starcandidate
        else:
            _log.debug("rejecting star %s - altitude ok: %s     moon separation ok: %s" % (starcandidate, altok, sepok))

    _log.debug("No viable star was found! full moon? returning None!")
    return None


def send_request_to_portal(requestgroup, dosubmit=False):
    _log.debug(json.dumps(requestgroup, indent=4))

    response = requests.post(
        'https://observe.lco.global/api/requestgroups/validate/',
        headers={'Authorization': 'Token {}'.format(VALHALLA_API_TOKEN)},
        json=requestgroup  # Make sure you use json!
    )
    print('API call return: {}'.format(response.content))

    if not dosubmit:
        print("Not submitting as per user request")
        return

    response = requests.post(
        'https://observe.lco.global/api/requestgroups/',
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
    print('View the observing request: https://observe.lco.global/requestgroups/{}/'.format(requestgroup_dict['id']))


def send_to_scheduler(user_request, dosubmit=False):
    """Submit a user request to LCO Scheduler via Valhalla interface
    """
    _log.fatal ("Not supported any more")
    exit (1)


def submit_observation(observation, dosubmit=False):
    """ Submit to LCO via an API call, direct submission"""

    if dosubmit:
        headers = {'Authorization': 'Token {token}'.format(token=VALHALLA_API_TOKEN)}
        response = requests.post(VALHALLA_URL + '/api/schedule/', json=observation, headers=headers)
        try:
            response.raise_for_status()
            _log.info(
                'Direct submitted observation with id: {0}. Check it at {1}/observations/{0}'.format(
                    response.json()['id'], VALHALLA_URL
                )
            )
        except Exception:
            _log.error(
                'Failed to submit observation: error code {}: {}'.format(response.status_code, response.content))
    else:
        _log.info("Not submitting block since CONFIRM  is not set")


import matplotlib.pyplot as plt
import matplotlib.dates as mdates


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
