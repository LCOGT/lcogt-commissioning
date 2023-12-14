import argparse
import copy
import json
import logging
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
import datetime as dt
import lcocommissioning.common.common as common

_logger = logging.getLogger(__name__)



def createRequestsForStar_scheduler(context):

    absolutestart = context.start
    windowend = context.start + dt.timedelta(hours=context.schedule_window)

    location = {'telescope': "2m0a",
                'telescope_class': "2m0",
                'enclosure' : "clma",
                'site': context.site, }

    requestgroup = {"name": context.title,
                    "proposal": context.proposal,
                    "ipp_value": context.ipp,
                    "operator": "SINGLE" ,
                    "observation_type": "NORMAL",
                    "requests": []
                    }
    requests = {'configurations': [],
               'windows': [{"start": str(absolutestart), "end": str(windowend)}, ],
               'location': location}

    pm_ra = context.pm[0]
    pm_dec = context.pm[1]
    target = {
        "type": "ICRS",
        "name": "{} {}".format(context.title, context.targetname),
        "epoch": "2000.0000000",
        "equinox": "2000.0000000",
        "ra": "%10f" % context.radec.ra.degree,
        "dec": "%10f" % context.radec.dec.degree,
        "proper_motion_ra": pm_ra,
        "proper_motion_dec": pm_dec,
    }

    configuration =      {'type': 'SPECTRUM',
                          'instrument_type': '2M0-FLOYDS-SCICAM',
                          'target': target,

                          'acquisition_config': {
                              'mode': 'BRIGHTEST'
                          },
                          'guiding_config': {
                              'mode': 'ON',
                              'optional': False
                          },
                          'constraints': {},
                          'instrument_configs': [{
                              'exposure_time': context.exptime,
                              'exposure_count': int(context.expcnt),
                              'mode': 'default',
                              'rotator_mode': 'VFLOAT',
                              'optical_elements': {
                                  'slit': context.slit
                              },
                              'extra_params': {
                                  'defocus': min(5, context.defocus),
                              }
                          }]
                          }

    requests['configurations'].append (configuration)
    requestgroup['requests'].append (requests)



    _logger.debug(json.dumps(requestgroup, indent=4))
    common.send_request_to_portal(requestgroup, context.opt_confirmed)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Submit an engineering Floyds observation.')

    parser.add_argument('--targetname', default='auto', type=str,
                        help='Name of star for Floyds test observation. if none is given, or auto, Software will try to'
                             ' find a cataloged flux stadnard star. Name must resolve via Simbad; if resolve failes,'
                             ' program will exit.')
    parser.add_argument('--pm', type=float, nargs=2, default = None, help="proper motion RA DEC in marcsec / year")
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--site', choices=common.lco_2meter_sites, required=True, help="To which site to submit")

    parser.add_argument('--exp-cnt', type=int, dest="expcnt", default=1)
    parser.add_argument('--exptime', type=float, default=120)
    parser.add_argument('--slit', type=str, default="slit_1.2as", choices=['slit_1.2as', 'slit_2.0as', 'slit_6.0as'])
    parser.add_argument('--defocus', type=float, default=0.0, help="Amount to defocus star.")
    parser.add_argument('--start', default=None,
                        help="When to start Floyds observation. If not given, defaults to \"NOW\"")
    parser.add_argument('--user', default='daniel_harbeck', help="Which user name to use for submission")
    parser.add_argument('--schedule-window', default=2, type=float)
    parser.add_argument('--title', default="Floyds commissioning")
    parser.add_argument('--proposal', default='LCOEngineering')
    parser.add_argument('--ipp', type=float, default=1.0, help="ipp value")
    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, block will be submitted.')

    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()

    args.instrument = 'floyds01' if "ogg" in args.site else 'floyds02'
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    if args.start is None:
        args.start = dt.datetime.utcnow()
    else:
        try:
            args.start = dt.datetime.strptime(args.start, "%Y%m%d %H:%M")
        except ValueError:
            _logger.error("Invalid start time argument: ", args.start)
            exit(1)

    if ('auto' in args.targetname):
        args.targetname = common.get_auto_target(common.goodFloydsFluxStandards, args.site, args.start)
        if args.targetname is None:
            print ("No viable target found; giving up")
            exit(1)

    try:
        _logger.debug("Resolving target name")
        args.radec = SkyCoord.from_name(args.targetname)
    except:
        print("Resolving target name failed, giving up")
        exit(1)
    if args.pm is None:
        if args.targetname in common.listofpm:
            args.pm = common.listofpm[args.targetname]
        else:
            args.pm = [0,0]
    print(f"Resolved target {args.targetname} at corodinates {args.radec} deg w/ pm= {args.pm} marcsec/yr")
    return args


def main():
    args = parseCommandLine()
    createRequestsForStar_scheduler(args)



if __name__ == '__main__':
    main()
