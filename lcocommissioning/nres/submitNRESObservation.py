import argparse
import json
import logging

import astropy
from astropy.coordinates import SkyCoord
import datetime as dt
import lcocommissioning.common.common as common

_logger = logging.getLogger(__name__)

defaultconstraints = {"max_airmass": 2.5,
                      "min_lunar_distance": 45.0, }


def createNRESRequestsConfiguration(args):
    configuration = {
        'type': 'NRES_SPECTRUM',
        'instrument_type': '1M0-NRES-SCICAM' if not args.commissioning else '1M0-NRES-COMMISSIONING',
        'guiding_config': {'mode': 'ON', 'optional' : False},
        'acquisition_config': {'mode': 'WCS' if args.forcewcs else 'BRIGHTEST'},
        'instrument_configs': [{
            'exposure_time': args.exptime,
            'exposure_count': args.expcnt,
            'mode': 'default'
        },]
    }
    return configuration


def createRequest(args):
    requestgroup = {"name": f'NRES engineering {args.targetname}',
                    "proposal": args.proposalid,
                    "ipp_value": args.ipp,
                    "operator": "SINGLE",  # "MANY" if args.dither else "SINGLE",
                    "observation_type": "NORMAL",
                    "requests": []
                    }

    absolutestart = args.start
    windowend = args.start + dt.timedelta(hours=args.schedule_window)

    location = {
        'telescope_class': '1m0',
    }
    if args.site is not None:
        location['site'] = args.site
    if args.dome is not None:
        location['enclosure'] = args.dome
    if args.telescope is not None:
        location['telescope'] = args.telescope



    request = {'configurations': [],
               'windows': [{"start": absolutestart.isoformat(), "end": windowend.isoformat()}, ],
               'location': location}

    nresconfiguration = createNRESRequestsConfiguration(args)

    pm_ra = args.pm[0]
    pm_dec = args.pm[1]
    target = {
        "type": "ICRS",
        "name": f"NRES Commissioning {args.targetname}",
        "epoch": "2000.0000000",
        "equinox": "2000.0000000",
        "ra": "%10f" % args.radec.ra.degree,
        "dec": "%10f" % args.radec.dec.degree,
        "proper_motion_ra": pm_ra,
        "proper_motion_dec": pm_dec,
    }

    nresconfiguration['target'] = target
    nresconfiguration['constraints'] = common.default_constraints
    request['configurations'].append(nresconfiguration)
    requestgroup['requests'].append(request)
    return requestgroup


def convert_request_for_direct_submission(nres, args):
    """ Need to fill in some extra fields for a direct submission"""
    SLEWTIME = 180
    ACQUISITIONTIME = 300

    READOUTTIME = 60
    # Calculate the4 end tim efor this reuquest
    endtime = args.start + dt.timedelta(seconds=SLEWTIME + ACQUISITIONTIME)
    endtime += dt.timedelta(seconds=args.expcnt * (args.exptime + READOUTTIME))

    data = {
        'name': f'NRES ENGINEERING {args.targetname}',
        'proposal': args.proposalid,
        'site': args.site,
        'enclosure': args.dome,
        'telescope': args.telescope,
        'start': args.start.isoformat(),
        'end': endtime.isoformat(),
        'request': {
            'acceptability_threshold': 90,
            'configurations': nres['requests'][0]['configurations']
        }

    }
    return data


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Submit an engineering NRES observation to SCHEDULER.')

    parser.add_argument("--proposalid", default="LCOEngineering")

    parser.add_argument('--targetname', type=str, default=None,
                        help='Name of star for NRES test observation. if none is given, or auto, Software will try to'
                             ' find a cataloged flux stadnard star. Name must resolve via Simbad; if resolve failes,'
                             ' program will exit.')

    parser.add_argument('--site', choices=common.lco_nres_sites, default=None,
                               help="To which site to submit")
    parser.add_argument('--dome', choices = ['doma','domb','domc'], default = None, help="Which dome. Important for direct submission or non-schedulable override")
    parser.add_argument('--telescope', choices = ['1m0a',], default = '1m0a', help="Which telescope. Important for direct submission or non-schedulable override")
    parser.add_argument('--commissioning', action='store_true', help="request observation to instrument 1m0-NRES-Commissioning" )

    parser.add_argument('--expcnt', type=int, dest="expcnt", default=1)
    parser.add_argument('--exptime', type=float, default=120)
    parser.add_argument('--forcewcs', action='store_true',
                        help='Force WCSW based acquistion')
    parser.add_argument('--window', default=3, type=int, help="scheduling window length")
    parser.add_argument('--ipp', default=1.0, help="IPP priority for block")
    parser.add_argument('--pm', type=float, nargs=2, default = None, help="proper motion RA DEC in marcsec / year")

    parser.add_argument('--start', default=None,
                        help="When to start NRES observation. If not given, defaults to \"NOW\"")
    parser.add_argument('--schedule-window', default=3, type=float,
                        help="How long after start should request be schedulable?")
    parser.add_argument('--direct', action='store_true',
                    help='If set, make a direct submission instead of a scheduler submission.')
    parser.add_argument('--user', default='daniel_harbeck', help="Which user name to use for submission")
    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, block will be submitted.')
    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()
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

    if args.targetname is None:
        args.targetname = common.get_auto_target(common.goodNRESFluxTargets, sitecode=args.site, starttime=args.start)
        _logger.info(f"Auto selecting target {args.targetname}")

    if args.targetname is None:
        _logger.error("No target given, giving up.")
        exit(1)

    if args.pm is None:
        if args.targetname in common.listofpm:
            args.pm = common.listofpm[args.targetname]
        else:
            args.pm = [0,0]



    astropy.coordinates.name_resolve.sesame_database.set("simbad")
    try:
        _logger.debug("Resolving target name")

        args.radec = SkyCoord.from_name(args.targetname, parse=True)
    except Exception as e:
        print("Resolving target name failed, giving up {}".format(e))
        exit(1)

    print("Resolved target %s at corodinates %s %s" % (args.targetname, args.radec.ra, args.radec.dec))
    return args


def main():
    args = parseCommandLine()
    requstgroup = createRequest(args)
    if not args.direct:
        common.send_request_to_portal(requstgroup, args.opt_confirmed)
    else:
        directrequest = convert_request_for_direct_submission(requstgroup,args)
        _logger.info(f"Attempting direct submission {directrequest['start']} {directrequest['end']}")
        _logger.debug(json.dumps(directrequest, indent=2))
        common.submit_request_group(directrequest, args.opt_confirmed)
    exit(0)


if __name__ == '__main__':
    main()
