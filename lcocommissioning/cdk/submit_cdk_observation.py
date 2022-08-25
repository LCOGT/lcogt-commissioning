''' Tool to directly submit observation for Delta Rho Commissioning.'''
import argparse
import datetime as dt
import json
import logging
import sys
import astropy.coordinates
from astropy.coordinates import SkyCoord

import lcocommissioning.common.common as common

_log = logging.getLogger(__name__)

def create_cdk_request_configuration(args):

    configuration = {
        'type': None,
        'instrument_type': '0m4-SCICAM-QHY600',
        'guiding_config': {'mode': 'ON', 'optional': True} if args.selfguide else {'mode': 'OFF'},
        'acquisition_config': {},
        'instrument_configs': [{
            'exposure_count': 1 if args.exp_cnt is None else args.exp_cnt,
            'exposure_time': args.exp_time,
            'mode': 'default' if args.readmode is None else args.readmode,
            'optical_elements': {
                'filter': args.filter,
            },
            'extra_params': {
                'bin_x': 1,
                'bin_y': 1,
                'defocus': args.defocus,
                'offset_ra': args.offsetRA,
                'offset_dec': args.offsetDec,
            }

        }, ]
    }

    # if args.selfguide:
    #     # several implemtnations tried out:
    #     # initially: one channal as guider camera:
    #     # configuration['guiding_config']['mode'] = f'MUSCAT_{args.selfguide.upper()}'
    #     # or: self-guide with an appropiate channel.
    #     configuration['guiding_config']['mode'] = 'OFF'

    if args.exp_cnt:
        configuration['type'] = 'EXPOSE'
        configuration['instrument_configs'][0]['exposure_count'] = args.exp_cnt

    elif args.filltime:
        configuration['type'] = 'REPEAT_EXPOSE'
        configuration['repeat_duration'] = args.filltime

    return configuration


def create_request(args):
    '''
    Create a schedulable delta Rho / CDK request
    :param args:
    :return:
    '''

    requestgroup = {"name": args.title,
                    "proposal": "DeltaRho Commissioning",
                    "ipp_value": args.ipp,
                    "operator": "SINGLE",  # "MANY" if args.dither else "SINGLE",
                    "observation_type": "NORMAL",
                    "requests": []
                    }

    absolutestart = args.start
    windowend = args.start + dt.timedelta(hours=args.schedule_window)

    location = {'telescope': '0m4c',
                'telescope_class': '0m4',
                'enclosure': 'clma',
                'site': args.site, }

    request = {'configurations': [],
               'windows': [{"start": absolutestart.isoformat(), "end": windowend.isoformat()}, ],
               'location': location}

    cdkconfiguration = create_cdk_request_configuration(args)

    target = {
        "type": "ICRS",
        "name": "{} {}".format(args.title, args.targetname),
        "epoch": "2000.0000000",
        "equinox": "2000.0000000",
        "ra": "%10f" % args.radec.ra.degree,
        "dec": "%10f" % args.radec.dec.degree,
        'proper_motion_ra': args.pm[0],
        'proper_motion_dec': args.pm[1],
    }

    cdkconfiguration['target'] = target
    cdkconfiguration['constraints'] = common.default_constraints
    request['configurations'].append(cdkconfiguration)
    requestgroup['requests'].append(request)

    return requestgroup


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='CDK17 @ LCO engineering commissioning submission tool')

    parser.add_argument('--targetname', default='auto', type=str,
                        help='Name of object to observe; will beresolved via simbad. Can be coordinates in the form of Jhhmmss+ddmmss')
    parser.add_argument('--title', default="CDK17 commissioning", help="Descriptive title for observation request")

    parser.add_argument('--site', default='ogg', choices=common.lco_2meter_sites,
                        help="To which site to submit")

    parser.add_argument('--start', default=None,
                        help="Time to start observation. If not given, defaults to \"NOW\". Specify as YYYYMMDD HH:MM")

    parser.add_argument('--schedule-window', default=2, type=float,
                        help="How long after start should request be schedulable?")

    # parser.add_argument('--dither', action='store_true', help='Dither the exposure in a 5 point pattern.')

    parser.add_argument('--defocus', type=float, default=0.0, help="Amount to defocus star.")

    parser.add_argument('--filter', default='rp', choices=['opaque', 'w', 'up', 'gp', 'rp', 'ip', 'zs', 'U', 'B', 'V'],
                        help="Select optical element filter")

    parser.add_argument('--exp-time', type=float, default=10,
                        help='Exposure time')

    parser.add_argument('--ipp', type=float, default=1.0, help="ipp value")

    parser.add_argument('--offsetRA', default=0, help="Extra pointing offset to apply R.A.")
    parser.add_argument('--offsetDec', default=0, help="Extra pointing offset to apply Dec")
    parser.add_argument('--pm', default=[0., 0.], nargs=2, type=float, help="Proper motions RA, Dec, in mas/yr")

    repeatgroup = parser.add_mutually_exclusive_group()
    repeatgroup.add_argument('--exp-cnt', type=int, help="How often to repeat each exposure")
    repeatgroup.add_argument('--filltime', type=float, help="How long to repeat exposures (seconds)")

    parser.add_argument('--readmode', type=str, default=None, )
    # parser.add_argument('--scheduler', action='store_true',
    #                     help='If set, submit to scheduler instead of trying a direct submission.')
    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, observation will be submitted. If omitted, nothing will be submitted.')
    parser.add_argument('--selfguide', action='store_true',
                        help='If set, selfguide.')
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
            _log.error(f"Invalid start time argument: {args.start}")
            sys.exit(1)

    if 'auto' in args.targetname:
        # automatically find the best target
        args.targetname = common.get_auto_target(common.goodXTalkTargets, args.site, args.start, moonseparation=40)
        if args.targetname is None:
            _log.error("Could not find a suitable auto target. Exiting.")
            sys.exit(1)


    try:
        if 'moon' in args.targetname:
            long, lat = common.lco_site_lonlat[args.site]
            alt = common.lco_site_alt[args.site]
            args.radec = astropy.coordinates.get_moon(time = astropy.time.Time(args.start), location = astropy.coordinates.EarthLocation.from_geodetic(lat=lat, lon=long, height=alt))
        else:
            _log.debug("Resolving target name")
            args.radec = SkyCoord.from_name(args.targetname, parse=True)
    except:
        _log.exception("Resolving target name failed, giving up")
        sys.exit(1)

    print(f"Resolved target >{args.targetname}< at coordinates {args.radec.ra} {args.radec.dec}")

    if not (args.exp_cnt or args.filltime):
        print("No exposure mode chosen, defaulting to EXPOSE")
        args.exp_cnt = 1

    return args


def ammend_request_for_direct_submission(cdk_request, args):
    """ Need to fill in some extra fields for a direct submission"""
    SLEWTIME = 120  # Driven by rotator which is excuiciatingly slow.
    READOUTTIME = 6
    # Calculate the end time for this request
    end_time = args.start + dt.timedelta(seconds=SLEWTIME)

    if args.filltime:
        end_time += dt.timedelta(seconds=args.filltime)
    if args.exp_cnt:
        end_time += dt.timedelta(seconds=args.exp_cnt * (float(args.exp_time) + READOUTTIME))

    data = {
        'name': args.title,
        'proposal': 'DeltaRho Commissioning',
        'site': args.site,
        'enclosure': 'clma',
        'telescope': '0m4c',
        'start': args.start.isoformat(),
        'end': end_time.isoformat(),
        'request': {
            'acceptability_threshold': 90,
            'configurations': cdk_request['requests'][0]['configurations']
        }

    }
    return data


def main():
    args = parseCommandLine()

    cdk = create_request(args)

    # if args.scheduler:
    #     _log.info("Submitting to scheduler")
    #     _log.debug(json.dumps(cdk, indent=2))
    #     common.send_request_to_portal(cdk, args.opt_confirmed)
    # else:
    cdk = ammend_request_for_direct_submission(cdk, args)
    _log.info(f"Attempting direct submission {cdk['start']} {cdk['end']}")
    _log.debug(json.dumps(cdk, indent=2))
    common.submit_request_group(cdk, args.opt_confirmed)


if __name__ == '__main__':
    main()
