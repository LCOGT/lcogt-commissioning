''' Tool to directly submit observation for Delta Rho Commissioning.'''
import argparse
import datetime as dt
import itertools
import json
import logging
import sys
import astropy.coordinates
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

import lcocommissioning.common.common as common

_log = logging.getLogger(__name__)


def starflat_dither_offsets(number_of_points, max_radius_arcsec, jitter_arcsec=0., random_seed=None):
    """ Dither offsets (east, north, in arcsec) for a star-flat / self-calibration data set.

    Points follow a Fermat (sunflower) spiral: uniform coverage of the disk of radius maxradius, with every point at a
    different radius and position angle, so that each star samples many distinct positions on the detector.
    The undithered pointing is visited at the start, middle and end to monitor transparency.

    :param number_of_points: number of spiral points, including the centre.
    :param max_radius_arcsec: radius of the outermost pointing [arcsec].
    :param jitter_arcsec: random extra offset per axis, uniform within +/- jitter_arcsec; breaks any regular grid.
    :param random_seed: random seed for the jitter, for a reproducible pattern.
    :return: list of (east, north) offsets in arcsec
    """
    random_generator = np.random.default_rng(random_seed)

    # Fermat spiral: point n sits at radius ~ sqrt(n) (equal area per point) and is rotated by the golden angle
    # relative to point n - 1, so no two points share a radius or a position angle
    golden_angle = np.pi * (3 - np.sqrt(5))
    spiral_point_index = np.arange(1, number_of_points)
    radius_arcsec = max_radius_arcsec * np.sqrt(spiral_point_index / (number_of_points - 1))
    position_angle = spiral_point_index * golden_angle

    # Position angle is measured from north towards east; add random jitter (east jitter drawn first, then north)
    number_of_spiral_points = len(spiral_point_index)
    east_offsets_arcsec = radius_arcsec * np.sin(position_angle) + random_generator.uniform(-jitter_arcsec, jitter_arcsec,
                                                                                             number_of_spiral_points)
    north_offsets_arcsec = radius_arcsec * np.cos(position_angle) + random_generator.uniform(-jitter_arcsec, jitter_arcsec,
                                                                                              number_of_spiral_points)
    spiral_offsets = list(zip(east_offsets_arcsec, north_offsets_arcsec))

    # Visit the undithered centre at the start, in the middle and at the end
    centre_offset = (0., 0.)
    midpoint = len(spiral_offsets) // 2
    return [centre_offset] + spiral_offsets[:midpoint] + [centre_offset] + spiral_offsets[midpoint:] + [centre_offset]


def offset_pointing(radec, east_offset_arcsec, north_offset_arcsec):
    """ Return the sky coordinate offset by (east, north) arcsec on the sky, tangent-plane correct incl. cos(dec)."""
    return radec.spherical_offsets_by(east_offset_arcsec * u.arcsec, north_offset_arcsec * u.arcsec)

def create_cdk_request_configuration(args):
    print ((args.filter, args.defocus, args.exp_time))
    configuration = {
        'type': None,
        'instrument_type': args.insttype,
        'guiding_config': {'mode': 'OFF'}  if args.noselfguide else {'mode': 'ON', 'optional': True},
        'acquisition_config': {},
        'instrument_configs': [ {
            'exposure_count': 1 if args.exp_cnt is None else args.exp_cnt,
            'exposure_time': exptime,
            'mode': 'full_frame' if args.readmode is None else args.readmode,
            'optical_elements': {
                'filter': filter,
            },
            'extra_params': {
                'bin_x': 1,
                'bin_y': 1,
                'defocus': defocus,
                'offset_ra': args.offsetRA,
                'offset_dec': args.offsetDec,
            }

        } for filter,defocus,exptime in itertools.product(args.filter, args.defocus, args.exp_time) ],
        'extra_params' : {}

    }


    # if args.stackexptime >0:
    #     configuration['extra_params']['sub_expose'] = True
    #     configuration['extra_params']['sub_exposure_time'] = (int) (args.stackexptime*1000)

    if args.instrument:
        configuration['instrument_name'] = args.instrument

    if args.exp_cnt:
        configuration['type'] = 'EXPOSE'
        configuration['instrument_configs'][0]['exposure_count'] = args.exp_cnt

    elif args.filltime:
        configuration['type'] = 'REPEAT_EXPOSE'
        configuration['repeat_duration'] = args.filltime

    return configuration


def create_icrs_target(args, radec, suffix=''):
    return {
        "type": "ICRS",
        "name": " ".join(filter(None, [args.title, args.targetname, suffix])),
        "epoch": "2000.0000000",
        "equinox": "2000.0000000",
        "ra": "%10f" % radec.ra.degree,
        "dec": "%10f" % radec.dec.degree,
        'proper_motion_ra': args.pm[0],
        'proper_motion_dec': args.pm[1],
    }


def create_request(args):
    '''
    Create a schedulable delta Rho / CDK request
    :param args:
    :return:
    '''

    requestgroup = {"name": args.title,
                    "proposal": args.proposalid,
                    "ipp_value": args.ipp,
                    "operator": "SINGLE",  # "MANY" if args.dither else "SINGLE",
                    "observation_type": "NORMAL",
                    "requests": []
                    }

    absolutestart = args.start
    windowend = args.start + dt.timedelta(hours=args.schedule_window)
    if args.filltime:
        windowend = absolutestart + dt.timedelta(seconds=args.filltime + 180)

    location = {'telescope': args.telescope,
                'telescope_class': '1m0',
                'enclosure': args.dome,
                'site': args.site, }

    request = {'configurations': [],
               'windows': [{"start": absolutestart.isoformat(), "end": windowend.isoformat()}, ],
               'location': location}

    if args.dither:
        # One configuration per dither point, each with its own target coordinate. Unlike the
        # offset_ra/offset_dec instrument_config extra_params, this is not limited to 40% of the guider field of view.
        dither_offsets = starflat_dither_offsets(args.dither_npoints, args.dither_radius, args.dither_jitter, args.dither_seed)
        for dither_index, (east_offset_arcsec, north_offset_arcsec) in enumerate(dither_offsets):
            dither_pointing = offset_pointing(args.radec, east_offset_arcsec, north_offset_arcsec)
            dither_label = f"dither {dither_index:02d} E{east_offset_arcsec:+.0f} N{north_offset_arcsec:+.0f}"
            configuration = create_cdk_request_configuration(args)
            configuration['target'] = create_icrs_target(args, dither_pointing, dither_label)
            configuration['constraints'] = common.default_constraints
            request['configurations'].append(configuration)
            _log.info(f"dither {dither_index:02d}: offset E {east_offset_arcsec:+7.1f}\" N {north_offset_arcsec:+7.1f}\" "
                      f"-> {dither_pointing.to_string('hmsdms', precision=1)}")
        requestgroup['requests'].append(request)
        return requestgroup

    cdkconfiguration = create_cdk_request_configuration(args)

    if not args.transientmode:
        target = create_icrs_target(args, args.radec)
    else:
        nh = args.site.upper() in'ELPOGGTFNTLVNGQ'
        print ('hemisphere',nh)
        target = {
            "type": "SATELLITE",
            "name": "{} {}".format(args.title, args.targetname),
            'altitude': 70.,
            'azimuth' :180. if nh else 0.,
            'diff_altitude_rate': 0.,
            'diff_pitch_rate': 0.,
            'diff_pitch_acceleration': 0,
            'diff_roll_rate': 0.,
            'diff_roll_acceleration': 0.,
            'diff_epoch_rate': 0.,
            'diff_azimuth_rate':0,
            'diff_epoch': 0,
            'diff_altitude_acceleration': 0,
            'diff_azimuth_acceleration': 0

        }

    cdkconfiguration['target'] = target
    cdkconfiguration['constraints'] = common.default_constraints
    request['configurations'].append(cdkconfiguration)
    requestgroup['requests'].append(request)

    return requestgroup


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Delta Rho @ LCO engineering commissioning submission tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    tgargetgroup  = parser.add_mutually_exclusive_group()
    tgargetgroup.add_argument('--targetname', default='auto', type=str,
                        help='Name of object to observe; will be resolved via simbad. Can be coordinates in the form of Jhhmmss+ddmmss')
    tgargetgroup.add_argument('--transientmode', action='store_true',

                              help='If set, observe at meridian only in drifting sky mode.')
    parser.add_argument('--title', default="Sophia commissioning", help="Descriptive title for observation request")
    parser.add_argument('--proposalid', default="ENG2026B-001", help="proposal ID")

    parser.add_argument('--site', default='elp', choices=['elp', 'cpt','tfn','coj','lsc'],
                        help="To which site to submit")

    parser.add_argument('--dome', default='doma', choices=['doma', 'domb', 'domc'])
    parser.add_argument('--tel', dest="telescope", default='1m0a', help="Telescope to use for the observation")
    parser.add_argument('--instrument', default=None, help="Camera to use, e.g., ep60. If not given, any camera of --insttype at the location.")
    parser.add_argument('--insttype', default='1M0-SCICAM-SINISTRO', help="Observation portal instrument type")

    parser.add_argument('--start', default=None, type=str,
                        help="Time to start observation. If not given, defaults to \"NOW\". Specify as YYYYMMDD HH:MM")

    parser.add_argument('--schedule-window', default=8, type=float,
                        help="How long after start should request be schedulable?")

    parser.add_argument('--dither', action='store_true',
                        help='Star-flat dither: observe the target at a spiral of pointings, one configuration per pointing.')
    parser.add_argument('--dither-npoints', type=int, default=25, help="Number of spiral dither points incl. centre; "
                        "the centre is visited at start, middle and end, so npoints+2 pointings in total.")
    parser.add_argument('--dither-radius', type=float, default=780, help="Radius of the outermost dither point [arcsec]")
    parser.add_argument('--dither-jitter', type=float, default=20, help="Random jitter per axis added to each dither point [arcsec]")
    parser.add_argument('--dither-seed', type=int, default=42, help="Random seed for the dither jitter")
    parser.add_argument('--pointing-overhead', type=float, default=45,
                        help="Slew/settle overhead per dither pointing [s], used for the duration of a direct submission")

    parser.add_argument('--defocus', type=float, default=[0.0], nargs="*",  help="Amount to defocus star.")

    parser.add_argument('--filter', default=['rp'], nargs="*", choices=['opaque', 'w', 'up', 'gp', 'rp', 'ip', 'zs', 'Y',  'U','B', 'V', 'I'],
                        help="Select optical element filter")

    parser.add_argument('--exptimes', type=float, dest='exp_time',  default=[10.0], nargs="*",
                        help='Exposure time[s]')

    parser.add_argument('--ipp', type=float, default=1.0, help="ipp value")

    parser.add_argument('--offsetRA', default=0, help="Extra pointing offset to apply R.A.")
    parser.add_argument('--offsetDec', default=0, help="Extra pointing offset to apply Dec")
    parser.add_argument('--pm', default=[0., 0.], nargs=2, type=float, help="Proper motions RA, Dec, in mas/yr")

    repeatgroup = parser.add_mutually_exclusive_group()
    repeatgroup.add_argument('--expcnt', dest='exp_cnt', type=int, help="How often to repeat each exposure")
    repeatgroup.add_argument('--filltime', dest='filltime', type=float, help="How long to repeat exposures (seconds)")

    parser.add_argument('--readmode', type=str, default='full_frame', )
    #parser.add_argument('--stackexptime', type=float, default=0, help=" Enter self-guided image stacking mode with a sub-exposure time defined in seconds. No stacking if time is 0 seconds. ")
    parser.add_argument('--direct', action='store_true',
                         help='If set, submit directly instead of via scheduler.')
    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, observation will be submitted. If omitted, nothing will be submitted.')
    parser.add_argument('--noselfguide', action='store_true',
                        help='If set, do not selfguide.')
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

    if not args.transientmode:
        if  ('auto' in args.targetname):
            # automatically find the best target
            args.targetname = common.get_auto_target(common.goodXTalkTargets, args.site, args.start, moonseparation=40)
            if args.targetname is None:
                _log.error("Could not find a suitable auto target. Exiting.")
                sys.exit(1)
        try:
            if ('moon' in args.targetname):
                long, lat = common.lco_site_lonlat[args.site]
                alt = common.lco_site_alt[args.site]
                args.radec = astropy.coordinates.get_moon(time = astropy.time.Time(args.start), location = astropy.coordinates.EarthLocation.from_geodetic(lat=lat, lon=long, height=alt))
            elif args.targetname:
                _log.debug("Resolving target name")
                args.radec = SkyCoord.from_name(args.targetname, parse=True)

            print(f"Resolved target >{args.targetname}< at coordinates {args.radec.ra} {args.radec.dec}")
        except:
            _log.exception("Resolving target name failed, giving up")
            sys.exit(1)
    else:
        args.targetname = "Meridian Sky Drift"


    print (f"Submitting to {args.site} {args.dome} {args.telescope}")

    if not (args.exp_cnt or args.filltime):
        print("No exposure mode chosen, defaulting to one single EXPOSE")
        args.exp_cnt = 1

    return args


def ammend_request_for_direct_submission(cdk_request, args):
    """ Need to fill in some extra fields for a direct submission"""
    SLEWTIME = 120  # Driven by rotator which is excruciatingly slow.
    READOUTTIME = 6
    # Calculate the end time for this request
    end_time = args.start + dt.timedelta(seconds=SLEWTIME)

    configurations = cdk_request['requests'][0]['configurations']
    if args.filltime:
        end_time += dt.timedelta(seconds=args.filltime)
    if args.exp_cnt:
        for configuration in configurations:
            end_time += dt.timedelta(seconds=sum(ic['exposure_count'] * (ic['exposure_time'] + READOUTTIME)
                                                 for ic in configuration['instrument_configs']))
    if args.dither:
        end_time += dt.timedelta(seconds=len(configurations) * args.pointing_overhead)

    data = {
        'name': args.title,
        'proposal': args.proposalid,
        'site': args.site,
        'enclosure': args.dome,
        'telescope': args.telescope,
        'start': args.start.isoformat(),
        'end': end_time.isoformat(),
        'request': {
            'acceptability_threshold': 90,
            'configurations': configurations
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
    if args.direct:
        cdk = ammend_request_for_direct_submission(cdk, args)
        _log.info(f"Modifying for direct submission {cdk['start']} {cdk['end']}")
        _log.debug(json.dumps(cdk, indent=2))
        common.submit_request_group(cdk, args.opt_confirmed)
    else:
        _log.debug(json.dumps(cdk, indent=2))
        common.send_request_to_portal(cdk, args.opt_confirmed, )

if __name__ == '__main__':
    main()
