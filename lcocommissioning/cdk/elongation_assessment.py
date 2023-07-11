import argparse
import datetime
import logging

import astropy.time as astt
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Angle
from astropy.table import Table
from opensearchpy import OpenSearch
from scipy.stats import gaussian_kde

log = logging.getLogger(__name__)
logging.getLogger('opensearch').setLevel(logging.WARNING)
logging.getLogger('connectionpool').setLevel(logging.WARNING)


def query_opensearch(opensearch_url='https://opensearch.lco.global', index='fitsheaders', site='cpt',enc='aqwa',instrument='sq38', telid='*', dateobs="2022-01-00", ndays=1):
    client = OpenSearch(opensearch_url)

    sourcecolumn = ['ALTITUDE', 'L1FWHM', 'L1ELLIP', 'AZIMUTH', 'DATE-OBS', 'ORIGNAME', 'FILTER',
                    'WINDSPEE', 'WINDDIR', 'INSTRUME', 'EXPTIME', 'L1ELLIPA',
                    'RA', 'DEC', 'LST', ]
    dateformat = '%Y%m%d'

    #ogg 0m4b preload change 2023-06-20
    query = f'INSTRUME:{instrument} AND SITEID:{site} AND ENCID:{enc} AND L1FWHM:* AND L1ELLIP:* AND TELID:{telid} AND' \
            f' DATE-OBS:["{dateobs.isoformat()}" TO "{ (dateobs+datetime.timedelta(days=ndays)).isoformat ()}"]'
    print (query)
    body = {
        'size': 10000,
        "_source": sourcecolumn,
        'query': {
            'query_string': {'query': query}
        },
    }

    r = client.search(body=body, index=index, request_timeout=1200)
    intermediate = [x['_source'] for x in r['hits']['hits']]
    t = [[item[col] for col in sourcecolumn] for item in intermediate]
    t = np.asarray(t)

    try:
        dtypes = [float for ii in range(len(sourcecolumn))]
        dtypes[sourcecolumn.index('DATE-OBS')] = str
        dtypes[sourcecolumn.index('ORIGNAME')] = str
        dtypes[sourcecolumn.index('FILTER')] = str
        dtypes[sourcecolumn.index('INSTRUME')] = str
        dtypes[sourcecolumn.index('RA')] = Angle
        dtypes[sourcecolumn.index('DEC')] = Angle
        dtypes[sourcecolumn.index('LST')] = Angle

        t = Table(t, names=sourcecolumn, dtype=dtypes)
        t['DATE-OBS'] = astt.Time(t['DATE-OBS'], scale='utc', format=None).to_datetime()
        t['RA'] = [Angle(x, unit=u.hourangle).to_value() for x in t['RA']]
        t['DEC'] = [Angle(x, unit=u.deg).to_value() for x in t['DEC']]
        t['LST'] = [Angle(x, unit=u.hourangle).to_value() for x in t['LST']]
        t['HA'] = t['RA'] - t['LST']
    except:
        log.exception("Cannot parse opensearch return")

        return None
    print (len(t))
    return t


def plotthings(t, site, enc, instrument, telid, dateobs, ndays):
    plt.clf()
    plt.plot(t['DATE-OBS'], t['L1ELLIP'], '.')
    plt.ylim([0, 1])
    plt.savefig('date_el.png')

    plt.clf()
    plt.scatter(t['WINDSPEE'], t['L1ELLIP'], c=t['WINDDIR'], s=2)
    plt.xlabel("Wind Speed")
    plt.ylabel("Ellipticity")
    plt.colorbar(label='Wind Direction [\deg]')
    plt.ylim([0, 1])
    plt.title (f'{site}-{enc}-{telid}--{instrument}-{dateobs.isoformat()}--{ndays}')

    plt.savefig(f'{site}-{enc}-{telid}--{instrument}-{dateobs.isoformat()}--{ndays}-windspee_el.png', bbox_inches='tight')

    plt.clf()

    plt.scatter(t['WINDDIR'], t['L1ELLIP'], c=t['WINDSPEE'], s=2)

    plt.xlabel("Wind Direction")
    plt.ylabel("Ellipticity")
    plt.colorbar(label='Windspeed')
    plt.ylim([0, 1])
    plt.title (f'{site}-{enc}-{telid}--{instrument}-{dateobs.isoformat()}--{ndays}')
    plt.savefig(f'{site}-{enc}-{telid}--{instrument}-{dateobs.isoformat()}--{ndays}-winddir_el.png', bbox_inches='tight')

    plt.clf()
    plt.plot(t['EXPTIME'], t['L1ELLIP'], '.')
    plt.xlabel("Exposure Time")
    plt.ylabel("Ellipticity")
    plt.ylim([0, 1])

    plt.savefig('exptime_el.png')

    plt.clf()
    plt.plot(t['L1ELLIPA'], t['L1ELLIP'], '.')
    plt.xlabel("Orientation of Ellipticity")
    plt.ylabel("Ellipticity")
    plt.ylim([0, 1])

    plt.savefig(f'{site}-{enc}-{telid}--{instrument}-{dateobs.isoformat()}--{ndays}-ellipa_el.png')

    plt.clf()
    plt.plot(t['DATE-OBS'], t['L1ELLIPA'], '.')
    plt.xlabel("DateOBS")
    plt.ylabel("Orientation of Ellipticity")
    plt.savefig(f'{site}-{enc}-{telid}--{instrument}-{dateobs.isoformat()}--{ndays}-dateobs-ellipa.png')

    # plt.clf()
    # xy = np.vstack([t['DEC'], t['L1ELLIP']])
    # z = gaussian_kde(xy)(xy)
    # plt.scatter(t['DEC'], t['L1ELLIP'], c=z, s=1)
    # plt.xlabel("DEC")
    # plt.ylabel("Ellipticity")
    # plt.ylim([0, 1])
    # plt.savefig('dec_el.png')

    plt.clf()
    xy = np.vstack([t['HA'], t['L1ELLIP']])
    z = gaussian_kde(xy)(xy)
    # plt.plot (t['HA'], t['L1ELLIP'], '.')
    plt.scatter(t['HA'], t['L1ELLIP'], c=z, s=1)
    plt.xlabel("HA")
    plt.ylabel("Ellipticity")
    plt.ylim([0, 1])
    plt.xlim([-6, 6])
    plt.savefig(f'{site}-{enc}-{telid}--{instrument}-{dateobs.isoformat()}--{ndays}-ha-el.png')

    plt.clf()
    plt.plot(t['AZIMUTH'], t['L1ELLIP'], '.')
    plt.xlabel("AZ")
    plt.ylabel("Ellipticity")
    plt.ylim([0, 1])
    plt.savefig(f'{site}-{enc}-{telid}--{instrument}-{dateobs.isoformat()}--{ndays}-az-el.png')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--site', default='ogg', choices=['ogg', 'elp', 'cpt'],
                    help="To which site to submit")

    parser.add_argument('--dome', default='clma', choices=['aqwa', 'aqwb', 'clma'])
    parser.add_argument('--telescope', default='0m4c', choices=['0m4a','0m4b','0m4c'])
    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                    help='Set the debug level')
    parser.add_argument('--date', type=lambda d: datetime.datetime.strptime(d, '%Y%m%d'), default=datetime.date.today(),
                    help="Date-OBS to start. If not given, defaults to \"NOW\". Specify as YYYYMMDD")
    parser.add_argument('--ndays', type = int, default=1)
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                    format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    instrument='*'
    t = query_opensearch(enc=args.dome,site=args.site,instrument=instrument,telid = args.telescope, dateobs=args.date, ndays = args.ndays)
    plotthings(t,enc=args.dome,site=args.site,instrument=instrument, telid = args.telescope, dateobs=args.date, ndays = args.ndays)



if __name__ == '__main__':
    main()