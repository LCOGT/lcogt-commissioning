from opensearchpy import OpenSearch
from opensearch_dsl import Search
import numpy as np
from astropy.table import Table
import sys
import matplotlib.pyplot as plt
import logging
import astropy.time as astt

log = logging.getLogger(__name__)
logging.getLogger('opensearch').setLevel(logging.WARNING)
logging.getLogger('connectionpool').setLevel(logging.WARNING)


def query_opensearch(opensearch_url='https://opensearch.lco.global', index='fitsheaders'):
    client = OpenSearch(opensearch_url)

    sourcecolumn = ['ALTITUDE', 'L1FWHM', 'L1ELLIP', 'AZIMUTH', 'DATE-OBS', 'ORIGNAME', 'FILTER',
                    'WINDSPEE', 'WINDDIR', 'INSTRUME', 'EXPTIME', 'L1ELLIPA']

    query = 'INSTRUME:ep03 AND SITEID:ogg AND ENCID:clma AND _exists_:L1FWHM'
    body = {
        'size': 10000,
        "_source": sourcecolumn,
        'query': {
            'query_string' : {'query': query}
        },
    }

    r = client.search(body=body, index=index,request_timeout=120 )
    intermediate= [r['_source'] for r in r['hits']['hits']]
    t = [[item[col] for col in sourcecolumn] for item in intermediate]
    t = np.asarray(t)

    try:
        dtypes = [float for ii in range (len(sourcecolumn))]
        dtypes[sourcecolumn.index('DATE-OBS')] = str
        dtypes[sourcecolumn.index('ORIGNAME')] = str
        dtypes[sourcecolumn.index('FILTER')] = str
        dtypes[sourcecolumn.index('INSTRUME')] = str
        t = Table(t,names=sourcecolumn, dtype=dtypes)
        t['DATE-OBS'] = astt.Time(t['DATE-OBS'], scale='utc', format=None).to_datetime()
    except:
        log.exception ("Cannot parse opensearch return")
        return None
    return t



def plotthings (data):

    plt.clf()
    plt.plot (t['DATE-OBS'], t['L1ELLIP'], '.')
    plt.ylim([0,1])
    plt.savefig ('date_el.png')


    plt.clf()
    plt.plot (t['WINDSPEE'], t['L1ELLIP'], '.')
    plt.xlabel("Wind Speed")
    plt.ylabel ("Ellipticity")
    plt.ylim([0,1])
    plt.savefig ('windspee_el.png')

    plt.clf()
    plt.plot (t['WINDDIR'], t['L1ELLIP'], '.')

    plt.xlabel("Wind Direction")
    plt.ylabel ("Ellipticity")
    plt.ylim([0,1])
    plt.savefig ('winddir_el.png')

    plt.clf()
    plt.plot (t['EXPTIME'], t['L1ELLIP'], '.')
    plt.xlabel("Exposure Time")
    plt.ylabel ("Ellipticity")
    plt.ylim([0,1])

    plt.savefig ('exptime_el.png')


    plt.clf()
    plt.plot (t['L1ELLIPA'], t['L1ELLIP'], '.')
    plt.xlabel("Orientation of Ellipticity")
    plt.ylabel ("Ellipticity")
    plt.ylim([0,1])

    plt.savefig ('ellipa_el.png')


    plt.clf()
    plt.plot (t['DATE-OBS'], t['L1ELLIPA'], '.')
    plt.xlabel("DateOBS")
    plt.ylabel ("Orientation of Ellipticity")
    plt.savefig ('dateobs-ellipa.png')


    plt.clf()
    plt.plot (t['ALTITUDE'], t['L1ELLIP'], '.')
    plt.xlabel("ALTITUDE")
    plt.ylabel ("Ellipticity")
    plt.ylim([0,1])
    plt.savefig ('alt_el.png')

    plt.clf()
    plt.plot (t['AZIMUTH'], t['L1ELLIP'], '.')
    plt.xlabel("AZ")
    plt.ylabel ("Ellipticity")
    plt.ylim([0,1])
    plt.savefig ('az-el.png')




t = query_opensearch()

plotthings(t)