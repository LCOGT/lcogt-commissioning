''' Utility to crawl end of focus sequence images and use their focus
    values to analyze focus compensation terms.
    the basic assumption is that at the end of the auto_focus sequence, an
    image is taken at the best possible curent fous for the current telescope's
    temperature and zenith distance.
    '''

import argparse
import datetime as dt
import json
import logging
import math

import astropy.time as astt
import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy.table import Table
from scipy.optimize import curve_fit

logging.getLogger('matplotlib').setLevel(logging.FATAL)

plt.rcParams["figure.figsize"] = (20, 34)
plt.style.use('ggplot')

log = logging.getLogger(__name__)

telescopedict = {
    'lsc': ['doma-1m0a', 'domb-1m0a', 'domc-1m0a', 'aqwa-0m4a', 'aqwb-0m4a'],
    'coj': ['clma-2m0a', 'doma-1m0a', 'domb-1m0a', 'clma-0m4a', 'clma-0m4b'],
    'ogg': ['clma-2m0a', 'clma-0m4b', 'clma-0m4c'],
    'elp': ['doma-1m0a', 'domb-1m0a', 'aqwa-0m4a'],
    'cpt': ['doma-1m0a', 'domb-1m0a', 'domc-1m0a', 'aqwa-0m4a'],
    'tfn': ['aqwa-0m4a', 'aqwa-0m4b', 'doma-1m0a', 'domb-1m0a'],
}

wellintonighttime = {
    'lsc': dt.timedelta(hours=26),
    'coj': dt.timedelta(hours=11),
    'ogg': dt.timedelta(hours=24 + 9),
    'elp': dt.timedelta(hours=28),
    'cpt': dt.timedelta(hours=19),
    'tfn': dt.timedelta(hours=23)
}


def robustfit(x, y, sigma=3, iterations=3, label=None):
    ''' Do a sigma clipping linear fit.
    TODO: Consider moving this over to an astropy routine that already does this
    '''
    _x = x
    _y = y
    if len(_x) < 3:
        return None
    for iteration in range(iterations):
        pf = np.polyfit(_x, _y, 1)
        p = np.poly1d(pf)
        residual = _y - p(_x)
        std = np.std(residual)
        good = np.abs(residual) < sigma * std
        _x = _x[good]
        _y = _y[good]
    xp = np.linspace(np.min(x), np.max(x), 10)
    plt.plot(xp, p(xp), label=p if label is not None else None)
    return p


def calulate_compensation(zdinput, tempinput, fitresult):
    ''' Convenience method to use a multi-linear fit result and apply it to data. '''
    return fitresult[0] + fitresult[1] * zdinput + fitresult[2] * tempinput


def simultaneousfit(temp, coszd, focus):
    ''' Do  A MULTI-LINEAR FIT TO SIMULTANEOUSLY FIT THE TEMP AND ZD DEPENDENCE
    OF THE focus position. '''
    def fn(x, a, cterm, tterm):
        return a + cterm * x[0] + tterm * x[1]

    np.set_printoptions(formatter={'float_kind': '{:6.4f}'.format})

    good_coszd = coszd
    good_temp = temp
    good_focus = focus

    for iter in range (3):
        # 3 sigma clipping of fit result to remove outliers
        popt, pcov, = curve_fit(fn, [good_coszd, good_temp], good_focus)
        errors = np.sqrt(np.diag(pcov))
        log.info(f'Multilinear fit result, N={len(good_focus)}: {popt} +/-  {errors}')

        residual = good_focus - calulate_compensation(good_coszd, good_temp, popt)
        goodvalues = np.abs(residual) < 3 * np.std (residual)
        good_coszd = good_coszd [goodvalues]
        good_temp = good_temp [goodvalues]
        good_focus = good_focus[goodvalues]

    return popt, errors


def get_focusStackData(args):
    """Get focus-relevant fits header keywords for end of focus exposures,
       as a proxy for the idealy focused telescope state """

    site = args.site
    enc = args.dome
    tel = args.tel

    sourcecolumn = ['FOCTEMP', 'ALTITUDE', 'FOCTELZP', 'FOCTOFF', 'FOCZOFF',
                    'FOCAFOFF', 'FOCINOFF', 'FOCFLOFF', 'L1FWHM', 'FOCOBOFF',
                    'AZIMUTH', 'REFHUMID', 'DATE-OBS', 'DAY-OBS', 'ORIGNAME',
                    'WMSTEMP', 'FILTER']

    query = f'SITEID:{site} AND ENCID:{enc} AND TELID:{tel} AND OBJECT:auto_focus' \
            f' AND FOCOBOFF:0 AND RLEVEL:91 AND OBSTYPE:EXPOSE AND _exists_:L1FWHM'
    log.debug(f"Query String: {query}")

    # due to cameras moving around, and mirrors being replaced, autofocus values are
    # informative only over a limited date range.
    bestaftertime = args.after
    log.info(f"Selecting after date: {bestaftertime}")

    body = {
        "size": 10000,
        "_source": sourcecolumn,
        "query": {
            "query_string": {
                "query": query
            }
        }
    }

    headers = {'Content-type': 'application/json'}
    response = requests.get('https://opensearch.lco.global:443/fitsheaders/_search',
                            data=json.dumps(body), headers=headers).json()
    # make a nice table out of that convoluted mess
    intermediate = [r['_source'] for r in response['hits']['hits']]
    t = [[item[col] for col in sourcecolumn] for item in intermediate]
    t = np.asarray(t)
    log.info(f"Number of raw records {len(t)}")
    if len(t) == 0:
        log.warning("Warning: empty return for {} {} {}".format(site, enc, tel))
        return None

    # reformatting boilerplate stuff
    dayobsidx = sourcecolumn.index('DAY-OBS')
    for ii in range(len(t)):
        t[ii, dayobsidx] = "-".join([t[ii, dayobsidx][0:4], t[ii, dayobsidx][4:6],
                                     t[ii, dayobsidx][6:8]])

    try:
        dtypes = [float for ii in range(len(sourcecolumn))]
        dtypes[sourcecolumn.index('DATE-OBS')] = str
        dtypes[sourcecolumn.index('DAY-OBS')] = str
        dtypes[sourcecolumn.index('ORIGNAME')] = str
        dtypes[sourcecolumn.index('FILTER')] = str

        t = Table(t, names=sourcecolumn,
                  dtype=dtypes)
    except:
        log.exception("Error parsing opensearch return")
        return None

    if '1m0' in tel:
        t['FOCAFOFF'] = t['FOCAFOFF'] / 11.24
    if '0m4' in tel:
        t['FOCAFOFF'] = t['FOCAFOFF'] / 1.
    if '2m0' in tel:
        t['FOCAFOFF'] = t['FOCAFOFF'] / 12.05

    mediantelfoc = np.mean(t['FOCTELZP'])

    t['ACTFOCUS'] = t['FOCTELZP'] + t['FOCTOFF'] + t['FOCZOFF'] + t['FOCAFOFF'] \
                    + t['FOCINOFF'] + t['FOCFLOFF']

    t['ZD'] = 90 - t['ALTITUDE']
    t['DATE-OBS'] = astt.Time(t['DATE-OBS'], scale='utc', format=None).to_datetime()
    t['DAY-OBS'] = astt.Time(t['DAY-OBS'], scale='utc', format=None).to_datetime()

    t['FOCAFOFF_CORRECTED'] = t['FOCAFOFF'] + (t['FOCTELZP'] - mediantelfoc)

    # some rejection of bad values
    limit = 10
    if '2m0' in tel:
        limit = 320

    t.sort('DATE-OBS')

    good = (t['DATE-OBS'] > bestaftertime)
    good = good & (np.abs(t['ACTFOCUS']) < limit) & (t['FOCOBOFF'] == 0)
    t = t[good] if np.sum(good) > 0 else None
    log.info(f"Number of sanitized records: {len(t)}")
    return t


def set_ylim(data, yrange):
    ''' Based on the inputdata stats, set the y range of the plt plot. '''
    med = np.median(data)
    plt.ylim([med - yrange / 2, med + yrange / 2])


def analysecamera(args, t=None, ):
    ''' Do a full thermal / compression analysi for a given telescope'''
    site = args.site
    enc = args.dome
    tel = args.tel

    focustermrange = 0.3
    if '0m4' in tel:
        focustermrange = 0.1
    if '2m0' in tel:
        focustermrange = 0.2
    focusvaluerange = 0.75
    if '0m4' in tel:
        focusvaluerange = 0.1
    if '2m0' in tel:
        focusvaluerange = 3

    if t is None:
        t = get_focusStackData(args)

    # Multi-linear fit to find dependency on temeprature and zenith distance.
    coszd = np.cos(t['ZD'] * math.pi / 180)
    temp = t['FOCTEMP']
    simfitresult, simfiterrors = simultaneousfit(temp, coszd, t['ACTFOCUS'])
    fitted_focus = calulate_compensation(coszd, temp, simfitresult)
    residual_focus = t['ACTFOCUS'] - fitted_focus

    # Now let's get some diagnostic plots.

    ### Temperature
    plt.clf()
    plt.figure(figsize=(12, 20))
    plt.subplot(5, 2, 1)
    xdata = temp
    ydata = t['ACTFOCUS']
    plt.plot(xdata, ydata, '.', label="Measurements")
    xdata = np.arange(-5, 35, 0.05)
    ydata = simfitresult[0] + simfitresult[2] * xdata + np.mean(coszd) * simfitresult[1]
    plt.plot(xdata, ydata, '-', label=f"Temperature term {simfitresult[2]:6.4f} +/- {simfiterrors[2]:6.4f}")
    plt.xlabel('FOCTEMP [deg C]')
    plt.ylabel('FOCUS [mm]')
    plt.xlim([-6, 35])
    set_ylim(t['ACTFOCUS'], focusvaluerange)
    plt.legend()
    plt.title("Temp vs absolute focus postion")

    ### ZD
    plt.subplot(5, 2, 2)
    ydata = t['ACTFOCUS']
    plt.plot(coszd, ydata, '.')
    plt.xlabel('cos (ZD)')
    plt.ylabel('FOCUS [mm], T corrected')
    plt.xlim([0, 1.05])
    set_ylim(t['ACTFOCUS'], focusvaluerange)
    xdata = np.arange(0, 4.1, 0.05)
    ydata = simfitresult[0] + simfitresult[1] * xdata + np.mean(temp) * simfitresult[2]
    plt.plot(xdata, ydata, '-', label=f"Compression term {simfitresult[1]:6.4f} +/- {simfiterrors[1]:6.4f}")
    plt.legend()
    plt.title("Zenith distance vs abs focus position")


    # Temperature Residual after multilinear fit corrections
    plt.subplot(5, 2, 3)
    plt.title("Residual after temeprature and ZD correction")
    plt.plot(temp, residual_focus, '.')
    plt.xlabel('FOCTEMP [deg C]')
    plt.ylabel('ZD & Temp corrected Focus')
    plt.xlim([-6, 35])
    plt.ylim([-focustermrange /2, focustermrange /2])

    # Compression Residual after multilinear fit corrections
    plt.subplot(5, 2, 4)
    plt.title("Residual after temeprature and ZD correction")
    plt.plot(coszd, residual_focus, '.')
    plt.xlabel('cos(ZD)')
    plt.ylabel('ZD & Temp corrected Focus')
    plt.xlim([0, 1.05])
    plt.ylim([-focustermrange /2, focustermrange/2 ])

    plt.subplot(5, 1, 3)
    # Time line plot - history of autofocus corrections
    plt.plot(t['DATE-OBS'], t['FOCAFOFF'], '.', label="AUTOFOCUS Correction")
    plt.xlabel('DATE-OBS')
    plt.ylabel('AUTO FOCUS TERM [mm]')
    plt.title("History of AF corrections")
    plt.setp(plt.gca().xaxis.get_minorticklabels(), rotation=25)
    plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=25)
    set_ylim(t['FOCAFOFF'], 2*focusvaluerange)

    # no wlook how hard autofocus has to work to compenssate for temp & ZD
    # compensation should be noise around constant offset value for good
    # compensation.
    plt.subplot(5, 2, 8)
    xdata = coszd
    ydata = t['FOCAFOFF']
    plt.xlabel('Cos ZD')
    plt.ylabel('Auto focus corrections')
    plt.xlim([0, 1.05])
    set_ylim(t['FOCAFOFF'], 2*focusvaluerange)
    plt.plot(xdata, ydata, '.')

    plt.subplot(5, 2, 7)
    xdata = t['FOCTEMP']
    ydata = t['FOCAFOFF']

    plt.xlabel('Temperature')
    plt.ylabel('Auto focus corrections')
    plt.xlim([-6, 35])
    set_ylim(t['FOCAFOFF'], 2*focusvaluerange)

    plt.plot(xdata, ydata, '.')

    plt.suptitle("{} {} {} \n".format(site, enc, tel), fontsize=24, y=1.05)

    # Lok at some inttersting residuals
    plt.subplot(5, 2, 9)
    xdata = t['REFHUMID']
    ydata = residual_focus

    plt.xlabel('Humidity')
    plt.ylabel('Residual after temeprature and ZD correction')
    plt.plot(xdata, ydata, '.')
    set_ylim(residual_focus, focusvaluerange)
    #plt.ylim([-focustermrange / 5, focustermrange / 5])

    plt.subplot(5, 2, 10)
    xdata = t['AZIMUTH']
    ydata = residual_focus

    plt.xlabel('AZ [deg]')
    plt.ylabel('Residual after temeprature and ZD correction')
    plt.plot(xdata, ydata, '.')
    set_ylim(residual_focus, focustermrange)
    #plt.ylim([-focustermrange / 5, focustermrange / 5])

    plt.suptitle("{} {} {} \n".format(site, enc, tel), fontsize=24, y=1.05)

    plt.tight_layout(pad=1)
    plt.savefig(f"focusstack_{site}_{enc}_{tel}.pdf")
    return t


def getargs():
    parser = argparse.ArgumentParser(
        description='Analyst focus dependency on temeprature and stuff')

    parser.add_argument('--site', default='ogg', type=str)
    parser.add_argument('--dome', default='clma', type=str)
    parser.add_argument('--tel', default='0m4c', type=str)
    parser.add_argument('--after', type=dt.datetime.fromisoformat)

    parser.add_argument('--loglevel', dest='log_level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARN'], help='Set the debug level')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    return args


def main():
    args = getargs()

    print(f"Analysing focus stack data for {args.site} {args.dome} {args.tel}")
    analysecamera(args)


if __name__ == '__main__':
    main()
