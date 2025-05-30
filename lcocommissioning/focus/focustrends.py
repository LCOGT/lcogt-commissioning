''' Utility to crawl end of focus sequence images and use their focus
    values to analyze focus compensation terms.
    the basic assumption is that at the end of the auto_focus sequence, an
    image is taken at the best possible current focus for the current telescope's
    temperature and zenith distance.




    '''

import argparse
import datetime
import datetime as dt
import json
import logging
import math

import astropy.time as astt
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

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

# wellintonighttime = {
#     'lsc': dt.timedelta(hours=26),
#     'coj': dt.timedelta(hours=11),
#     'ogg': dt.timedelta(hours=24 + 9),
#     'elp': dt.timedelta(hours=28),
#     'cpt': dt.timedelta(hours=19),
#     'tfn': dt.timedelta(hours=23)
# }


focus_temp_reference_points  = {
    '1m0': 15,
    '0m4': 15,
   
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


def calulate_compensation(zdinput, tempinput, fitresult, tempref = 0):
    ''' Convenience method to use a multi-linear fit result and apply it to data. '''
    deltaT = tempinput - tempref
    result =  fitresult[0] + fitresult[1] * zdinput + fitresult[2] * deltaT  + fitresult[3] * deltaT**2
    return result

def simultaneousfit(temp, coszd, focus, zdterm = None, focus_temp_reference=0, order = 1):
    ''' Do  A MULTI-LINEAR FIT TO SIMULTANEOUSLY FIT THE TEMP AND ZD DEPENDENCE
    OF THE focus position. '''
    def fn(x, a, cterm, tterm, tterm2):
        deltaT = x[1] - focus_temp_reference
        return  a + cterm * x[0] + tterm * deltaT + tterm2 * deltaT**2
    
    
    np.set_printoptions(formatter={'float_kind': '{:6.4f}'.format})

    zdmin = -math.inf if zdterm is None else np.nextafter(-zdterm, -1) 
    zdmax = math.inf if zdterm is None else -zdterm
    tterm2low = -math.inf if order == 2 else 0  
    tterm2high = math.inf if order == 2 else np.nextafter(0, 1)
       
    bounds = ([-math.inf,zdmin, -math.inf, tterm2low], [math.inf, zdmax, math.inf, tterm2high])
    
    goodvalues =[True for  ii in range (len(focus))]


    for iter in range (5):
        popt, pcov, = curve_fit(fn, [coszd[goodvalues], temp[goodvalues]], focus[goodvalues], bounds=bounds)
        errors = np.sqrt(np.diag(pcov))  
        residual = focus - calulate_compensation(coszd, temp, popt, tempref=focus_temp_reference)
        residualrms = np.std (residual[goodvalues])
        log.info(f'Multilinear fit result, N={np.sum (goodvalues)}: {popt} +/-  {errors}, measured rms {residualrms:5.2f}')
        goodvalues = np.abs(residual) < 2.2 * residualrms
   
    return popt, errors, goodvalues


def get_focusStackData(args):
    """Get focus-relevant fits header keywords for end of focus exposures,
       as a proxy for the ideally focused telescope state """

    site = args.site
    enc = args.dome
    tel = args.tel

    sourcecolumn = ['FOCTEMP', 'ALTITUDE', 'FOCTELZP', 'FOCTOFF', 'FOCZOFF',
                    'FOCAFOFF', 'FOCINOFF', 'FOCFLOFF', 'FOCOBOFF', 'FOCPOSN',
                    'AZIMUTH', 'REFHUMID', 'DATE-OBS', 'DAY-OBS', 'ORIGNAME',
                    'WMSTEMP', 'FILTER', 'CCDATEMP', 'REQNUM', 'L1FWHM', 'SUNALT']

    # due to cameras moving around, and mirrors being replaced, autofocus values are
    # informative only over a limited date range.
    bestaftertime = args.after
    log.info(f"Selecting after date: {bestaftertime}")

    query = f'SITEID:{site} AND ENCID:{enc} AND TELID:{tel} AND OBJECT:auto_focus' \
            f' AND FOCOBOFF:0 AND RLEVEL:91 AND OBSTYPE:EXPOSE AND FILTER:rp' \
            f' AND L1FWHM:[0.5 TO {args.maxfwhm}] AND SUNALT:[-180 TO {args.maxsunalt}]' \
            f' AND DATE-OBS:[{bestaftertime.strftime("%Y-%m-%dT%H:%M:%S")} TO {args.before.strftime("%Y-%m-%dT%H:%M:%S")}]'



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
        dtypes[sourcecolumn.index('REQNUM')] = str
        t = Table(t, names=sourcecolumn,
                  dtype=dtypes)
    except:
        log.exception("Error parsing opensearch return")
        return None

    if '1m0' in tel:
        magnification = 1+3.2**2
        #t['FOCAFOFF'] = t['FOCAFOFF'] / 11.24
    if '0m4' in tel:
        magnification = 1
        #t['FOCAFOFF'] = t['FOCAFOFF'] / 1.
    if '2m0' in tel:
        #t['FOCAFOFF'] = t['FOCAFOFF'] / 12.05
        magnification = 1+3.3239**2

    t['FOCAFOFF'] = t['FOCAFOFF'] / magnification
    t['FOCFLOFF'] = t['FOCFLOFF'] / magnification
    t['ACTFOCUS'] = t['FOCTELZP'] + t['FOCTOFF'] + t['FOCZOFF'] + t['FOCAFOFF']  \
                    + t['FOCINOFF'] + t['FOCFLOFF']

    t['ZD'] = 90 - t['ALTITUDE']
    t['DATE-OBS'] = astt.Time(t['DATE-OBS'], scale='utc', format=None).to_datetime()
    t['DAY-OBS'] = astt.Time(t['DAY-OBS'], scale='utc', format=None).to_datetime()

    # some rejection of bad values
    limit = 10
    if '2m0' in tel:
        limit = 320
    if '1m0' in tel:
        limit = 5

    fullsize = len(t)
    t.sort('DATE-OBS')
    good = (np.abs(t['ACTFOCUS']) < limit) & (t['FOCOBOFF'] == 0)
    good = good & (t['FOCTEMP'] != 0.0)

    t = t[good] if np.sum(good) > 0 else None
    log.info(f"Number of sanitized records: {len(t)} / {fullsize}")
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
        focusvaluerange = 0.3
    if '2m0' in tel:
        focusvaluerange = 3

    if t is None:
        t = get_focusStackData(args)

    # Multi-linear fit to find dependency on temeprature and zenith distance.
    coszd = np.cos(t['ZD'] * math.pi / 180)
    temp = t['FOCTEMP']
    simfitresult, simfiterrors, goodvalues = simultaneousfit(temp, coszd, t['ACTFOCUS'], zdterm=0.15 if "1m0" in tel else None, order=args.fitorder, focus_temp_reference=args.focus_temp_reference)
    fitted_focus = calulate_compensation(coszd, temp, simfitresult, tempref=args.focus_temp_reference)
    residual_focus = t['ACTFOCUS'] - fitted_focus

    # Now let's get some diagnostic plots.

    ### Temperature
    plt.clf()
    plt.figure(figsize=(12, 20))
    plt.subplot(5, 2, 1)
    xdata = temp
    ydata = t['ACTFOCUS']
    plt.plot(xdata, ydata, '.', c='grey')
    plt.plot(xdata[goodvalues], ydata[goodvalues], '.', c='red')
    for i, txt in enumerate (t['REQNUM']):
        plt.annotate(txt, (xdata[i], ydata[i]), size=0.2)
    xdata = np.arange(-5, 35, 0.05)
    ydata = simfitresult[0] + simfitresult[2] * (xdata-args.focus_temp_reference) + simfitresult[3] * (xdata - args.focus_temp_reference)**2 + np.mean(coszd) * simfitresult[1]
    plt.plot(xdata, ydata, '-', label=f"Temperature term {simfitresult[2]:5.4f}[{simfiterrors[2]:5.4f}]T\n + {simfitresult[3]:6.5f}[{simfiterrors[3]:6.5f}]T^2\nTref={args.focus_temp_reference} deg C")

    plt.xlabel('FOCTEMP [deg C]')
    plt.ylabel('FOCUS [mm]')
    plt.xlim([-6, 35])
    set_ylim(t['ACTFOCUS'], focusvaluerange)
    plt.legend()
    plt.title("Temp vs absolute focus position")

    ### ZD
    plt.subplot(5, 2, 2)
    ydata = t['ACTFOCUS']
    plt.plot(coszd, ydata, '.', c='grey')
    plt.plot(coszd[goodvalues], ydata[goodvalues], '.', c='red')
    plt.xlabel('cos (ZD)')
    plt.ylabel('FOCUS [mm], T corrected')
    plt.xlim([0, 1.05])
    set_ylim(t['ACTFOCUS'], focusvaluerange)
    xdata = np.arange(0, 4.1, 0.05)
    deltaTmean = np.mean (temp) -args.focus_temp_reference
    ydata = simfitresult[0] + simfitresult[1] * xdata + deltaTmean * simfitresult[2] + deltaTmean**2 * simfitresult[3]
    plt.plot(xdata, ydata, '-', label=f"Compression term {simfitresult[1]:5.4f} +/- {simfiterrors[1]:5.4f}")
    plt.legend()
    plt.title("Zenith distance vs abs focus position")

    # Temperature Residual after multilinear fit corrections
    plt.subplot(5, 2, 3)
    plt.title("Residual after temperature and ZD correction")
    plt.plot(temp, residual_focus, '.', c='grey')
    plt.plot(temp[goodvalues], residual_focus[goodvalues], '.', c='red')
    plt.xlabel('FOCTEMP [deg C]')
    plt.ylabel('ZD & Temp corrected Focus')
    plt.xlim([-6, 35])
    plt.ylim([-focustermrange /2, focustermrange /2])

    # Compression Residual after multilinear fit corrections
    plt.subplot(5, 2, 4)
    plt.title("Residual after temperature and ZD correction")
    plt.plot(coszd, residual_focus, '.', c='grey')
    plt.plot(coszd[goodvalues], residual_focus[goodvalues], '.', c='red')
    plt.xlabel('cos(ZD)')
    plt.ylabel('ZD & Temp corrected Focus')
    plt.xlim([0, 1.05])
    plt.ylim([-focustermrange /2, focustermrange/2 ])

    plt.subplot(5, 1, 3)
    # Time line plot - history of autofocus corrections
    plt.plot (t['DATE-OBS'], t['ACTFOCUS'], '.', c='grey')
    plt.plot (t['DATE-OBS'][goodvalues], t['ACTFOCUS'][goodvalues], '.', c='red', label="Actual Focus")
    plt.plot (t['DATE-OBS'], t['FOCAFOFF'], '.', c='grey', label="AUTOFOCUS Correction")
    plt.plot (t['DATE-OBS'][goodvalues], t['FOCAFOFF'][goodvalues], '.', c='blue', label="AUTOFOCUS Correction")

    plt.xlabel('DATE-OBS')
    plt.ylabel('AUTO FOCUS TERM [mm]')
    plt.title("History of AF corrections")
    plt.setp(plt.gca().xaxis.get_minorticklabels(), rotation=25)
    plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=25)
    plt.legend()
    set_ylim(t['FOCAFOFF'], 3*focusvaluerange)

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
    plt.plot(xdata, ydata, '.', c='grey')
    plt.plot(xdata[goodvalues], ydata[goodvalues], '.', c='red')

    plt.subplot(5, 2, 7)
    xdata = t['FOCTEMP']
    ydata = t['FOCAFOFF']

    plt.xlabel('Temperature')
    plt.ylabel('Auto focus corrections')
    plt.xlim([-6, 35])
    set_ylim(t['FOCAFOFF'], 2*focusvaluerange)

    plt.plot(xdata, ydata, '.', c='grey')
    plt.plot(xdata[goodvalues], ydata[goodvalues], '.', c='red')

    plt.suptitle("{} {} {} \n".format(site, enc, tel), fontsize=24, y=1.05)

    # Lok at some inttersting residuals
    plt.subplot(5, 2, 9)
    xdata = t['REFHUMID']
    ydata = residual_focus

    plt.xlabel('Humidity')
    plt.ylabel('Residual after temperature and ZD correction')
    plt.plot(xdata, ydata, '.')
    set_ylim(residual_focus, focusvaluerange)
    #plt.ylim([-focustermrange / 5, focustermrange / 5])

    plt.subplot(5, 2, 10)
    xdata = t['CCDATEMP']
    ydata = residual_focus

    plt.xlabel('CCD Temp [deg C]')
    plt.ylabel('Residual after temperature and ZD correction')
    plt.plot(xdata, ydata, '.')
    set_ylim(residual_focus, focusvaluerange)
    #plt.ylim([-focustermrange / 5, focustermrange / 5])

    plt.suptitle("{} {} {} \n".format(site, enc, tel), fontsize=24, y=1.05)

    plt.tight_layout(pad=1)
    plt.savefig(f"focusstack_{site}_{enc}_{tel}.pdf")
    plt.clf()
    return t


def getargs():
    parser = argparse.ArgumentParser(
        description='Analyse focus dependency on temperature and stuff')

    parser.add_argument('--site', default='ogg', type=str)
    parser.add_argument('--dome', default='clma', type=str)
    parser.add_argument('--tel', default='0m4c', type=str)
    parser.add_argument('--after', type=dt.datetime.fromisoformat, help="Consider data only after the given date, in ISO format (e.g., 2022-09-01 00:00:00")
    parser.add_argument('--before',  type=dt.datetime.fromisoformat, default=dt.datetime.utcnow().isoformat())
    parser.add_argument('--maxfwhm', type=float, default=5.0, help='Maximum allowable seiing in \'\'')
    parser.add_argument('--maxsunalt', type=float, default=-18, help='Maximum allowable sun altitude')
    parser.add_argument('--fitorder', type=int, default=1, choices=[1,2], help='Fit order')
    parser.add_argument('--focus_temp_reference', type=float, default=None , help='Focus temp reference point')
    parser.add_argument('--loglevel', dest='log_level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARN'], help='Set the debug level')

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    if args.after is None:
        args.after = datetime.datetime.fromisoformat("2023-01-01")


    if (args.fitorder > 1) and (args.tel[:3] not in focus_temp_reference_points.keys() and (args.focus_temp_reference is None)):
        log.error(f"Polynomial fit order > 1 not supported for this telescope class {args.tel[:3]}. Need to specifying different temp reference points with --focus_temp_reference")
        args.focus_temp_reference = 0
        exit(1) 
    elif args.focus_temp_reference is None:
        args.focus_temp_reference = focus_temp_reference_points[args.tel[:3]]
        log.info(f"Using {args.focus_temp_reference} as focus temp reference point for {args.tel[:3]} telescope class")
    else:
        pass

    return args


def main():
    args = getargs()

    print(f"Analysing focus stack data for {args.site} {args.dome} {args.tel}")
    analysecamera(args)


if __name__ == '__main__':
    main()
