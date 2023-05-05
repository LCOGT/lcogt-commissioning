import argparse
import logging

from astropy.io import fits
from astropy import stats
from astropy import time

log = logging.getLogger(__name__)

def AnalyseFitsObject(FitsObject,where, width=None, ext=0):

    dimx = FitsObject[ext].header['NAXIS1']
    dimy = FitsObject[ext].header['NAXIS2']
    dateobs = FitsObject[ext].header['DATE-OBS']
    timeobj = time.Time (dateobs, format='isot', scale='utc')
    datetimeobj = timeobj.to_datetime()
    if width is not None:
        assert (width < dimx)
        assert (width < dimy)
    else:
        width = int (min(dimx, dimy)/8)

    centerx = int (dimx // 2)
    centery = int (dimx // 2)

    minx = int(centerx - width//2)
    maxx = int(centerx + width//2)
    miny = int(centery - width//2)
    maxy = int(centery + width//2)

    if 'top' in where:
        maxy = dimy-10
        miny = maxy-width
    if 'bottom' in where:
        miny = 11
        maxy = miny + width

    log.debug (f"NAXIS[12] {dimx} {dimy} ROI {minx}<x< {maxx}  {miny}<y< {maxy}")

    roi = FitsObject[ext].data[miny:maxy, minx:maxx]
    mean,median, std = stats.sigma_clipped_stats(roi, sigma=3.0, maxiters=5, axis=None)
    fracsecond = datetimeobj.microsecond/1000000.

    return (fracsecond, mean[0],std)






def AnalyseFitsFile (fname,  where = 'center', width = None):
    fitsobject = fits.open(fname)
    fraseconds, mean, std =  AnalyseFitsObject(fitsobject,where,width)
    fitsobject.close()
    return fraseconds, mean,std


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Extract mean and fractional seconds from a series of FITS input files.')

    parser.add_argument('outputfile', type=str,nargs=1,help="OUtput file")
    parser.add_argument('inputfiles', type=str, nargs="+", help="FITS files to process")
    parser.add_argument('--where', type=str, default='center', choices=['center', 'top', 'bottom'], help="FITS files to process")

    parser.add_argument('--loglevel', dest='log_level', default='DEBUG', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args



args= parseCommandLine()

with open(args.outputfile[0], "a") as outputfile:

    for fitsfile in args.inputfiles:
        fracsec, mean, std= AnalyseFitsFile(fitsfile, where=args.where)

        log.info (f"{fitsfile}  {fracsec} {mean} {std}")
        outputfile.write (f"{fitsfile}  {fracsec} {mean} {std}\n")