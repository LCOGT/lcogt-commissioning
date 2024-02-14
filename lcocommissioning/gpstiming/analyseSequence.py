import argparse
import logging
import numpy as np
from astropy.io import fits
from astropy import stats
from astropy import time
from astropy.nddata import block_reduce

log = logging.getLogger(__name__)

def AnalyseFitsObject(FitsObject,where, width=None, ext='SCI'):

    dimx = FitsObject[ext].header['NAXIS1']
    dimy = FitsObject[ext].header['NAXIS2']
    dateobs = FitsObject[ext].header['DATE-OBS']
    exptime = FitsObject[ext].header['EXPTIME']
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

    return (fracsecond, mean, std)



def AnalyseFitsFile (fname,  where = 'center', width = None):
    fitsobject = fits.open(fname)
    fraseconds, mean, std =  AnalyseFitsObject(fitsobject,where,width)
    fitsobject.close()
    return fraseconds, mean,std


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Extract mean and fractional seconds from a series of FITS input files.')

    parser.add_argument('--out', type=str, help="Output file")
    parser.add_argument('inputfiles', type=str, nargs="+", help="FITS files to process")

    parser.add_argument('--where', type=str, default='center', choices=['center', 'top', 'bottom', 'block'], help="FITS files to process")
    parser.add_argument('--width', type=int, default=16)
    parser.add_argument('--loglevel', dest='log_level', default='DEBUG', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args



def SimpleregionAnalysis(args):
    with open(args.out, "a") as outputfile:

        for fitsfile in args.inputfiles:
            fracsec, mean, std= AnalyseFitsFile(fitsfile, where=args.where, width=args.width)

            log.info (f"{fitsfile}  {fracsec} {mean} {std}")
            outputfile.write (f"{fitsfile}  {fracsec} {mean} {std}\n")



def BlockAveAnalysis(args, ext='SCI'):
    fitstemplate = fits.open (args.inputfiles[0])
    dimX = fitstemplate[ext].header['NAXIS1']
    dimY = fitstemplate[ext].header['NAXIS2']
    fitstemplate.close()

    meanarray = np.zeros (( len(args.inputfiles), dimY// args.width, dimX // args.width))
    stdarray = np.zeros (( len(args.inputfiles), dimY// args.width, dimX // args.width))
    fracsec = np.zeros (len(args.inputfiles))


    idx = 0
    for fitsfile in args.inputfiles:
        log.debug (fitsfile)
        f = fits.open (fitsfile)
        meanarray[idx, :,:] = block_reduce(f[ext].data, args.width, func=np.mean)
        stdarray[idx, :,:] = block_reduce(f[ext].data, args.width, func=np.std)
        dateobs = f[ext].header['DATE-OBS']
        timeobj = time.Time (dateobs, format='isot', scale='utc')
        fracsec[idx] = timeobj.to_datetime().microsecond/1000000.
        f.close()
        idx = idx+1

    outf = fits.PrimaryHDU(meanarray)
    outf.header['BLK'] = args.width
    avg = fits.ImageHDU(stdarray)
    hdu = fits.BinTableHDU.from_columns([fits.Column(name='fracsec', array=fracsec, format='E'),])

    hdul = fits.HDUList([outf, avg, hdu])
    hdul.writeto(args.out, overwrite=True)

args= parseCommandLine()
if 'block' not in args.where:
    SimpleregionAnalysis(args)
else:
    BlockAveAnalysis(args)





