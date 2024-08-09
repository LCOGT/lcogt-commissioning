import argparse
import datetime
import logging
import threading
import time
from ctypes import *
import astropy.io.fits as fits
import numpy as np
from lcocommissioning.cmostest.qhyccdpython.libqhy import *
from lcocommissioning.common.lco_ccdlab import LED_Illuminator

_logger = logging.getLogger(__name__)

class QHYCCD:

    def __init__(self):
        self.overscan = None
        self.qhyccd = CDLL('/usr/local/lib/libqhyccd.so')
        self.qhyccd.GetQHYCCDParam.restype = c_double
        self.qhyccd.OpenQHYCCD.restype = ctypes.POINTER(c_uint32)
        ret = -1
        self.qhyccd.InitQHYCCDResource()
        self.qhyccd.ScanQHYCCD()
        type_char_array_32 = c_char * 32
        self.id = type_char_array_32()
        self.qhyccd.GetQHYCCDId(c_int(0), self.id)  # open the first camera
        _logger.info(f"Found camera: {self.id.value.decode('UTF-8')}")
        self.cam = self.qhyccd.OpenQHYCCD(self.id)
        self.qhyccd.SetQHYCCDStreamMode(self.cam, 0)  # 0 for single frame
        self.qhyccd.InitQHYCCD(self.cam)
        self.chipw = c_double()
        self.chiph = c_double()
        self.w = c_uint()
        self.h = c_uint()
        self.pixelw = c_double()
        self.pixelh = c_double()
        self.bpp = c_uint()
        self.gain = c_double()
        self.channels = c_uint32(1)
        self.qhyccd.GetQHYCCDChipInfo(self.cam, byref(self.chipw), byref(self.chiph), byref(self.w), byref(self.h),
                                      byref(self.pixelw), byref(self.pixelh), byref(self.bpp))
        self.setBinning(1, 1)
        self.setROI(self.w.value, self.h.value)
        self.gain = c_uint(int(self.getGain()))

        _logger.info(self)

        numreadmodes = c_uint32()
        self.qhyccd.GetQHYCCDNumberOfReadModes(self.cam, byref(numreadmodes))
        print (f"BPP:{self.bpp}")
        print(f"Number of readmodes: {numreadmodes.value}")

        for readmode in range(numreadmodes.value):
            readmodename = (c_char * 64)()
            self.qhyccd.GetQHYCCDReadModeName(self.cam, c_uint(readmode), readmodename)
            _logger.info(f"\treadmode #{readmode} : {readmodename.value.decode('UTF-8')} ")

    """ Relase camera and close sdk """

    def close(self):
        self.qhyccd.CloseQHYCCD(self.cam)
        self.qhyccd.ReleaseQHYCCDResource()

    def __str__(self):
        w, h = self.getPixelSize()
        return f" {str(self.id.value.decode('UTF-8'))} Pixels: {w} x {h} pixels"

    def setReadMode(self, readmodenumber):
        numreadmodes = c_uint32()
        self.qhyccd.GetQHYCCDNumberOfReadModes(self.cam, byref(numreadmodes))
        if readmodenumber < numreadmodes.value:
            ret = self.qhyccd.SetQHYCCDReadMode(self.cam, c_uint32(readmodenumber))
            if ret == ERR.QHYCCD_SUCCESS:
                self.readmode = c_uint32(readmodenumber)
                _logger.debug(f"Set readmode to {readmodenumber}")
            else:
                _logger.error(f"Error setting readmode {readmodenumber}")
        else:
            _logger.error("Cannot set readmode larger than available read mode number")

    def setBinning(self, wbin, hbin):
        self.wbin = c_uint(wbin)
        self.hbin = c_uint(hbin)
        ret = self.qhyccd.SetQHYCCDBinMode(self.cam, self.wbin, self.hbin)
        if ret != ERR.QHYCCD_SUCCESS:
            _logger.error(f"Error while setting binning {ret}")
        else:
            _logger.debug(f"Set binning to {wbin} {hbin}")

    def setGain(self, gain):
        """ Set camera gain """

        _logger.debug(f"current gain is {self.getGain()}")

        ret = self.qhyccd.IsQHYCCDControlAvailable(self.cam, CONTROL_ID.CONTROL_GAIN)
        if ret == ERR.QHYCCD_SUCCESS:
            self.gain = c_short(gain)
            self.qhyccd.SetQHYCCDParam(self.cam, CONTROL_ID.CONTROL_GAIN, self.gain)
            _logger.debug(f"Setting gain to {gain}")
        else:
            _logger.warning("Sett gain is not supported on this camera")

    def getGain(self):
        return self.qhyccd.GetQHYCCDParam(self.cam, CONTROL_ID.CONTROL_GAIN)

    def SetBit(self, bpp):
        """ Set camera depth """
        self.bpp.value = c_double(bpp)
        self.qhyccd.SetQHYCCDParam(self.cam, CONTROL_ID.CONTROL_TRANSFERBIT, self.bpp)

    def setTemp(self, temperature):
        self.setpoint = c_double(temperature)

        settempsupported = self.qhyccd.IsQHYCCDControlAvailable(self.cam, CONTROL_ID.CONTROL_COOLER)
        if settempsupported != ERR.QHYCCD_SUCCESS:
            print("Setting temeprature is not supported on this camera.")
            return

        ret = self.qhyccd.SetQHYCCDParam(self.cam, CONTROL_ID.CONTROL_COOLER, self.setpoint)
        if ret == ERR.QHYCCD_SUCCESS:
            print(f"Set temperature control to {temperature}")
        else:
            print(f"Setting temeprature failed with error {ret}")

    def getTemp(self):
        ret = self.qhyccd.GetQHYCCDParam(self.cam, CONTROL_ID.CONTROL_CURTEMP)
        print(f"CCD temperature is {ret}")
        rh = self.getRH()
        return ret

    def setSensorChamberCycle(self, on):
        onoff = 1 if on else 0
        if self.qhyccd.IsQHYCCDControlAvailable(self.cam,
                                                CONTROL_ID.CONTROL_SensorChamberCycle_PUMP) != ERR.QHYCCD_SUCCESS:
            _logger.info("Chamber pump not available")
        else:
            ret = self.qhyccd.SetQHYCCDParam(self.cam, CONTROL_ID.CONTROL_SensorChamberCycle_PUMP, onoff)
            if ret != ERR.QHYCCD_SUCCESS:
                _logger.error(f"Error while  toggling SensorchamberCyclePUmp {ret}")
            else:
                _logger.info(f"Turned SnesorchamberCyclePump to {onoff.value}")

    def getRH(self):
        ret = self.qhyccd.GetQHYCCDParam(self.cam, CONTROL_ID.CAM_HUMIDITY)
        _logger.info(f"Relative Humidity is {ret}")
        return ret

    def setROI(self, x, y):

        if (x <= self.w.value) & (y <= self.h.value):
            self.roi_w = c_uint(x)
            self.roi_h = c_uint(y)

            ret = self.qhyccd.SetQHYCCDResolution(self.cam, 0, 0, self.roi_w, self.roi_h)
            if ret == ERR.QHYCCD_SUCCESS:
                _logger.debug(f'Setting ROI to {x} x {y}')
            else:
                _logger.error("Error while setting ROI")

    def getChipDimensions(self):
        return self.chipw.value, self.chiph.value

    def getPixelSize(self):
        return self.w.value, self.h.value

    def getframe(self, exptime, filename, args=None):

        self.qhyccd.SetQHYCCDParam(self.cam, CONTROL_ID.CONTROL_EXPOSURE,
                                   c_double((exptime * 1000. * 1000.)))  # unit: us
        _logger.info(f"Starting exposure {exptime} seconds")


        start_expose = datetime.datetime.utcnow()
        ret = self.qhyccd.ExpQHYCCDSingleFrame(self.cam)
        if ret != ERR.QHYCCD_SUCCESS:
            _logger.error(f"Failure while exposing image {ret}")
        _logger.debug(f"Starting readout")

        binned_w = c_uint(int(self.roi_w.value / self.wbin.value))
        binned_h = c_uint(int(self.roi_h.value / self.hbin.value))
        self.imgdata = (ctypes.c_uint16 * binned_w.value * binned_h.value)()

        start_readout = datetime.datetime.utcnow()

        ret = self.qhyccd.GetQHYCCDSingleFrame(self.cam, byref(binned_w), byref(binned_h), byref(self.bpp),
                                               byref(self.channels), self.imgdata)
        if ret != ERR.QHYCCD_SUCCESS:
            _logger.error(f"failure while downloading image data {ret}")
        _logger.debug(f"Starting fits write to {filename}")
        start_fitswrite = datetime.datetime.utcnow()
        x = np.asarray(self.imgdata)
        print (f"mean of image: {np.mean(x)}")
        object = None
        if args is not None:
            object = f"led {args.ledvoltage} nburst {args.nburstcycles}"
        prihdr = fits.Header()
        prihdr['OBJECT'] = object
        prihdr['EXPTIME'] = exptime
        prihdr['FILTER'] = 'None'
        prihdr['AIRMASS'] = 1.0
        dateobs = start_readout + datetime.timedelta(seconds=exptime)
        prihdr['DATE-OBS'] = dateobs.strftime('%Y-%m-%dT%H:%M:%S.%f')
        prihdr['GAIN'] = self.gain.value
        prihdr['CCDSUM'] = f"[{self.wbin.value} {self.hbin.value}]"
        prihdr['READMODE'] = self.readmode.value
        prihdr['CCDTEMP'] = self.getTemp()

        hdu = fits.PrimaryHDU(x, header=prihdr)
        hdul = fits.HDUList([hdu])

        hdul.writeto(filename, overwrite=True)
        end_fitswrite = datetime.datetime.utcnow()
        _logger.info("done taking image")
        _logger.debug(f"Timing statistics:\n"
                      f"\tTotal     : {end_fitswrite - start_expose} s\n" \
                      f"\tExposing  : {start_readout - start_expose} s\n" \
                      f"\tReadout   : {start_fitswrite - start_readout} s\n" \
                      f"\tFITS write: {end_fitswrite - start_fitswrite} s\n")


def main():
    args = parseCommandLine()

    if args.testled:
        lab = LED_Illuminator()
        print("LED ON")
        lab.expose(exptime=args.exptime[0], overhead=1, block=True, voltage=args.ledvoltage)
        print("LED OFF")
        exit(0)
    qhyccd = QHYCCD()

    if args.settemp is not None:
        qhyccd.setTemp(args.settemp)
        qhyccd.close()
        exit(0)

    if args.chamberpump is not None:
        qhyccd.setSensorChamberCycle(args.chamberpump)
        qhyccd.close()
        exit(0)

    if args.gettemp:
        qhyccd.getTemp()
        qhyccd.close()
        exit(0)

    qhyccd.setReadMode(args.readmode)
    qhyccd.setGain(args.gain)
    qhyccd.setBinning(args.binning, args.binning)

    suffix = None
    if args.bias:
        suffix = 'b00'
    if args.dark:
        suffix = 'd00'
    if args.flat:
        suffix = 'f00'
        lab = None if args.noled else LED_Illuminator()

    for exptime in args.exptime:
        _logger.info(f"taking exposures for exptime {exptime}")
        for ii in range(args.expcnt):
            imagename = f"{args.outputpath}/qhytest-{datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%S')}.{suffix}.fits"
            actexptime = exptime
            if args.flat and lab is not None:

                if args.prewarmled and (exptime < 60):
                    _logger.info ("Prewarming LED. Stand by")
                    lab.expose(exptime=60, overhead=0, block=True)

                if args.nburstcycles is None:
                    # This is a conventional exposure where we ensure the LED is on befor we open the shutter and stays on until shutter closes.
                    _logger.info ("Starting conventional shutter-defined exposure")
                    lab.expose(exptime=exptime, overhead=2, block=False)
                    time.sleep (0.25)

                else:
                    # Here we open the shutter, and then turn the LED on for a determined amount of time. it takes a few seconds from requesting an exposure
                    # until the shutter actually opens. Hence we are putting the LED con command into a background thread that starts its working day with sleeping.

                    _logger.info (f"Starting frequency generator defined exposure for {args.nburstcycles} cycles.")
                    th =threading.Thread ( target=lab.expose_burst, kwargs={'exptime':exptime, 'ncycles':args.nburstcycles, 'overhead':7, 'voltage':args.ledvoltage, 'block':False})
                    actexptime= exptime+5.5
                    th.start()
                    time.sleep (0.25)

            qhyccd.getframe(actexptime, imagename, args=args)

    qhyccd.close()
    if args.flat:
        lab.close()
    exit(0)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='QHYCCD command line control')

    actions = parser.add_mutually_exclusive_group()
    actions.add_argument("--bias", action="store_true", help="Take bias exposures. Set number via --expcnt")
    actions.add_argument("--dark", action="store_true",
                         help="Take dark exposures. Set number via --expcnt and darktime via --exptime")
    actions.add_argument("--flat", action="store_true",
                         help="Take flat exposures. Set number via --expcnt and darktime via --exptime")
    actions.add_argument("--settemp", type=float, help="Set CCD target temperature")
    actions.add_argument("--gettemp", action="store_true", help="get CCD target temperature")
    actions.add_argument("--testled", action="store_true", help="testled")

    actions.add_argument("--chamberpump", type=bool, help="cycle detector chaber decissitant")


    parser.add_argument("--prewarmled", action="store_true", help="testled")
    parser.add_argument('--expcnt', type=int, dest="expcnt", default=1)
    parser.add_argument('--exptime', type=float, nargs='*', default=[0, ])
    parser.add_argument('--gain', type=int, default=5)
    parser.add_argument('--binning', type=int, default=1)
    parser.add_argument('--readmode', type=int, default=0)
    parser.add_argument('--ledvoltage', type=float, default=5.0)
    parser.add_argument('--nburstcycles', type=int, default=None)
    actions.add_argument("--noled", action="store_true", help="do not use the lab led for illumination")


    parser.add_argument('--outputpath', type=str, default="data", help="outputpath")
    parser.add_argument('--loglevel', dest='log_level', default='DEBUG', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    return args


if __name__ == '__main__':
    main()
