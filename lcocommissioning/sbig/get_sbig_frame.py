import argparse
import datetime
import logging
import math
import tempfile
import threading
from ctypes import *
import astropy.io.fits as fits
import requests
import time
from lcocommissioning.common.lco_ccdlab import LED_Illuminator

_logger = logging.getLogger(__name__)

requests_logger = logging.getLogger('connectionpool')
requests_logger.setLevel(logging.ERROR)




class restcam:


    def __init__(self, ipaddr):
        self.ipaddr = ipaddr
        self.bin_x = 1
        self.bin_y = 1

    """ Relase camera and close sdk """
    def close(self):
        pass


    def __str__(self):
        w,h = self.getPixelSize()
        return f" {str(self.id.value.decode('UTF-8'))} Pixels: {w} x {h} pixels"

    def setReadMode (self, readmodenumber):
        pass


    def sendSettings (self, key, value):
        r = requests.get(f"http://{self.ipaddr}/api/ImagerSetSettings.cgi?{key}={value}")
        return r

    def getSettings (self, datum):

        r = requests.get(f"http://{self.ipaddr}/api/ImagerGetSettings.cgi?{datum}")
        _logger.debug (f'get Settings text: {r.text}')
        return r.text

    def setBinning(self, wbin, hbin):
        self.bin_x = wbin
        self.bin_y = hbin
        r = self.sendSettings('BinX', self.bin_x)
        r = self.sendSettings('BinY', self.bin_y)
        _logger.info(f"Binning returned {r} {r.text}")


    def setGain(self, gain):
        """ Set camera gain """
        self.gain = gain

    def getGain (self):
        return 1.0

    def setTemp (self, temperature):
        self.setpoint = float(temperature)
        assert (self.setpoint > -30)
        assert (self.setpoint < 30)
        self.sendSettings('CCDTemperatureSetpoint', self.setpoint)
        self.sendSettings('CoolerState',1)


    def getTemp (self):
        try:
            r = float(self.getSettings('CCDTemperature'))
        except:
            r = math.nan
        return r

    def getTempSetpoint (self):
        try:
            r = float(self.getSettings('CCDTemperatureSetpoint'))
        except:
            r = math.nan
        return r


    def setROI(self,x,y):

        pass



    def getChipDimensions(self):
        return self.chipw.value, self.chiph.value

    def getPixelSize(self):
        return self.w.value, self.h.value

    def getStatus (self):
        try:
            status = requests.get (f"http://{self.ipaddr}/api/ImagerState.cgi").json()
        except:
            logging.error(f"Error while asking for SBIG status : {status}")
        return status

    def getframe(self, exptime, filename, args=None):

        _logger.info (f"Starting exposure {exptime} seconds")

        self.sendSettings('OverScan',1)
        frametype = 1
        if ('b00' in filename) or ('d00') in filename:
            frametype=2

        params = {
            'Duration': exptime,
            'FrameType' : frametype,
        }
        r = requests.get (f"http://{self.ipaddr}/api/ImagerStartExposure.cgi", params = params  )
        _logger.debug(f"Started exposure {r}")
        time.sleep (0.5)

        while ( self.getStatus() > 0):
            time.sleep (0.5)
        start_readout = datetime.datetime.utcnow()
        _logger.info ("Downloading image from camera")
        fitsdata = requests.get (f"http://{self.ipaddr}/api/Imager.FIT", params ={}  )

        _logger.info (f"writing to fits file {filename}")
        with open(filename, "wb") as myfile:
            myfile.write (fitsdata.content)
            myfile.close()

        hdul = fits.open (filename)

        object = None
        if args is not None:
            object = f"led {args.ledvoltage} nburst {args.nburstcycles}"

        prihdr = hdul[0].header
        prihdr['OBJECT'] = object
        prihdr['EXPTIME'] = exptime
        prihdr['FILTER'] = 'None'
        prihdr['AIRMASS'] = 1.0
        prihdr['DATE-OBS'] = start_readout.strftime('%Y-%m-%dT%H:%M:%S')
        #prihdr['GAIN'] = self.gain.value
        prihdr['CCDSUM'] = f"[{self.bin_x} {self.bin_y}]"
        prihdr['READMODE'] = 'default'
        prihdr['CCDTEMP'] = prihdr['CCD-TEMP']

        hdul.writeto (filename,  overwrite=True )
        end_fitswrite = datetime.datetime.utcnow()
        _logger.info("done taking image")




def main():
    args = parseCommandLine()


    if args.testled is not None:
        lab = LED_Illuminator()
        print ("LED ON")
        lab.expose(exptime = args.testled, overhead = 1, block=True, voltage=args.ledvoltage)
        print ("LED OFF")
        exit (0)

    if args.testledburst:
        lab = LED_Illuminator()
        print ("LED ON")
        lab.expose_burst(exptime = args.testledburst, overhead = 1, block=True, voltage=args.ledvoltage)
        print ("LED OFF")
        exit (0)

    qhyccd = restcam("10.6.11.64:8080")


    if args.settemp is not None:
        print (f"Setting target temperature to {args.settemp: 6.2f}")
        qhyccd.setTemp(args.settemp)
        qhyccd.close()
        exit(0)


    if args.gettemp:
        t = qhyccd.getTemp()
        setpoint = qhyccd.getTempSetpoint()
        print (f"detector temperature is {t: 6.2f} deg C, setpoint is {setpoint: 6.2f}")
        qhyccd.close()
        exit(0)

    qhyccd.setBinning (args.binning, args.binning)

    suffix = None
    if args.bias:
        suffix = 'b00'
    if args.dark:
        suffix = 'd00'
    if args.flat:
        suffix = 'f00'
        lab = LED_Illuminator()


    for exptime in args.exptime:
        _logger.info (f"taking exposures for exptime {exptime}")
        for ii in range (args.expcnt):
            imagename=f"{args.outputpath}/restcam-{datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%S')}.{suffix}.fits"
            if args.flat and lab is not None:
                if args.nburstcycles is None:
                    # This is a conventional exposure where we ensure the LED is on befor we open the shutter and stays on until shutter closes.
                    _logger.info ("Starting conventional shutter-defined exposure")
                    lab.expose(exptime = exptime, overhead = 2, block=False, voltage=args.ledvoltage)
                else:
                    # Here we open the shutter, and then turn the LED on for a determined amount of time. it takes a few seconds from requesting an exposure
                    # until the shutter actually opens. Hence we are putting the LED con command into a background thread that starts its working day with sleeping.

                    _logger.info (f"Starting frequencey generator defined exposure for {args.nburstcycles} cycles.")
                    th =threading.Thread ( target=lab.expose_burst, kwargs={'exptime':exptime, 'ncycles':args.nburstcycles, 'overhead':7, 'voltage':args.ledvoltage, 'block':False})
                    th.start()

            qhyccd.getframe(exptime, imagename, args=args)

            if args.flat:
                time.sleep (1)


    qhyccd.close()
    if args.flat:
        lab.close()
    exit(0)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='QHYCCD command line control')

    actions = parser.add_mutually_exclusive_group()
    actions.add_argument ("--bias", action = "store_true", help="Take bias exposures. Set number via --expcnt" )
    actions.add_argument ("--dark", action = "store_true", help="Take dark exposures. Set number via --expcnt and darktime via --exptime" )
    actions.add_argument ("--flat", action = "store_true", help="Take flat exposures. Set number via --expcnt and darktime via --exptime" )
    actions.add_argument ("--settemp", type=float, help="Set CCD target temperature" )
    actions.add_argument ("--gettemp", action = "store_true",  help="get CCD target temperature" )
    actions.add_argument ("--testled", type=float,  help="testled" )
    actions.add_argument ("--testledburst", type=float,  help="testled in burst mode" )

    actions.add_argument ("--chamberpump", type = bool,  help="cycle detector chaber decissitant" )


    parser.add_argument ("--ledvoltage", type=float, default = None)
    parser.add_argument('--expcnt', type=int, dest="expcnt", default=1)
    parser.add_argument('--exptime', type=float, nargs='*', default=[0,])
    parser.add_argument('--gain', type=int, default=5)
    parser.add_argument('--binning', type=int, default=1)
    parser.add_argument('--readmode', type=int, default=0)
    parser.add_argument('--nburstcycles', type=int, default = None)



    parser.add_argument('--outputpath', type=str, default="data", help="outputpath")
    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    return args





if __name__ == '__main__':
    main()
