import sys
from time import sleep
import lcocommissioning.archon.archon_v02_plus
import numpy as np
from astropy.io import fits
import logging
import argparse
import astropy.time
import requests

log = logging.getLogger(__name__)




class LCOTelescope:
    def __init__(self, dome='doma', site='bpl'):
        self.dome = dome
        self.site = site
        self.shutterurl = 'http://fws.1m0a.%s.%s/cgi-bin/expose' % (self.dome, self.site)

    def expose (self, exptime, block = True):
        data = { 'time' : float (exptime), 'light': '1'}
        log.debug (data)
        log.info ("Exposing w/ shutter at %s %s for %6.2f sec" % (self.dome, self.site, exptime))
        r = requests.post(self.shutterurl,data )
        if (block):
            log.debug ("Sleeping for %f seconds" % (exptime))
            sleep (exptime)
            log.debug ("Exposure time over.")


class archonexposure:


    def __init__(self, host, configfile):
        self.host = host
        self.archon = archon_v02_plus.Archon()
        self.instrumentstatus = {}

        log.info ("Connecting to Archon at ip " + host)
        self.open()
        if configfile:
            try:
                log.info ("Loading config file: " + configfile)
                cfg = self.archon.load_configuration(configfile)
            except archon_v02_plus.CommandFailure as ex:
                log.error (ex.message)
                logmsg = self.archon.command('FETCHLOG').decode()
                while logmsg != '(null)':
                    log.error(logmsg)
                    logmsg = self.archon.command('FETCHLOG').decode()
                sys.exit()

        self.readinstrumentStatus()


    def readinstrumentStatus (self):
        status = self.archon.get_status()
        for item in status:

            if item.strip() == 'MOD10/TEMPB':
                self.instrumentstatus['TEMPB'] = float(status[item])

            if item.strip() == 'MOD10/TEMPA':
                self.instrumentstatus['TEMPA'] = float(status[item])

            #print ("% 12s -> %s" % (item, status[item]))


    def fetchlog (self):
        logmsg = self.archon.command('FETCHLOG').decode()
        while logmsg != '(null)':
            log.error("archon Fetchlog: " + logmsg)
            logmsg = self.archon.command('FETCHLOG').decode()


    def open (self ):
        try:
            log.debug ("Attempt to open Archon")
            self.archon.open(ip=self.host)
        except Exception as ex:
            log.error ("Could not connect to Archon: " + ex.message)

    def close (self):
        log.debug ("Closing connection to archon")
        self.archon.close()


    def pon (self):
        """
        Power on the detector
        :return:
        """
        log.info ("PON")
        self.archon.command('POWERON')
        sleep(2)
        self.fetchlog()


    def savereadconfigparam (self, paramname):
        log.warning("Readback is disfunctional, i.e., backplane memory is not yet updated from state machine. Readback can be invalid. ")
        nexp_p = list(self.archon.find_parameter_in_configuration(self.archon.cfg, paramname))
        result = self.archon.read_configuration_line(nexp_p[0])
        return result

    def saveupdateconfigparam (self, paramname, paramline):
        nexp_p = list(self.archon.find_parameter_in_configuration(self.archon.cfg, paramname))
        nexp_p[2] = paramline
        log.debug ("SaveUpdateConfig: %s" % (nexp_p) )
        self.archon.set_configuration_line(*nexp_p)
        self.fetchlog()
        self.archon.load_parameter(paramname)
        self.fetchlog()


    def backgroundCleanOn  (self):
        log.info ("Activate background cleaning")
        self.saveupdateconfigparam('BackgroundClean', 'BackgroundClean=1')
        self.fetchlog()
        self.saveupdateconfigparam('Flushes', 'Flushes=0')
        self.fetchlog()

    def backgroundCleanOff (self):
        log.info ("Deactivate background cleaning")
        self.saveupdateconfigparam('BackgroundClean', 'BackgroundClean=0')
        self.fetchlog()

    def pof (self):
        log.info ("POF")
        self.archon.command('POWEROFF')
        sleep(2)
        self.fetchlog()


    def cleanAndIntegrate (self, numflushes = 1):

        self.callFlush(numflushes)



    def callFlush (self, numflushes = 5):
        self.backgroundCleanOff()
        sleep(0.2)
        for ii in range (numflushes):

            self.saveupdateconfigparam('Flushes', 'Flushes=1')
            log.info('waiting for flush %d to complete...' % (ii))
            sleep(1.2*4)  # anything better than just trying to wait for the next flush iteration to complete?

        log.info ("Flush complete")






    def readoutSingleFrame (self):

        lastframe, buf, framew, frameh, samplemode = self.archon.newest();
        log.info('waiting for next frame to be ready; current frame # is %d, buffer %d' % (lastframe, buf))

        ### Readout one frame from detctor
        self.saveupdateconfigparam('Exposures', 'Exposures=1')

        waittime = 0
        while True:
            frame, buf, framew, frameh, samplemode = self.archon.newest();
            if (frame != lastframe):

                break
            sleep(0.1)
            waittime = waittime + 0.1
        log.debug ("Image readout complete (frame number %d in buffer %d), waited % 5.1f seconds" % (frame, buf, waittime))
        # download it
        framesize = framew * frameh
        if samplemode:  # if samplemode, it's 32-bit data, otherwise it's 16-bit data
            framesize *= 4
            datatype = np.dtype('<u4')    # 32-bit little endian
            bpp = 32
        else:
            framesize *= 2
            datatype = np.dtype('<u2')    # 16-bit little endian
            bpp = 16

        linesize = self.archon.BURST_LEN
        lines = (framesize + linesize - 1) // linesize
        log.info('Frame is {:d}x{:d}, {:d}bpp'.format(framew, frameh, bpp))

        self.archon.command('LOCK{:d}'.format(buf+1 ))
        msgref = self.archon.msgref
        self.archon.send('FETCH{:08X}{:08X}'.format(((buf+1) | 4) << 29, lines))
        log.debug('downloading')
        rawdata = bytearray()
        bytesremaining = framesize
        for i in range(lines):
            self.archon.msgref = msgref
            rawdata.extend(self.archon.binrecv()[0:min(linesize, bytesremaining)])
            bytesremaining -= linesize

        self.archon.msgref = (self.archon.msgref + 1) % 256
        self.archon.command('LOCK0')

        pixels = np.frombuffer(rawdata, np.uint8).view(datatype)
        image = np.reshape(pixels, (frameh, framew))

        hdu = fits.PrimaryHDU(image)
        hdul = fits.HDUList([hdu])

        self.readinstrumentStatus()
        for item in self.instrumentstatus:
            hdul[0].header[item] = self.instrumentstatus[item]

        filename = 'test{:05d}x{:05d}_{:d}.fits'.format(framew, frameh, frame)
        print('writing image to file %s' % (filename))
        hdul.writeto(filename)

        print('image write complete.')


    def loadconfig (self):
        pass

    def setPreExpHeader (self):
        now =  astropy.time.Time.now()
        now.format = 'isot'
        self.instrumentstatus['DATE-OBS'] = now.value

    def bias (self):
        archon.cleanAndIntegrate()
        self.setPreExpHeader()
        self.instrumentstatus['EXPTIME'] = float(0)
        self.instrumentstatus['OBSTYPE'] = 'BIAS'

        log.debug ("Bias readout")
        archon.readoutSingleFrame()


    def dark (self, darktime=0):
        archon.cleanAndIntegrate()
        self.setPreExpHeader()
        self.instrumentstatus['EXPTIME'] = float(darktime)
        self.instrumentstatus['OBSTYPE'] = 'DARK'

        log.info("Dark exposure % 5.2f seconds" % (darktime))
        sleep (darktime)

        archon.readoutSingleFrame(args)

    def flat(self, exptime=1):

        print ("flat field exposure %f seconds" % (exptime))



        ins = LCOTelescope()
        self.setPreExpHeader()
        self.instrumentstatus['EXPTIME'] = float(exptime)
        self.instrumentstatus['OBSTYPE'] = 'FLAT'

        archon.cleanAndIntegrate()
        ins.expose(exptime)


        archon.readoutSingleFrame()
        log.info ("flat field read out")

def parseCommandLine( ):
    """ Read command line parameters
    """


    parser = argparse.ArgumentParser(
        description='Fetch bias, dark, exposure from Archon with clean operation')
    parser.add_argument('--log_level', dest='log_level', default='DEBUG', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')
    parser.add_argument('--archonip', default='127.0.0.1',
                        help='IP for Archon controller')
    parser.add_argument('--loadconfig', dest='configfile', nargs='+', default =  "lco2_500kHz_DRH.acf")

    parser.add_argument('--exptime', default=0, )


    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    return args



def main():

    args = parseCommandLine()

    # telescope = LCOTelescope()
    # telescope.expose(1)
    # exit (0)



    archon = archonexposure(args.archonip, args.configfile)
    archon.readinstrumentStatus()

    archon.pon()
    #archon.callFlush(10)
    

    archon.bias()
    archon.bias()

    for ii in range(2):

        archon.flat(10)

    for ii in range(0):

        archon.bias(2)
        archon.backgroundCleanOn()


    # for int in range (5):
    #     archon.dark (60)
    #


    archon.pof()
    archon.close()



if __name__ == '__main__':
    main()