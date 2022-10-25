
import sys
from time import sleep

import datetime

import archon_v02_plus
import numpy as np
from astropy.io import fits
import logging
import argparse
import astropy.time
import time

log = logging.getLogger(__name__)


import vxi11

class LCOLab:

    def __init__(self):
        self.ins = vxi11.Instrument('172.16.4.6')
        log.info ("VXI11 interface: %s" % (self.ins.ask("*IDN?")))

    def expose(self, exptime, block=True):
        log.info ("Lab exposing for % 5.2f s" % (exptime))
        self.ins.write ("puls:per %f" % (exptime+1))

        self.ins.write ("puls:widt %f" % (exptime))

        self.ins.trigger()
        if block:
            sleep (exptime)
        log.info ("Done exposing")

    def close(self):
        self.ins.close()


import requests
class LCOTelescope:
    def __init__(self, dome='doma', site='bpl'):
        self.dome = dome
        self.site = site
        self.shutterurl = 'http://scicam-fw.1m0a.%s.%s/cgi-bin/expose' % (self.dome, self.site)

    def expose (self, exptime, block = True):
        data = { 'time' : float (exptime), 'light': '1'}
        log.debug (data)
        log.info ("Exposing w/ shutter at %s %s for %6.2f sec" % (self.dome, self.site, exptime))
        r = requests.post(self.shutterurl,data )
        log.debug ("shutter url:       %s\nshutter request answer:  %s " % (r.url, r.text))
        if (block):
            log.debug ("Sleeping for %f seconds" % (exptime))
            sleep (exptime)
            log.debug ("Exposure time over.")

        log.debug ("waiting for shutter to close...")
        sleep (0.3)


    def openShutter(self):
        log.info ("Shutter open")
        self.expose (-2, False)

    def closeShutter(self):
        log.info ("Shutter close")
        self.expose (-1, False)


    def deFocus (self, defocusMMfocalplane):

        pass


class archonexposure:


    def __init__(self, host, configfile,  offline=False):
        self.host = host

        if offline is not True:

            self.archon = archon_v02_plus.Archon()
            self.instrumentstatus = {}

            log.info ("Connecting to Archon at ip " + host)
            self.open()
            if configfile:
                try:
                    log.info ("Loading config file: " + configfile)
                    cfg = self.archon.load_configuration(configfile)
                    self.pon()
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

            if item.strip() == 'MOD8/TEMPB':
                self.instrumentstatus['TEMPB'] = float(status[item])

            if item.strip() == 'MOD8/TEMPA':
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
            log.error ("Could not connect to Archon: " + str(ex))

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


    def flashconfig(self):
        log.info ("Store active config")
        self.archon.command ('FLASHACTIVECONFIG')
        sleep (2)
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
        log.debug ("Deactivate background cleaning")
        self.saveupdateconfigparam('BackgroundClean', 'BackgroundClean=0')
        self.fetchlog()

    def pof (self):
        log.info ("POF")
        self.archon.command('POWEROFF')
        sleep(2)
        self.fetchlog()


    def cleanAndIntegrate (self, numflushes = 3):

        self.callFlush(numflushes)



    def callFlush (self, numflushes = 5):
        self.backgroundCleanOff()
        sleep(0.05)
        for ii in range (numflushes):

            self.saveupdateconfigparam('Flushes', 'Flushes=1')
            log.info('waiting for flush %d to complete...' % (ii))
            sleep(2.8)  # anything better than just trying to wait for the next flush iteration to complete?

        log.info ("Flush complete")



    def callSingleFlush(self):
        self.backgroundCleanOff()
        self.saveupdateconfigparam('SingleFlush', 'SingleFlush=1')
        sleep(0.05)
        log.debug ("Single flush complete")


    def contniuousreadout (self, extraheader = None, filename = None, exptimems=0, nframes = 100):
        """ Read out cotninuously and add up to a data cube"""



        lastframe, buf, framew, frameh, samplemode, lasttimestamp = self.archon.newest();
        frame = lastframe
        # Start readout loop
        log.info ("Set video exposure time to %d ms" % exptimems)
        self.saveupdateconfigparam('IntegrationTimeMS', 'IntegrationTimeMS=%d' %exptimems)

        log.info ("Starting readout loop")
        self.saveupdateconfigparam('Readouts', 'Readouts=1')
        self.saveupdateconfigparam('IntegrationTimeMS', 'IntegrationTimeMS=%d' %exptimems)

        datacube= None
        deltats = []
        try:
            # wait for next frame to be ready
            while nframes > 0:

                startwait = time.time()
                while True:
                    # wait for next image to be ready from archon
                    frame, buf, framew, frameh, samplemode, timestamp = self.archon.newest()
                    if (frame != lastframe):
                        lastframe = frame
                        deltat = timestamp - lasttimestamp
                        lasttimestamp = timestamp
                        deltats.append (deltat)
                        break

                framesize = framew * frameh * 2 # 2 byte per pixels
                datatype = np.dtype('<u2')    # 16-bit little endian
                bpp = 16


                linesize = self.archon.BURST_LEN
                lines = (framesize + linesize - 1) // linesize
                log.info('Frame is {:d}x{:d}, {:d}bpp {:d}'.format(framew, frameh, bpp, timestamp))

                self.archon.command('LOCK{:d}'.format(buf+1 ))
                msgref = self.archon.msgref
                self.archon.send('FETCH{:08X}{:08X}'.format(((buf+1) | 4) << 29, lines))

                rawdata = bytearray()
                bytesremaining = framesize
                for i in range(lines):
                    self.archon.msgref = msgref
                    rawdata.extend(self.archon.binrecv()[0:min(linesize, bytesremaining)])
                    bytesremaining -= linesize

                self.archon.msgref = (self.archon.msgref + 1) % 256
                self.archon.command('LOCK0')

                pixels = np.frombuffer(rawdata, np.uint8).view(datatype)
                image = np.asarray (np.reshape(pixels, (frameh, framew)))

                if datacube is None:
                    datacube = np.zeros((1,image.shape[0],image.shape[1]))
                    datacube[0,:,:] = image[:,:]

                else:
                    datacube = np.concatenate ( (datacube, np.array([image])))

                nframes = nframes -1

        except:
            log.exception()

        finally:
            # end readout loop
            log.info ("Ending readout loop")
            self.saveupdateconfigparam('Readouts', 'Readouts=0')
            if filename is not None:
                hdu = fits.PrimaryHDU(datacube)
                hdul = fits.HDUList([hdu])
                log.info('writing image to file %s' % (filename))
                hdul.writeto(filename, overwrite=True)
            deltats = np.asarray (deltats) / 100. / 1000.
            log.info ("Time stamps: %f mean pm %f ms " % (np.mean (deltats[1:]), np.std (deltats[1:])))


    def readoutSingleFrame(self, extraheader=None, filename=None):
        lastframe, buf, framew, frameh, samplemode, timestamp = self.archon.newest()
        log.info('waiting for next frame to be ready; current frame # is %d, buffer %d' % (lastframe, buf))

        ### Readout one frame from detctor
        self.saveupdateconfigparam('Readouts', 'Readouts=1')

        waittime = 0
        while True:
            frame, buf, framew, frameh, samplemode, timestamp = self.archon.newest()
            if (frame != lastframe):
                break
            sleep(0.1)
            waittime = waittime + 0.1
        log.debug("Image readout complete (frame number %d in buffer %d), waited % 5.1f seconds" % (frame, buf, waittime))
        # download it
        framesize = framew * frameh
        if samplemode:  # if samplemode, it's 32-bit data, otherwise it's 16-bit data
            framesize *= 4
            datatype = np.dtype('<u4')  # 32-bit little endian
            bpp = 32
        else:
            framesize *= 2
            datatype = np.dtype('<u2')  # 16-bit little endian
            bpp = 16

        linesize = self.archon.BURST_LEN
        lines = (framesize + linesize - 1) // linesize
        log.info('Frame is {:d}x{:d}, {:d}bpp'.format(framew, frameh, bpp))

        self.archon.command('LOCK{:d}'.format(buf + 1))
        msgref = self.archon.msgref
        self.archon.send('FETCH{:08X}{:08X}'.format(((buf + 1) | 4) << 29, lines))
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

        if extraheader is not None:
            for item in extraheader:
                hdul[0].header[item] = extraheader[item]

        if not filename:
            filename = 'test{:05d}x{:05d}_{:d}.fits'.format(framew, frameh, frame)

        print('writing image to file %s' % (filename))
        hdul.writeto(filename)



    def loadconfig (self):
        pass

    def setPreExpHeader (self):
        now =  astropy.time.Time.now()
        now.format = 'isot'
        self.instrumentstatus['DATE-OBS'] = now.value

    def bias (self, archon):
        archon.cleanAndIntegrate()
        self.setPreExpHeader()
        self.instrumentstatus['EXPTIME'] = float(0)
        self.instrumentstatus['OBSTYPE'] = 'BIAS'

        log.debug ("Bias readout")
        archon.readoutSingleFrame(filename=f"archon-{datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%S')}.b00.fits")


    def dark (self, darktime=0, archon):
        archon.cleanAndIntegrate()
        self.setPreExpHeader()
        self.instrumentstatus['EXPTIME'] = float(darktime)
        self.instrumentstatus['OBSTYPE'] = 'DARK'

        log.info("Dark exposure % 5.2f seconds" % (darktime))
        sleep (darktime)

        archon.readoutSingleFrame(filename=f"archon-{datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%S')}.d00.fits")

    def flat(self, archon, exptime=1):

        print ("flat field exposure %f seconds" % (exptime))



        ins = LCOTelescope()
        self.setPreExpHeader()
        self.instrumentstatus['EXPTIME'] = float(exptime)
        self.instrumentstatus['OBSTYPE'] = 'FLAT'

        archon.cleanAndIntegrate()
        ins.expose(exptime)


        archon.readoutSingleFrame(filename=f"archon-{datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%S')}.f00.fits")
        log.info ("flat field read out")

    def focussequence(self, archon, focusstep=0.5, nSteps=5, expTime = 10):
        first = True
        ins = LCOTelescope()
        archon.cleanAndIntegrate()

        for ii in range ( int(-nSteps /2),int(nSteps/2)+1,1):
            defocus = focusstep * ii
            ins.defocus (defocus)
            print ("focus step %d  -> % 5.2f " % (ii, defocus))

            ins.expose(expTime)

            self.callSingleFlush()
            if first is True:
                first = False
                self.callSingleFlush()

        ins.deFocus(0)
        archon.readoutSingleFrame()
        ins.close()


        pass

    def ptcSequence (self, archon, texpSat=96, pixelsPerShift=32, nshifts=4, pixelsPerDetector=4096):


        nSteps = int (pixelsPerDetector / pixelsPerShift/ nshifts)
        tExp = texpSat / nSteps

        log.info ("PTC: doing  %i steps, exposed for %5.2f seconds each" % (nSteps, tExp))

        ins = LCOTelescope()
        self.cleanAndIntegrate()


        for ii in range(0, nSteps):
            log.info("  PTC iteration %i of %i" %(ii, nSteps))
            ins.expose(tExp)
            for jj in range (0,nshifts):
                self.callSingleFlush()

        config={}
        config['PTCSTEPS'] = nSteps
        config['PTCEXPT'] = tExp

        archon.readoutSingleFrame(config)


def parseCommandLine( ):
    """ Read command line parameters
    """


    parser = argparse.ArgumentParser(
        description='Fetch bias, dark, exposure from Archon with clean operation')
    parser.add_argument('--log_level', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')
    parser.add_argument('--archonip', default='10.0.0.2',
                        help='IP for Archon controller')
    parser.add_argument('--configfile', dest='configfile',  default =  None)
    parser.add_argument('--loadconfig', action='store_true')
    parser.add_argument('--exptime', default=1, type=float, help="flat or dark exposure time")
    parser.add_argument('--nexp', type=int, default = 1)
    parser.add_argument('--exptimems', type=int, default=0)
    parser.add_argument('--nvideo', type=int, default=20)

    parser.add_argument('--site',   choices=['bpl', 'lsc', 'elp', 'cpt', 'coj',], default=None)
    parser.add_argument('--dome',   choices=['doma', 'domb', 'domc'], default=None)


    imagetype = parser.add_mutually_exclusive_group()

    imagetype.add_argument('--bias', action = "store_true")
    imagetype.add_argument('--dark',  action = "store_true")
    imagetype.add_argument('--flat',  action = "store_true")
    imagetype.add_argument('--continuous', action = "store_true")
    imagetype.add_argument('--pof', action = "store_true")
    imagetype.add_argument('--pon', action = "store_true")
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    return args




def main():
    args = parseCommandLine()


    telescope = None
    if (args.site is not None) and (args.dome is not None):
        log.info ("Setting up telescope")
        telescope = LCOTelescope(dome=args.dome, site=args.site)



    archon = archonexposure(args.archonip, args.configfile, )
    archon.readinstrumentStatus()

    archon.callFlush(1)

    if args.continuous:

        if telescope is not None:
            telescope.openShutter()

        archon.contniuousreadout(None, "continuous.fits", exptimems=args.exptimems, nframes=args.nvideo)

        if telescope is not None:
            telescope.closeShutter()

    if args.flat:
        for ii in range(args.nexp):
            archon.flat(archon, exptime=args.exptime)
    if args.bias:
        for ii in range(args.nexp):
            archon.bias(archon, )
    if args.dark:
        for ii in range(args.nexp):
            archon.dark (archon, args.exptime)

    if args.pon:
        archon.pon()
    if args.pof:
        archon.pof()

    archon.close()





if __name__ == '__main__':
    main()