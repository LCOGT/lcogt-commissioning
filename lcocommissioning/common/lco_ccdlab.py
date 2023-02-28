import vxi11
import logging
import time

_logger = logging.getLogger(__name__)


class LED_Illuminator:

    def __init__(self, address = '10.6.249.99'):
        self.ins = vxi11.Instrument(address)
        _logger.debug ("VXI11 interface: %s" % (self.ins.ask("*IDN?")))

    def expose(self, exptime, voltage = None, block=True, overhead=0):
        ''' Expose via constant illumination.
        Sets the function generator into pulse mode'''
        _logger.debug ("Lab exposing for % 5.2f s" % (exptime))
        self.ins.write ("puls:mode TRIG")
        self.ins.write (f"burst:ncycles 1")
        self.ins.write ("puls:state ON")
        self.ins.write (f"PULSe:DCYC 100.0" )
        self.ins.write ("PULSe:per %fs" % (exptime+overhead+1))
        self.ins.write ("PULSe:widt %fs" % (exptime+overhead))
        self.ins.write ("PULSe:DELay 0s")
        self.ins.write ("burst:DELay 0s")

        # self.ins.write (f"PULSe:DCYCs 100" )
        if (voltage is not None) and (voltage >=0.) and (voltage <= 5.):
            _logger.debug (f"Setting LED voltage to {voltage}")
            self.ins.write(f"voltage:level:imm:high {voltage} V")

        self.ins.trigger()
        if block:
            _logger.info("Blocking during exposure time")
            time.sleep (exptime)
        _logger.debug ("Done exposing")


    def expose_burst (self, exptime,  frequency=100, ncycles = 10, voltage=None, block=True, overhead=5):
        _logger.info (f"Lab burst exposing for {exptime}, led {voltage}, overhead {overhead}")
        time.sleep(overhead)
        _logger.info ("Done sleeping, firing up the LED")
        if ncycles == 0:
            return
        self.ins.write ("burst:state ON")
        self.ins.write ("burst:mode TRIG")
        self.ins.write (f"burst:ncycles {ncycles}")
        self.ins.write (f"freq:fixed {frequency}Hz")
        self.ins.write (f"PULSe:DCYC 50.0" )
        self.ins.write (f"burst:DELay {overhead}s")
        self.ins.write (f"pulse:DELay {overhead}s")

        if (voltage is not None) and (voltage >=0.) and (voltage <= 5.):
            _logger.debug (f"Setting LED voltage to {voltage}")
            self.ins.write(f"voltage:level:imm:high {voltage} V")

        self.ins.trigger()
        if block:
            _logger.info("Blocking during exposure time")
            time.sleep (exptime)
            _logger.debug ("Done exposing")

    def close(self):
        self.ins.close()
