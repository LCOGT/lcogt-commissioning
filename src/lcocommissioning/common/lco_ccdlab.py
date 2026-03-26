import vxi11
import logging
import time

_logger = logging.getLogger(__name__)


class LED_Illuminator:

    def __init__(self, address = 'functiongenerator.wtf.lco.gtn'):
        self.ins = vxi11.Instrument(address)
        _logger.debug ("VXI11 interface: %s" % (self.ins.ask("*IDN?")))

    def expose(self, exptime, voltage = None, block=True, overhead=0):
        ''' Expose via constant illumination.
        Sets the function generator into pulse mode'''
        if exptime > 950:
            _logger.fatal ("cannot expose for more than 950 seconds. not supported")
            exit(1)
        _logger.debug ("Lab exposing for % 5.2f s" % (exptime))

        self.ins.write ("*RST")
        self.ins.write ("SOURce1:FUNCtion:Shape PULSe")
        self.ins.write ("TRIGger:Source Ext")
        self.ins.write ("VOLTage:LEVel:IMMediate:LOW 0V")
        self.ins.write ("SOURce1:burst:state ON")
        self.ins.write ("SOURce1:puls:mode TRIG")
        self.ins.write (f"SOURce1:burst:ncycles 1")
        self.ins.write ("SOURce1:PULSe:per %fs" % (exptime+overhead+1))
        self.ins.write ("SOURce1:PULSe:widt %fs" % (exptime+overhead))
        self.ins.write ("SOURce1:PULSe:DELay 0s")
        self.ins.write ("SOURce1:burst:DELay 0s")

        # self.ins.write (f"PULSe:DCYCs 100" )
        if not ((voltage is not None) and (voltage >=0.) and (voltage <= 5.)):
            voltage = 5.
        _logger.debug (f"Setting LED voltage to {voltage}")
        self.ins.write(f"voltage:level:imm:high {voltage} V")
        self.ins.write ("OUTPut1:STATe ON")

        self.ins.trigger()
        if block:
            _logger.info("Blocking during exposure time")
            time.sleep (exptime)
        _logger.debug ("Done setting up exposing. ")


    def expose_burst (self, exptime,  frequency=100, ncycles = 10, voltage=None, block=True, overhead=5):
        _logger.info (f"Lab burst exposing for {exptime}, led {voltage}, overhead {overhead}")
        time.sleep(overhead)
        _logger.info ("Done sleeping, firing up the LED")
        if ncycles == 0:
            return

        self.ins.write ("*RST")
        self.ins.write ("SOURce1:FUNCtion:Shape PULSe")
        self.ins.write ("SOURce1:VOLTage:LEVel:IMMediate:LOW 0V")
        self.ins.write ("TRIGger:Source Ext")


        self.ins.write ("SOURCE1:burst:mode TRIG")
        self.ins.write ("SOURCE1:burst:state ON")
        self.ins.write (f"SOURCE1:freq:fixed {frequency}Hz")
        self.ins.write (f"SOURCE1:burst:ncycles {ncycles}")
        self.ins.write (f"SOURCE1:PULSe:DCYC 50.0" )
        self.ins.write (f"SOURCE1:burst:DELay {overhead}s")
        self.ins.write (f"SOURCE1:pulse:DELay {overhead}s")

        # self.ins.write (f"PULSe:DCYCs 100" )
        if not ((voltage is not None) and (voltage >=0.) and (voltage <= 5.)):
            voltage = 5.
        _logger.debug (f"Setting LED voltage to {voltage}")
        self.ins.write(f"voltage:level:imm:high {voltage} V")
        self.ins.write ("OUTPut1:STATe ON")

        self.ins.trigger()
        if block:
            _logger.info("Blocking during exposure time")
            time.sleep (exptime)
            _logger.debug ("Done exposing")

    def close(self):
        self.ins.close()
