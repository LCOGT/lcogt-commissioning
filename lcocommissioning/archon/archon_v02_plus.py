# Refactored for python3

import socket
# import configparser
import select
from time import sleep

class CommandFailure(Exception):
    """Exception raised when an Archon command fails.

    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

class Archon(object):
    def __init__(self, update_rate=2.0, status_filename=None):
        self.BURST_LEN = 1024
        # Message reference
        self.msgref = 0
        self.msgbuf = b''
        self.status_filename = status_filename
        self.update_rate = update_rate
        
    def open(self, ip='10.0.0.2', port=4242):
        self.controller = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.controller.connect((ip, port))
        self.controller.settimeout(10)

    def close(self):
        self.controller.close()

    def send(self, cmd):
        self.controller.sendall(str.encode('>%02X%s\n' % (self.msgref, cmd)))
        self.msgref = (self.msgref + 1) % 256
        return

    def recv(self):
        while not (b'\n' in self.msgbuf):
            self.msgbuf = self.msgbuf + self.controller.recv(4096)
        (reply, self.msgbuf) = self.msgbuf.split(b'\n', 1)
        if reply[0:3].decode() != '<%02X' % self.msgref:
            raise CommandFailure('Archon recv failed')
        self.msgref = (self.msgref + 1) % 256
        return reply[3:]
    
    def binrecv(self):
        binlen = self.BURST_LEN + 4
        while len(self.msgbuf) < binlen:
            self.msgbuf = self.msgbuf + self.controller.recv(4096)
        reply = self.msgbuf[0:binlen]
        self.msgbuf = self.msgbuf[binlen:]
        if reply[0:4].decode() != '<%02X:' % self.msgref:
            raise CommandFailure('Archon binrecv failed')
        self.msgref = (self.msgref + 1) % 256
        return reply[4:]
    
    def command(self, c):
        cmd = str.encode('>%02X%s\n' % (self.msgref, c))
        #print(cmd)
        self.controller.sendall(cmd)
        reply = b'';
        while not (b'\n' in reply):
            if select.select([self.controller], [], [], 0.01)[0]:
                reply = reply + self.controller.recv(1)
        reply = reply.splitlines()[0]
        #print(reply)
        if reply[0:3].decode() != '<%02X' % self.msgref:
            raise CommandFailure("Archon command '{:s}' failed".format(c))
        self.msgref = (self.msgref + 1) % 256
        return reply[3:]

    def read_configuration_line(self, line_number):
        return self.command('RCONFIG{:04X}'.format(line_number)).decode()
    
    def set_configuration_line(self, line_number, k, v):
        self.command('WCONFIG{:04X}{:s}={:s}'.format(line_number, k, v))
        s = self.read_configuration_line(line_number)
        return s.split('=',1)

    def get_status(self):
        status = {}
        for pair in self.command('STATUS').split():
            d = pair.decode().split('=')
            status[d[0]] = d[1]
        return status

    def append_status(self, filename):
        txt = self.command('STATUS').decode()
        with open(filename, 'a') as fh:
            print(txt, file=fh)
            fh.close()

    def newest(self):
        framestatus = {}
        for pair in self.command('FRAME').split():
            d = pair.decode().split('=')
            framestatus[d[0]] = d[1]
        rbuf = int(framestatus['RBUF']) - 1
        frames = []
        framecomplete = []
        for i in range(1,4):
            frames.append(int(framestatus['BUF{:d}FRAME'.format(i)]))
            framecomplete.append(int(framestatus['BUF{:d}COMPLETE'.format(i)]) == 1)
        if rbuf >= 0 and rbuf <= 2:
            newestframe = frames[rbuf]
            newestbuf = rbuf
        else:
            newestframe = -1
            newestbuf = 0
        for i in range(0, 3):
            if frames[i] > newestframe and framecomplete[i]:
                newestframe = frames[i]
                newestbuf = i
        framew = int(framestatus['BUF{:d}WIDTH'.format(newestbuf + 1)])
        frameh = int(framestatus['BUF{:d}HEIGHT'.format(newestbuf + 1)])
        samplemode = int(framestatus['BUF{:d}SAMPLE'.format(newestbuf + 1)])
        return (newestframe, newestbuf, framew, frameh, samplemode)

    def load_configuration(self, filename):
        lines = []
        print('reading configuration from {:s}...'.format(filename))
        with open(filename, 'r') as f:
            line = f.readline()
            while not line.startswith('[CONFIG]') and line != '':
                line = f.readline()
            if line == '':
                print('no [CONFIG] section found!')
                return None
            line = f.readline()
            while '=' in line:
                kvpair = line.split('=', 1)
                kvpair[0] = kvpair[0].replace('\\','/')
                kvpair[1] = kvpair[1].strip('"\n')
                lines.append(kvpair)
                line = f.readline()
        if len(lines) == 0:
            print('[CONFIG] section is empty!')
            return None
        print('writing configuration to camera...')
        self.command('POLLOFF')
        self.command('CLEARCONFIG')
        line_number = 0
        for k, v in lines:
            self.set_configuration_line(line_number, k, v)
            line_number += 1
        self.command('POLLON')
        print('\rapplying configuration....')
        self.command('LOADTIMING')
        self.command('APPLYALL')
        print('camera configuration complete.')
        self.cfg = lines
        return lines
 
    def get_configuration(self):
        lines = []
        while True:
            line = self.read_configuration_line(len(lines))
            if line == '':
                break
            else:
                lines.append(line.split('=', 1))
        self.cfg = lines
        return lines
    
    def load_parameters(self):
        return self.command('LOADPARAMS').decode()
        
    def load_parameter(self, p):
        """Load a single parameter.
        
        e.g.:
          ccd = Archon()
          ccd.open()
          ccd.set_configuration_line(601, 'PARAMETER1', 'Exposures=1')
          ccd.load_parameter('Exposures')
        """
        return self.command('LOADPARAM {:s}'.format(p))
    
    def expose_frame(self, int_s, no_int_s, fn_prefix):
        """Expose a frame, all times in floating point seconds
        
        Light time = int_s s
        Dark time = int_s + no_int_s s
        
        Assumes ContinuousExposures=0 and no exposure currently running.
        Use print_relevant_parameters and newest to check
        """
        int_ms = int_s * 1000
        int_us = (int_ms - int(int_ms)) * 1000
        int_ms = int(int_ms)
        int_us = int(int_us)
        no_int_ms = no_int_s * 1000
        no_int_us = (no_int_ms - int(no_int_ms)) * 1000
        no_int_ms = int(no_int_ms)
        no_int_us = int(no_int_us)
        lastframe, lastbuf, _, _, _ = self.newest()
        ims_p = list(self.find_parameter_in_configuration(self.cfg, 'int_ms'))
        nims_p = list(self.find_parameter_in_configuration(self.cfg, 'noint_ms'))
        ius_p = list(self.find_parameter_in_configuration(self.cfg, 'int_us'))
        nius_p = list(self.find_parameter_in_configuration(self.cfg, 'noint_us'))
        ims_p[2] = 'int_ms={:d}'.format(int(int_ms))
        nims_p[2] = 'noint_ms={:d}'.format(int(no_int_ms))
        ius_p[2] = 'int_us={:d}'.format(int(int_us))
        nius_p[2] = 'noint_us={:d}'.format(int(no_int_us))
        nexp_p = list(self.find_parameter_in_configuration(self.cfg, 'Exposures'))
        nexp_p[2] = 'Exposures=0'
        self.set_configuration_line(*nexp_p)
        self.load_parameter('Exposures')
        self.set_configuration_line(*ims_p)
        if nims_p[1] is not None:
            self.set_configuration_line(*nims_p)
        if ius_p[1] is not None:
            self.set_configuration_line(*ius_p)
        if nius_p[1] is not None:
            self.set_configuration_line(*nius_p)
        # load_parameters triggers an incorrect exposure, use load_parameter
        self.load_parameter('int_ms')
        if nims_p[1] is not None:
            self.load_parameter('noint_ms')
        if ius_p[1] is not None:
            self.load_parameter('int_us')
        if nius_p[1] is not None:
            self.load_parameter('noint_us')
        nexp_p[2] = 'Exposures=1'
        self.set_configuration_line(*nexp_p)
        # Trigger an actual exposure
        self.load_parameter('Exposures')
        print('exposing/reading', end='')
        while True:
            print('.', end='', flush=True)
            if self.status_filename is not None:
                self.append_status(self.status_filename)
            frame, buf, framew, frameh, samplemode = self.newest()
            if frame != lastframe:
                break
            sleep(self.update_rate)
        self.command('LOCK{:d}'.format(buf+1))
        framesize = framew * frameh
        if samplemode:  # if samplemode 32-bit data, else 16-bit data
            framesize *= 4
        else:
            framesize *= 2
        linesize = self.BURST_LEN
        lines = (framesize + linesize - 1) // linesize
        ref = self.msgref
        self.send('FETCH{:08X}{:08X}'.format(((buf+1) | 4) << 29, lines))
        print(' fetching', end='')
        with open('{:s}_{:0.6f}s_{:0.6f}s_{:05d}x{:05d}_{:d}.raw'.format(fn_prefix, int_s, no_int_s, framew, frameh,
                                                                         frame), 'wb') as f:
            bytesremaining = framesize
            for i in range(lines):
                self.msgref = ref
                f.write(self.binrecv()[0:min(linesize, bytesremaining)])
                bytesremaining -= linesize
                if i % 10000 == 0:
                    print('.', end='')
        self.msgref = (self.msgref+1) % 256


    def find_parameter_in_configuration(self, cfg, param):
        i = 0
        for k, v in cfg:
            if  k.startswith('PARAMETER') and v.startswith(param+'='):
                return i, k, v
            else:
                i += 1
        return i, None, None
        
    def print_exposure_parameters(self, cfg):
        for p in ['IntMS', 'NoIntMS', 'ContinuousExposures', 'Exposures']:
            print(self.find_parameter_in_configuration(cfg, p))

                  
if __name__ == '__main__':
    
    # Connect to Archon
    ccd = Archon(2.0, 'statusfile.txt')
    ccd.open()
    cfg = ccd.get_configuration()
    ccd.print_exposure_parameters(cfg)
    
    # Test script here
    # but usually a separate script is run importing this module
    
    # close Archon
    ccd.close()
