import logging
import sys
from copy import deepcopy

import astropy.time as astt
import numpy as np
from astropy.table import Table
from sqlalchemy import Column, Integer, String, Float, create_engine, or_
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

assert sys.version_info >= (3, 5)
_logger = logging.getLogger(__name__)

Base = declarative_base()

class MuscatFocusMeasurement(Base):
    __tablename__ = 'muscatnoise'

    requestid = Column(Integer, primary_key=True)
    muscat = Column(String, index=True)
    dateobs = Column(String)
    temperature = Column(Float)

    camera_g = Column(String, index=True)
    focus_g = Column(Float)
    seeing_g = Column(Float)
    error_g = Column(Float)

    camera_r = Column(String, index=True)
    focus_r = Column(Float)
    seeing_r = Column(Float)
    error_r = Column(Float)

    camera_i = Column(String, index=True)
    focus_i = Column(Float)
    seeing_i = Column(Float)
    error_i = Column(Float)

    camera_z = Column(String, index=True)
    focus_z = Column(Float)
    seeing_z = Column(Float)
    error_z = Column(Float)

    def __repr__(self):
        return f'{self.requestid} {self.muscat} {self.dateobs} {self.temperature} {self.focus_g}{self.focus_r} {self.focus_i} {self.focus_z} '


class muscatfocusdb():
    def __init__(self, fname):
        _logger.debug("Open data base file %s" % fname)
        self.engine = create_engine(fname, echo=False)
        MuscatFocusMeasurement.__table__.create(bind=self.engine, checkfirst=True)
        self.session = sessionmaker(bind=self.engine)()

    def close(self):
        """ Close the database safely"""
        _logger.debug("Closing data base session")
        self.session.close()

    def exists(self, request_id):
        entry =  self.session.query(MuscatFocusMeasurement).filter_by(requestid=int(request_id)).first()
        _logger.info (f"Entry query {request_id} returns {entry}")
        return entry

    def addMeasurement(self, m, commit=True):
        existingEntry = self.exists(m.requestid)
        if (existingEntry):
            pass
            # TODO: etc...
        else:
            self.session.add(m)
        if commit:
            self.session.commit()

    def getMeasurementsForMuscat(self, muscat=None):

        q = self.session.query(MuscatFocusMeasurement).filter(MuscatFocusMeasurement.muscat == muscat)

        allrows = [[e.requestid, e.muscat, e.dateobs, e.temperature, e.camera_g, e.focus_g, e.seeing_g, e.error_g,
                    e.camera_r, e.focus_r, e.seeing_r, e.error_r,
                    e.camera_i, e.focus_i, e.seeing_i, e.error_i,
                    e.camera_z, e.focus_z, e.seeing_z, e.error_z]
                   for e in q.all()]

        t = Table(np.asarray(allrows),
                  names=['requestid', 'muscat', 'dateobs', 'temperature', 'camera_g', 'focus_g', 'seeing_g', 'error_g',
                         'camera_r', 'focus_r', 'seeing_r', 'error_r',
                         'camera_i', 'focus_i', 'seeing_i', 'error_i',
                         'camera_z', 'focus_z', 'seeing_z', 'error_z', ])


        t['dateobs'][t['dateobs'] == None] = '1900-01-01'
        _logger.info (t['dateobs'])
        t['dateobs'] = t['dateobs'].astype(str)
        t['dateobs'] = astt.Time(t['dateobs'], scale='utc', format=None).to_datetime()
        t['temperature'] = t['temperature'].astype (float)
        t['focus_g'] = t['focus_g'].astype (float)
        t['focus_r'] = t['focus_r'].astype (float)
        t['focus_i'] = t['focus_i'].astype (float)
        t['focus_z'] = t['focus_z'].astype (float)
        t['seeing_r'] = t['seeing_r'].astype (float)
        t['requestid'] = t['requestid'].astype (int)
        return t


if __name__ == '__main__':
    # TODO: Move this stuff into a test routine
    logging.basicConfig(level=getattr(logging, 'DEBUG'),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    c = muscatfocusdb('sqlite:///muscatfocus.sqlite')

    print(c.getMeasurementsForMuscat('mc04'))
    c.close()
