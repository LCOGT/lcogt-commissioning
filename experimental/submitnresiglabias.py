import lcocommissioning.common.common as common
import datetime as dt
import logging
_logger = logging.getLogger(__name__)
logging.basicConfig(level=getattr(logging, 'DEBUG'),
                    format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
start = dt.datetime.utcnow() + dt.timedelta(minutes=1)
end = start + dt.timedelta(minutes=2)
start = str(start).replace(' ', 'T')
end = str(end).replace(' ', 'T')


request = {

    "molecules": [
        {

            "tag_id": "LCOGT",
            "user_id": "daniel_harbeck",
            "prop_id": "engineering test",
            "group": "testflat",
            "exposure_count": 1,
            "exposure_time": 30,
            "bin_x": 1,
            "bin_y": 1,
            "inst_name": "nres02",
            "priority": 1,
            "type": "LAMP_FLAT",
            "ag_filter": "",
            "exposure_time": "60",
            "spectra_slit": "slit_1.6as",

        }
    ],
    "start": start,
    "end": end,
    "site": "elp",
    "observatory": "igla",
    "telescope": "1m0a",
    "instrument_class": "1M0-NRES-SCICAM",
    "length": 0,
    "priority": 10,


}

print ("Need to update this program for direct submission. Quitting.")
exit(1)

common.submit_request_group (request, dosubmit=True)
