import requests
from lcocommissioning.common.common import VALHALLA_API_TOKEN

data = {
    'name': 'NRES ThAr  test',
    'proposal': 'LCOEngineering',
    'site': 'tlv',
    'enclosure': 'doma',
    'telescope': '1m0a',
    'start': '2023-08-01T18:45:00',
    'end': '2023-08-01T19:32:00',
    'request': {
        'acceptability_threshold': 100,
        'configurations': [
            {
                'type': 'ARC',
                'instrument_type': '1M0-NRES-SCICAM',
                'target': {
                    'type': 'HOUR_ANGLE',
                    'name': 'none',
                    'hour_angle': 1.0,
                    'dec': 0.0
                },
                'acquisition_config': {
                    'mode': 'OFF'
                },
                'guiding_config': {
                    'mode': 'ON',
                    'optional': True
                },
                'constraints': {},
                'extra_params': {
                #    'lamp': 'ThAr Master'
                },
                'instrument_configs': [{
                    'exposure_time': 300,
                    'exposure_count': 3,
                    'mode': '',
                    'bin_x': 1,
                    'bin_y': 1,
                    'rotator_mode': '',
                    'optical_elements': {}
                }]
            }
        ]
    }
}



url ='http://internal-observation-portal.lco.gtn/api/schedule/'
headers =  {'Authorization': 'Token {}'.format(VALHALLA_API_TOKEN)}
print (headers)
response = requests.post(url, json=data, headers=headers)
print (response.content)

