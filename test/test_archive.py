from lcocommissioning.common.lco_archive_utilities import get_auto_focus_frames
import logging
logging.basicConfig(level=getattr(logging, 'DEBUG'),
                    format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')


def test_image():
    result = get_auto_focus_frames(2896481)

    efs = [x for x in result if 'ef' in x['basename']]
    fas = [x for x in result if 'ef' in x['basename']]

    assert len (result) == 18

    assert len (efs) == len(fas)

    for i in efs:
        print (i['basename'])
    print ('---------------------')
    for i in fas:
        print (i['basename'])
