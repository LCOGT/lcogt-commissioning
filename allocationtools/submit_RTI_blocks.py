import datetime
import json
from tracemalloc import start
from astropy.table import Table, unique
from matplotlib import pyplot as plt
import numpy as np
import requests
import matplotlib.dates as mdates

downtime_url = "http://downtime.lco.gtn"
thetoken = os.getenv('DOWNTIME_TOKEN')
startdate = datetime.date(2026, 2, 1)
enddate = datetime.date(2026, 7, 31)


valid_weekdays = { 'coj': [0, 1, 2, 3, 4, 5 ],  # Mon-Sa
                   'ogg': [0, 1, 2, 3, 4],      # Mon-Fri
                 }

holidays = {
    "springbreak2026": {
        'start': datetime.date(2026, 3, 30),
        'end': datetime.date(2026, 4, 11),
    }
}

def readslots(site,):
    filename = f"allocationtools/reduced_startdates_{site}.csv"
    time_templates = Table.read(filename, format='csv')
    time_templates['commonstartdate'] = [datetime.date.fromisoformat(x) for x  in time_templates['commonstartdate']]
    time_templates['starttime'] = [datetime.datetime.fromisoformat(x) for x  in time_templates['starttime']]
    return time_templates

def isHoliday(date, holidays):
    for holiday in holidays.values():
        if holiday['start'] <= date <= holiday['end']:
            return True
    return False


def submit_RTI_block(site, enclosure, telescope, starttime, endtime, reason, submit=False):
    print('Will submit:  ', starttime, '  to  ', endtime, '  at  ', site, enclosure, telescope, '  where downtime is', (endtime-starttime   ).seconds/60, 'minutes long.')
    if submit:
        data = {
                    'site': site,
                    'enclosure': enclosure,
                    'telescope': telescope,
                    'start': starttime.isoformat(),
                    'end': endtime.isoformat(),
                    'reason': reason
                }
        response = requests.post(
                    f'{downtime_url}/api/',
                    json=data,
                    headers={'Authorization': 'Token {thetoken}'}
                   )
                
        response.raise_for_status()

def define_RTI_blocks (site, templatetimes, startdate, enddate, holidays):
    totaltime = 0
    datetosubmit = startdate
    while datetosubmit <= enddate:
        if  isHoliday(datetosubmit, holidays):
            print (f"{datetosubmit} - Skipping a holiday day")
         
        elif  datetosubmit.weekday() not in valid_weekdays[site]:
            print (f"{datetosubmit} - Skipping a non-valid weekday {datetosubmit.weekday()} for site {site}")
        
        else:

       

            reduceddate = datetime.date(2024, datetosubmit.month, datetosubmit.day)
            indexintemplate = np.where(templatetimes['commonstartdate'].astype('datetime64[D]').astype(datetime.date) == reduceddate)[0][0]
       
            starttime = templatetimes['starttime'][indexintemplate]
            starttime = datetime.datetime(year=datetosubmit.year, month=datetosubmit.month, day=datetosubmit.day,
                                          hour=starttime.hour, minute=starttime.minute, second=starttime.second)
            endtime = starttime + datetime.timedelta(minutes=30)
            print (f"{datetosubmit} - Submitting RTI block for {datetosubmit.strftime('%a')} {starttime} to {endtime}")
            totaltime += 30
        datetosubmit += datetime.timedelta(days=1)

    print (f"Total allocated time for {site}: {totaltime/60} hours")
        



if __name__ == "__main__":
    for site in ['coj', 'ogg']:
        #template_times_ogg = readslots(site,)
        #define_RTI_blocks(site, template_times_ogg, startdate, enddate, holidays)
        pass
    now = datetime.datetime.utcnow()
    endtime = now + datetime.timedelta(minutes=30)
    submit_RTI_block('lsc', 'aqwb', '0m4a', now, endtime, 'Daniel Test', submit=True)