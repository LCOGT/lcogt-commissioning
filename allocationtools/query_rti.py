import datetime
import json
from  astropy.table import Table, unique
from matplotlib import pyplot as plt
import numpy as np
import requests
import matplotlib.dates as mdates



def list_of_dicts_to_dict_of_lists(data):
    if not data:
        return {}
    print (data[0])
    keys = data[0].keys()
    reformatted_data = {key: [row[key] for row in data] for key in keys}
    return reformatted_data

def query_api(endpoint_url="http://downtime.lco.gtn/api/"):
    """
    Query an API endpoint with a GET request.
    
    Args:
        endpoint_url (str): The API endpoint URL
        
    Returns:
        dict: JSON response from the API
    """

    params = {'reason': 'RTI',
              'site':'ogg',
              'enclosure': 'clma',
              'telescope': '2m0a',
              'starts_after' : "2023-01-01T00:00:00",}
              

    try:
        response = requests.get(endpoint_url, params=params)
        response.raise_for_status()
        data = response.json()['results'] if response is not None else None
        data = list_of_dicts_to_dict_of_lists(data)
        data['start'] = [datetime.datetime.fromisoformat(dt[:-1]) for dt in data['start']]
        myTable = Table(data)
        myTable['startdate'] = np.asarray([dt.date() for dt in myTable['start']])
        myTable['starttime'] = np.asarray([datetime.datetime.combine(datetime.date(1900, 1, 1), dt.time()) for dt in myTable['start']])
        myTable['commonstartdate'] = np.asarray([datetime.date(2024, dt.month, dt.day) for dt in myTable['start']])
        return myTable
    except requests.exceptions.RequestException as e:
        print(f"Error querying API: {e}")
        return None



def print_reduced_datatable(datatable):

    reduced_table = unique(datatable, keys='commonstartdate')
  
    # interpolate
    check_date = datetime.date(2024, 1, 1)
    while check_date <= datetime.date(2024,12,31):
        if check_date not in reduced_table['commonstartdate']:

            differences = (reduced_table['commonstartdate'] - check_date)
            
            closest_index = np.argmin(np.abs(differences))
            
            new_row = {'commonstartdate': check_date,
                       'starttime': reduced_table['starttime'][closest_index],  # placeholder, will adjust below
                       'site': reduced_table['site'][0].upper()}  # assuming site is the same for all rows
            print ("interpolating:", new_row)
            reduced_table.add_row(new_row)
        check_date += datetime.timedelta(days=1)

    reduced_table.sort('commonstartdate') 


    site = reduced_table['site'][0]
    reduced_table.write(f'reduced_startdates_{site}.csv', format='csv', overwrite=True, include_names=['commonstartdate', 'starttime', 'site'])



def plot_allstarts(datatable):
    # x: dates (date objects), y: times represented on a fixed reference date so we can format them as HH:MM
    fig, ax = plt.subplots(figsize=(10, 6))
    for obs in ['coj', 'ogg']:
        mask = datatable['site'] == obs          
        ax.scatter(datatable['commonstartdate'][mask], datatable['starttime'][mask], marker='o', label = obs)

    # format x axis as dates
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d'))
    fig.autofmt_xdate(rotation=45)

    # format y axis to show time of day
    ax.yaxis.set_major_locator(mdates.HourLocator(interval=1))
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    ax.set_xlabel('Date')
    ax.set_ylabel('Time (UTC)')
    ax.set_title('Start times by date')
    ax.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('startimes.png')



if __name__ == "__main__":
    # Example usage
    data = query_api()
    plot_allstarts(data)
    print_reduced_datatable(data)
  