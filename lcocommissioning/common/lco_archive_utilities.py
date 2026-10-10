import datetime
import glob
import tempfile
import logging
import os
import time

import numpy as np
from numpy import rec
import requests
from astropy.io import fits
from astropy.table import Table
from opensearch_dsl import Search
from opensearchpy import OpenSearch
from opensearchpy import exceptions as opensearch_exceptions
from opensearchpy.helpers import ScanError

from lcocommissioning.common.common import lco_site_lonlat

log = logging.getLogger(__name__)
# The client logs every failed request with a traceback at WARNING level; failures are reported by scan_opensearch.
logging.getLogger('opensearch').setLevel(logging.ERROR)
logging.getLogger('connectionpool').setLevel(logging.WARNING)

ARCHIVE_ROOT = "/archive"
ARCHIVE_API_TOKEN = os.getenv('ARCHIVE_API_TOKEN', '')


class ArchiveDiskCrawler:
    ''' Legacy code from the good old times (2019) when the /archive mount was accessible, and everybody was happy in
    the file system land.

    Now everybody moved on to databases, object stores, and other witchcraft.

    '''

    archive_root = None

    def __init__(self, rootdirectory=ARCHIVE_ROOT):
        self.archive_root = rootdirectory

    def find_cameras(self, sites=lco_site_lonlat, cameras=["fa??", "fs??", "kb??"]):
        sitecameras = []
        for site in sites:
            for camera in cameras:
                dir = "{}/{}/{}".format(self.archive_root, site, camera)
                candiates = glob.glob(dir)
                sitecameras.extend(candiates)
        return sitecameras

    @staticmethod
    def get_last_n_days(lastNdays):
        ''' Utility method to return the last N days as string YYYYMMDD, nicely arranged in an array.'''

        date = []
        today = datetime.datetime.utcnow()
        for ii in range(lastNdays):
            day = today - datetime.timedelta(days=ii)
            date.append(day.strftime("%Y%m%d"))
        return date[::-1]

    @staticmethod
    def findfiles_for_camera_dates(sitecamera, date, raworprocessed, filetempalte, prefix=""):
        dir = "{}{}/{}/{}/{}".format(prefix, sitecamera, date, raworprocessed, filetempalte)
        files = glob.glob(dir)
        if (files is not None) and (len(files) > 0):
            myfiles = np.asarray([[f, "-1"] for f in files])
            return Table(myfiles, names=['FILENAME', 'FRAMEID'])
        return None


def make_opensearch(index, filters, queries=None, exclusion_filters=None, range_filters=None, prefix_filters=None,
                    terms_filters=None,
                    opensearch_url='https://opensearch.lco.global'):
    """
    Make an ElasticSearch query

    Parameters
    ----------
    index : str
            Name of index to search
    filters : list of dicts
              Each dict has a criterion for an OpenSearch "filter"
    queries : list of dicts
              Each dict has a "type" and "query" entry. The 'query' entry is a dict that has a criterion for an
              ElasticSearch "query"
    exclusion_filters : list of dicts
                        Each dict has a criterion for an OpenSearch "exclude"
    range_filters: list of dicts
                   Each dict has a criterion an openSearch "range filter"
    opensearch_url : str
             URL of the openSearch host

    Returns
    -------
    search : opensearch_dsl.Search
             The OpenSearch object. Note the query is only sent to the server when the search is executed, e.g.,
             via .scan() or .execute(); connection and server side errors are raised there.

    Raises
    ------
    ValueError
             If the index or any of the filter / query definitions is malformed.
    """
    if not isinstance(index, str) or not index:
        raise ValueError(f"OpenSearch index must be a non-empty string, got {index!r}")
    if not isinstance(opensearch_url, str) or not opensearch_url:
        raise ValueError(f"OpenSearch URL must be a non-empty string, got {opensearch_url!r}")

    filters = _validate_opensearch_filters('filters', filters)
    terms_filters = _validate_opensearch_filters('terms_filters', terms_filters)
    range_filters = _validate_opensearch_filters('range_filters', range_filters)
    prefix_filters = _validate_opensearch_filters('prefix_filters', prefix_filters)
    exclusion_filters = _validate_opensearch_filters('exclusion_filters', exclusion_filters)
    queries = [] if queries is None else list(queries)
    for q in queries:
        if not isinstance(q, dict) or 'type' not in q or not isinstance(q.get('query'), dict):
            raise ValueError(f"OpenSearch query must be a dict with 'type' and a 'query' dict, got {q!r}")

    try:
        opensearch = OpenSearch(opensearch_url)
        s = Search(using=opensearch, index=index)
        for f in filters:
            s = s.filter('term', **f)
        for f in terms_filters:
            s = s.filter('terms', **f)
        for f in range_filters:
            s = s.filter('range', **f)
        for f in prefix_filters:
            s = s.filter('prefix', **f)
        for f in exclusion_filters:
            s = s.exclude('term', **f)
        for q in queries:
            s = s.query(q['type'], **q['query'])
    except Exception as e:
        log.error(f"Failed to build OpenSearch query on index {index} at {opensearch_url}: {e}")
        raise

    log.debug(f"OpenSearch query on index {index} at {opensearch_url}: {s.to_dict()}")
    return s


class OpenSearchQueryError(Exception):
    """ An OpenSearch query could not be executed; the message describes the query and the cause."""


# HTTP status codes of transient server side conditions worth retrying.
_OPENSEARCH_RETRY_STATUS = (429, 502, 503, 504)


def scan_opensearch(search, description, opensearch_url='https://opensearch.lco.global', retries=2, retry_delay=10):
    """
    Execute an OpenSearch search (e.g., from make_opensearch) and return all matching documents.

    Transient errors (connection problems, timeouts, server overloaded) are retried up to `retries` times, waiting
    `retry_delay` seconds between attempts. All other errors, and transient errors that persist, are logged and
    raised as OpenSearchQueryError with a human readable explanation.

    Parameters
    ----------
    search : opensearch_dsl.Search
             The search to execute
    description : str
             What the query is for, used in log and error messages, e.g. "noise/gain frames for 20260817"
    opensearch_url : str
             URL of the OpenSearch host, used in log and error messages

    Returns
    -------
    hits : list
           All documents matching the query
    """
    for attempt in range(retries + 1):
        try:
            # scan() is a lazy generator; errors can occur at any page, so read all results here.
            hits = list(search.scan())
            log.debug(f"OpenSearch query for {description} returned {len(hits)} documents")
            return hits
        except Exception as e:
            reason, transient = _explain_opensearch_error(e)
            message = f"OpenSearch query for {description} at {opensearch_url} failed: {reason}"
            if transient and attempt < retries:
                log.warning(f"{message}; retrying in {retry_delay} s (attempt {attempt + 1} of {retries})")
                time.sleep(retry_delay)
                continue
            log.error(message)
            raise OpenSearchQueryError(message) from e


def _explain_opensearch_error(e):
    """ Return (human readable reason, is transient) for an exception raised while executing an OpenSearch query."""
    if isinstance(e, opensearch_exceptions.ConnectionTimeout):
        return f"connection timed out ({e.error})", True
    if isinstance(e, opensearch_exceptions.SSLError):
        return f"SSL error ({e.error})", False
    if isinstance(e, opensearch_exceptions.ConnectionError):
        return f"cannot connect to server ({e.error})", True
    if isinstance(e, opensearch_exceptions.AuthenticationException):
        return "authentication failed (HTTP 401)", False
    if isinstance(e, opensearch_exceptions.AuthorizationException):
        return "access denied (HTTP 403)", False
    if isinstance(e, opensearch_exceptions.NotFoundError):
        return f"index not found (HTTP 404: {e.error})", False
    if isinstance(e, opensearch_exceptions.RequestError):
        return f"query rejected by server (HTTP 400: {e.error} {e.info})", False
    if isinstance(e, opensearch_exceptions.TransportError):
        return f"server error (HTTP {e.status_code}: {e.error})", e.status_code in _OPENSEARCH_RETRY_STATUS
    if isinstance(e, ScanError):
        return f"scan incomplete, some shards failed ({e})", True
    return f"{type(e).__name__}: {e}", False


def _validate_opensearch_filters(name, filters):
    """ Check that filters is a list of non-empty dicts without None values; return it as a list."""
    if filters is None:
        return []
    if isinstance(filters, dict):
        raise ValueError(f"OpenSearch {name} must be a list of dicts, not a single dict: {filters!r}")
    filters = list(filters)
    for f in filters:
        if not isinstance(f, dict) or len(f) == 0:
            raise ValueError(f"OpenSearch {name} entries must be non-empty dicts, got {f!r}")
        nonevalues = [key for key, value in f.items() if value is None]
        if nonevalues:
            raise ValueError(f"OpenSearch {name} entry {f!r} has no value for {nonevalues}")
    return filters


def get_frames_for_noisegainanalysis(dayobs, site=None, cameratype=None, camera=None, readmode=['full_frame'],
                                     obstype=['BIAS', 'SKYFLAT'], rlevel = None, filter=None, isMaster=False, opensearch_url='https://opensearch.lco.global'):
    """ Queries for a list of processed LCO images that are viable to get a photometric zeropoint in the griz bands measured.

        Selection criteria are by DAY-OBS, site, by camera type (fs,fa,kb), what filters to use, and minimum exposure time.
        Only day-obs is a mandatory fields, we do not want to query the entire archive at once.
     """
    
    # terms filters need lists; accept single values as well.
    if isinstance(readmode, str):
        readmode = [readmode]
    if isinstance(obstype, str):
        obstype = [obstype]
    if isinstance(filter, str):
        filter = [filter]

    query_filters = [{'DAY-OBS': dayobs}, {'RLEVEL': 0 if not rlevel else rlevel}, ]
    range_filters = []
    terms_filters = [{'OBSTYPE': obstype}, {'CONFMODE': readmode}]
    prefix_filters = []

    if site is not None:
        query_filters.append({'SITEID': site})
    if camera is not None:
        query_filters.append({'INSTRUME': camera})
    if cameratype is not None:
        prefix_filters.append({'INSTRUME': cameratype})
    if isMaster:
        query_filters.append({'ISMASTER': True})
    if filter is not None:
        terms_filters.append({'FILTER': filter})
    
    
    queries = []
    search = make_opensearch('fitsheaders', query_filters, queries, exclusion_filters=None,
                             opensearch_url=opensearch_url,
                             range_filters=range_filters, prefix_filters=prefix_filters,
                             terms_filters=terms_filters)
    records = scan_opensearch(search, f"frames of DAY-OBS {dayobs} (camera {camera}, camera type {cameratype}, "
                                      f"site {site}, readmode {readmode}, obstype {obstype})", opensearch_url)

    records_sanitized = [[record['filename'], record['SITEID'], record['INSTRUME'], record['RLEVEL'], record['DAY-OBS'],
                          record['frameid'],  record['FILTER'], record['CONFMODE']] for record in records]

    t = Table(np.asarray(records_sanitized), names=('FILENAME', 'SITEID', 'INSTRUME', 'RLEVEL', 'DAY-OBS', 'frameid', 'FILTER', 'CONFMODE'))
    return t


def get_muscat_focus_request_ids(muscat, before=None, after=None, opensearch_url='https://opensearch.lco.global'):
    log.debug("Starting opensearch query for Muscat focus requests")

    query_filters = [{'INSTRUME': 'ep07' if muscat == 'mc04' else 'ep02'},
                     {'OBJECT':'auto_focus'},
                     {'OBSTYPE': 'EXPOSE'},
                     ]
    range_filters = [{'DATE-OBS' :  {
        "gte": after,
        "lte": before
        } } ]
    terms_filters = []
    prefix_filters = []

    queries = []
    search = make_opensearch('fitsheaders', query_filters, queries, exclusion_filters=None,
                             opensearch_url=opensearch_url,
                             range_filters=range_filters, prefix_filters=prefix_filters,
                             terms_filters=terms_filters)
    records = scan_opensearch(search, f"auto focus requests of {muscat} between {after} and {before}", opensearch_url)
    records_sanitized = [int(record['REQNUM']) if record['REQNUM'] is not None else 0  for record in records]
    records_sanitized = np.unique(records_sanitized)
    records_sanitized = records_sanitized[ records_sanitized > 0]
    return records_sanitized


def filename_to_archivepath_dict(filenametable, rootpath=ARCHIVE_ROOT):
    ''' Return a dictionary with camera -> list of FileIO-able path of imagers from an elastic search result.
        We are still married to /archive file names here - because reasons. Long term we should go away from that.
    '''

    cameras = set(filenametable['INSTRUME'])
    returndict = {}
    for camera in cameras:
        returndict[camera] = [['{}/{}/{}/{}/{}/{}'.format(rootpath, record['SITEID'], record['INSTRUME'],
                                                          record['DAY-OBS'], 'raw', record['FILENAME']),
                               record['frameid']] for record in filenametable[filenametable['INSTRUME'] == camera]]

        returndict[camera] = Table(np.asarray(returndict[camera]), names=['FILENAME', 'frameid'])
    return returndict


def get_frames_from_request(requestid):
    url = f'https://archive-api.lco.global/frames/'
    params = {'request_id': requestid,
              'reduction_level': 0,
              }
    headers = {'Authorization': 'Token {}'.format(ARCHIVE_API_TOKEN)}
    response = requests.get(url, headers=headers, params=params)
    response.raise_for_status()
    return response.json()


def get_auto_focus_frames(requestid):
    candidates = get_frames_from_request(requestid)
    focusimagelist = []
    for imageinfo in candidates['results']:
        if 'x00' in imageinfo['basename']:
            focusimagelist.append(
                {'basename': imageinfo['basename'], 'id': imageinfo['id'], 'INSTRUME': imageinfo['INSTRUME']})

    return focusimagelist


def download_from_archive(frameid):
    """
    Download a file from the LCO archive by frame id.
    param frameid: Archive API frame ID
    return: Astropy HDUList
    """
    url = f'https://archive-api.lco.global/frames/{frameid}'
    log.debug("Downloading image frameid {} from URL: {}".format(frameid, url))
    headers = {'Authorization': 'Token {}'.format(ARCHIVE_API_TOKEN)}
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    response_dict = response.json()
    if response_dict == {}:
        log.warning("No file url was returned from id query")
        raise Exception('Could not find file remotely.')
    frame_url = response_dict['url']
    log.debug(frame_url)
    # Spool the file to disk instead of holding it in memory: fits files can be large, and astropy would keep
    # another copy of the (compressed) data in memory when reading from an in-memory buffer.
    with requests.get(frame_url, stream=True) as file_response:
        file_response.raise_for_status()
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as tmp:
            for chunk in file_response.iter_content(chunk_size=1 << 20):
                tmp.write(chunk)
    try:
        f = fits.open(tmp.name)
    finally:
        # The open file handle keeps the data accessible until the HDUList is closed; the file is removed then.
        os.unlink(tmp.name)
    return f


if __name__ == '__main__':

    camera = 'fa15'
    dates = ArchiveDiskCrawler.get_last_n_days(3)

    for dayobs in dates:
        listofframes = get_frames_for_noisegainanalysis(dayobs, camera='ep60', readmode=['full_frame'])
        filelist = filename_to_archivepath_dict(listofframes)
        print("{} {} ".format(dayobs, filelist.keys()))

    reqids = get_muscat_focus_request_ids('mc04')
    print (reqids)
    # c = ArchiveCrawler()
    # for dayobs in dates:
    #     listofframes = c.findfiles_for_camera_dates("/archive/engineering/lsc/{}".format(camera), dayobs, 'raw', "*[xbf]00.fits*")
    #     #print ("{} {} ".format (dayobs, listofframes))
