import numpy as np
import requests
from astropy.io import ascii
from astropy.stats import sigma_clip
from glob import glob
from astropy.table import Table, unique
from astropy.timeseries import TimeSeries, aggregate_downsample
from astropy.time import Time
import astropy.units as u

data_dir = '/home/soichiro/'

# Functions for ZTF
def to_ztf_lc(time, mag, magerr=None, time_col_name='mjd', time_col_format='mjd'):
    time, mag = np.array(time), np.array(mag)
    if magerr is None:
        magerr = np.full(time.size, np.nan)
    lc = Table([time*u.day, mag*u.mag, magerr*u.mag], names=(f'{time_col_name}', 'mag', 'magerr'))  # pylint: disable=no-member
    lc.add_column(Time(lc[f'{time_col_name}'], format=f'{time_col_format}'), name='time', index=0)
    lc = TimeSeries(data=lc)
    return lc


def load_ztf_lc(source=None, ra=None, dec=None, data_dir=data_dir, 
               sort=True, verbose=False, add_tess_time_column=True, use_tess_time_as_index=False,
               time_format='mjd', time_scale='tdb'):
    """
    """
    if source is not None:
        if (ra is not None) or (dec is not None):
            print(f'Please only supply either the ra, dec arguments OR the source argument, but not both.')
            return
        ra, dec = source['ra'][0], source['dec'][0]
    search_path = f'{data_dir}/light_curves/ra_{ra}_dec_{dec}'
    if verbose:
        print(f'Checking for file \'{search_path}\'')
    path = glob(search_path)
    assert len(path) == 1
    
    lc = ascii.read(path[0])
    if sort == True:
        sorted_idx = lc['mjd'].argsort()
        lc = lc[sorted_idx]
    if add_tess_time_column == True:
        lc = insert_tess_time_column(lc)
    
    time_column = 'mjd'
    # scale = 'tdb' solves the warning issue. I do not know why but whatever.
    if use_tess_time_as_index == True:
        time_column = 'tess_time'
    lc.add_column(Time(lc[f'{time_column}'], format=f'{time_format}', scale=f'{time_scale}'), name='time', index=0)
    lc = TimeSeries(data=lc)
#     lc = TimeSeries(data=lc, time=Time(lc['mjd'], format='mjd'))
    return lc


def insert_tess_time_column(lc, index=None, name='tess_time'):
    """FIXME: I need to understand the nuances of the date conversion
              since ZTF is ground-based while TESS is space-based."""
    """
    Add a column that contains the date in TESS convention (BJD - 2457000)
    """
    
    offset = 2457000 * u.day  # pylint: disable=no-member
    tess_time = lc['hjd'] - offset
    if index is None:
        index = lc.colnames.index('mjd')+1
    lc.add_column(tess_time, index=index, name=name)
    return lc


def add_flux_column(lc, name='flux', mag_zp=0):
    """
    Add a column that converts the recorded ZTF magnitude to (relative) flux.
    """
    
    mag = lc['mag'].value
    exponent = (mag_zp-mag) / 2.0
    flux = 10**exponent
    index = lc.colnames.index('magerr')+1
    lc.add_column(flux, index=index, name=name)
    return lc

def objectID_filter(lc):
    """
    Filters the ZTF light curve based on Object ID.
    
    As ZTF light curve files can contain data points from multiple sources this function filters
    the light curve file so that only data points associated with the ObjectID with the most number 
    of measurements are returned. 
    """
    
    oid, n_obs = np.unique(lc['oid'], return_counts=True)
    oid_with_max_obs = oid[n_obs.argmax()]
    return lc[lc['oid'] == oid_with_max_obs]


def split_into_seasons(lc):
    """
    Split ZTF light curve into it's three observational 'seasons'.
    """
    
    mjd = lc['mjd']
    s1 = mjd < 58500*u.day  # pylint: disable=no-member
    s2 = (mjd > 58500*u.day) * (mjd < 58850*u.day)  # pylint: disable=no-member
    s3 = mjd > 58850*u.day  # pylint: disable=no-member
    season_idx = [s1, s2, s3]
    season_lcs = [lc[s] for s in season_idx]
    return season_lcs


def preprocess_ztf_lc(lc, remove_outliers=True, sigma=3, remove_trend=True, degree=1, time_column='mjd'):
    """
    Sigma clip outliers and fit & subtract a polynomial to remove trends in the mag.
    """
    
    outliers = None
    med_mag = None
    p = None
    if remove_outliers == True:
        outliers = sigma_clip(lc['mag'], sigma=sigma)
        lc = lc[~outliers.mask]
    if remove_trend == True:
        x, y = lc[f'{time_column}'], lc['mag']
        med_mag = np.median(y)
        p = np.polynomial.Polynomial
        p = p.fit(x.value, y.value, deg=degree)
        lc['mag'] = (y.value - p(x.value) + med_mag.value) * u.mag
    return (lc, {'med_mag': med_mag, 'outliers': outliers, 'p': p})


def bin_lc(lc, time_bin_size, cols=['time', 'hjd', 'mjd', 'mag', 'magerr'], aggregate_func=np.nanmean):
    """
    I should probably warn the user that only the date column used as the index will make sense.
    """
    lc_copy = lc.copy()
    colnames = lc_copy.colnames
    if 'tess_time' in colnames:
        cols.append('tess_time')
    if 'flux' in colnames:
        cols.append('flux')
    lc_copy.keep_columns(cols)
    binned_lc = aggregate_downsample(lc_copy, time_bin_size=time_bin_size, aggregate_func=aggregate_func)
    return binned_lc