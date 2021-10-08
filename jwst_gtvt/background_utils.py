"""Utils script to facilitate use of jwst-backgrounds tool in jwst-gtvt/mtvt
"""

import math

from astropy.time import Time
from astropy.table import Table
from datetime import date, datetime, timedelta
from jwst_backgrounds import jbt
import numpy as np


def build_background_table(bkg_data):
    """Build table compatible with visibility tool.

    Parameters
    ----------
    bkg_data : obj
        jwst-backgrounds object

    Returns
    -------
    bkg_table : astropy.table.Table
        Table with date(MM-DD-YYYY) and background value
    """

    bkg_table_data = {'Date':[], 'bkg':[]}

    for i, calendar_day in enumerate(bkg_data.bkg_data['calendar']):
        for year in ["2021", "2022", "2023", "2024"]:
            date = convert_day_to_date(str(calendar_day), year) # Convert date to MM-DD-YYYY
            bkg_table_data['Date'].append(date)
            bkg_table_data['bkg'].append(bkg_data.bathtub['total_thiswave'][i])

    bkg_table = Table(data=bkg_table_data)
    return bkg_table


def get_moving_target_background(ra, dec, date):
    """Calculate and get background for each RA and DEC of moving target.
    """
    bkg_table = construct_target_background_for_gtvt(ra, dec)

    year, month, day = jd_to_date(date)
    date_to_match = datetime(year, month, int(day), 0, 0)

    if date_to_match in bkg_table['Date']:
        index = np.where(bkg_table['Date'] == date_to_match)
        bkg_value = bkg_table['bkg'][index].data[0]
    else:
        bkg_value = 0.0

    return date_to_match, bkg_value


def compare_and_match_visibility_and_background_data(bkg_table, visibility_table):
    """Iterate dates in visibility table and match if possible to dates in
    backgroud table. 
    """

    bkg_column = []

    # If target is not visibile, all columns will contain nan for that date
    # so to obtain the index, I am just selecting one column to find the indices
    # of the non nan values.
    
    # no_nan_index = ~np.isnan(visibility_table['NIRCam min'])
    # mytable_no_nan = visibility_table[no_nan_index]

    # mytable_no_nan = visibility_table

    for date in visibility_table['Date']:
        if date in bkg_table['Date']:
            index = np.where(bkg_table['Date'] == date)[0]
            bkg_value_at_date = bkg_table['bkg'][index][0]
            bkg_column.append(bkg_value_at_date)
        else:
            bkg_column.append(0.0)

    visibility_table['bkg'] = bkg_column

    return visibility_table


def construct_target_background_for_gtvt(ra, dec, wavelength=2.15, thresh=1.1):
    """Calculate and construct table with background data for target.

    Parameters
    ----------
    ra : float
        Right ascention of target
    dec : float
        Declination of target
    wavelength : float
        Wavelength at which the bathtub curve is calculated, in micron
    thresh : float
        the background threshold, relative to the minimum.  Default=1.1, which corresponds to <5% above the minimum background noise.
        Note that the actual noise difference will be even smaller, as there are often other significant sources of noise than just the
        background (source shot noise, detector noise, etc.).

    Returns
    -------
    bkg_table : astropy.table.Table
        Astropy table object with columns date and bkg
    """

    bkg_data = jbt.background(ra, dec, wavelength, thresh)
    bkg_table = build_background_table(bkg_data)

    return bkg_table


def convert_day_to_date(day_number, year):
    """Convert day number of the year to MM-DD-YYYY.

    day_number : str
        Day number of the year

    year : str
        Year of date conversion
    """
    # Initializing start date
    strt_date = date(int(year), 1, 1)
    
    # converting to date
    res_date = strt_date + timedelta(days=int(day_number) - 1)
    res = res_date.strftime(("%Y-%m-%d"))

    formatted_date = datetime.strptime(res, '%Y-%m-%d')

    return formatted_date


# Taken from jdutil.py
def jd_to_date(jd):
    """
    Convert Julian Day to date.
    
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
        4th ed., Duffet-Smith and Zwart, 2011.
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
        
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    
    """
    jd = jd + 0.5
    
    F, I = math.modf(jd)
    I = int(I)
    
    A = math.trunc((I - 1867216.25)/36524.25)
    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    
    D = math.trunc((C - 122.1) / 365.25)
    
    E = math.trunc(365.25 * D)
    
    G = math.trunc((C - E) / 30.6001)
    
    day = C - E + F - math.trunc(30.6001 * G)
    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715

    return year, month, day

