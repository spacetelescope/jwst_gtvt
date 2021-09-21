"""Utils script to facilitate use of jwst-backgrounds tool in jwst-gtvt/mtvt
"""

from astropy.time import Time
from astropy.table import Table
from datetime import date, timedelta
from jwst_backgrounds import jbt
import numpy as np
import pandas as pd


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
            bkg_table_data['Date'].append(Time(date))
            bkg_table_data['bkg'].append(bkg_data.bathtub['total_thiswave'][i])

    bkg_table = Table(data=bkg_table_data)
    return bkg_table


def compare_and_match_visibility_and_background_data(bkg_table, visibility_table):
    """Iterate dates in visibility table and match if possible to dates in
    backgroud table. 
    """

    bkg_column = []

    has_nan = np.zeros(len(visibility_table), dtype=bool)
    for col in visibility_table.itercols():
        if col.info.dtype.kind == 'f':
            has_nan |= np.isnan(col)
    mytable_no_nan = visibility_table[~has_nan]

    for date in mytable_no_nan['Date']:
        if date in bkg_table['Date']:
            index = np.where(bkg_table['Date'] == date)[0]
            bkg_value_at_date = bkg_table['bkg'][index][0]
            bkg_column.append(bkg_value_at_date)
        else:
            bkg_column.append(0.0)

    mytable_no_nan['bkg'] = bkg_column
    
    return mytable_no_nan


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
    
    return res
