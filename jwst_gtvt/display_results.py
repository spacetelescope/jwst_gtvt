import datetime
from collections import defaultdict
from importlib.metadata import version
import pandas as pd
from tabulate import tabulate


def display_results(ephemeris):
    """print out results to screen"""

    now = datetime.datetime.now()
    version = version("jwst_gtvt")
    welcome_string = "JWST General Target Visibility Tool"
    date_string = "Runtime/Date: {}".format(now)
    version_string = "Version Number: {}".format(version)

    print(tabulate([[welcome_string],
                    [date_string],
                    [version_string]], tablefmt='grid'))

    print('')

    pa_columns = {'NIRCAM_max_pa_angle': 'NRC Max PA', 
                  'NIRCAM_min_pa_angle': 'NRC Min PA',
                  'NIRSPEC_max_pa_angle': 'NRS Max PA', 
                  'NIRSPEC_min_pa_angle': 'NRS Min PA', 
                  'NIRISS_max_pa_angle': 'NIS Max PA',
                  'NIRISS_min_pa_angle': 'NIS Min PA', 
                  'MIRI_max_pa_angle': 'MIR Max PA', 
                  'MIRI_min_pa_angle': 'MIR Min PA',
                  'FGS_max_pa_angle': 'FGS Max PA', 
                  'FGS_min_pa_angle': 'FGS Min PA', 
                  'V3PA_max_pa_angle': 'V3 Max PA',
                  'V3PA_min_pa_angle': 'V3 Min PA'}

    df = ephemeris.dataframe

    if ephemeris.fixed is True:
        ra = max(df['ra'])
        dec = max(df['dec'])
        target_info_string = 'RA: %-*s  Dec: %-*s  Ecliptic Latitude: %s' % (10, ra, 10, dec, 24.284867)
    else:
        target_info_string = 'Target Name: %-*s Ecliptic Latitude: %s' % (10, ephemeris.target_name, 24.284867)

    df = df.loc[df['in_FOR']==True]

    interval_string = 'Checked Interval [{}, {}]'.format(ephemeris.start_date, ephemeris.end_date)

    print(target_info_string)
    print('-' * len(target_info_string))
    print('')
    print(interval_string)
    print('')

    window_indices = get_visibility_windows(df.index.tolist())
    window_dict = defaultdict(list)

    for start, end in window_indices:
        window_dict['Window Start'].append(ephemeris.dataframe.iloc[start]['Calendar Date (TDB)'])
        window_dict['Window End'].append(ephemeris.dataframe.iloc[end]['Calendar Date (TDB)'])
        window_dict['Window Duration'].append(end - start)
        window_dict['V3 Angle Start'].append(ephemeris.dataframe.iloc[start]['V3PA'])
        window_dict['V3 Angle End'].append(ephemeris.dataframe.iloc[end]['V3PA'])

    window_df = pd.DataFrame(window_dict)
    print(tabulate(window_df, headers='keys', tablefmt='psql', showindex=False))

    df = df[list(pa_columns.keys())]
    print(tabulate(df, headers=list(pa_columns.values()), tablefmt='psql', showindex=False))


def get_visibility_windows(df_indices):
    """Indices represent the spacing in positional measurements (1 day).
    This function captures the indices of visible windows for JWST observers.
    """

    window_indices = []
    window_start_index = df_indices[0]
    for x,y in zip(df_indices[::], df_indices[1::]):
        difference = y - x
        if difference > 1:
            window_end_index = x
            window_indices.append((window_start_index, window_end_index))
            window_start_index = y
        if y == df_indices[-1]:
            window_indices.append((window_start_index, y))
        else:
            continue

    return window_indices

