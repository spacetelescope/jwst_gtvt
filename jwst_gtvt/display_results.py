import datetime
from tabulate import tabulate


def display_results(ephemeris):
    """print out results to screen"""

    now = datetime.datetime.now()

    welcome_string = """JWST General Target Visibility Tool
                        {}
                     """.format(now)

    print(tabulate([[welcome_string]], tablefmt='grid'))

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

    ra = max(df['ra'])
    dec = max(df['dec'])

    df = df.loc[df['in_FOR']==True]

    ra_dec_string = 'RA: %-*s  Dec: %-*s  Ecliptic Latitude: %s' % (10, ra, 10, dec, 24.284867)
    interval_string = 'Checked Interval [{}, {}]'.format(ephemeris.start_date, ephemeris.end_date)
    print(ra_dec_string)
    print('-' * len(ra_dec_string))
    print(interval_string)
    print('-' * len(interval_string))

    window_indices = get_visibility_windows(df.index.tolist())

    for start, end in window_indices:
        start_date = ephemeris.dataframe.iloc[start]['Calendar Date (TDB)']
        end_date = ephemeris.dataframe.iloc[end]['Calendar Date (TDB)']
        window_duration = end - start
        start_v3 = ephemeris.dataframe.iloc[start]['V3PA']
        end_v3 = ephemeris.dataframe.iloc[end]['V3PA']

        print('%-*s  %-*s  %-*s  %-*s  %-*s' % (10, start_date, 10, end_date, 
                                                10, window_duration,10, start_v3,
                                                10, end_v3))

    df = df[list(pa_columns.keys())]
    print(tabulate(df, headers=list(pa_columns.values()), tablefmt='psql'))


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