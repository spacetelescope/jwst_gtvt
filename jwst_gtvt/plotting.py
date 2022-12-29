import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from astropy.time import Time

from jwst_gtvt.display_results import get_visibility_windows 

def plot_visibility(ephemeris, instrument=None):
    # Just incase dataframe hasn't been sorted yet
    dataframe = ephemeris.dataframe
    df = dataframe.loc[dataframe['in_FOR']==True]

    times = Time(df['MJD'], format='mjd').datetime

    # These indices allow us to get the visible regions regions
    window_indices = get_visibility_windows(df.index.tolist())

    if instrument:
        # Set plotting configs
        plt.figure(figsize=(18, 8))
        plt.grid(color='k', linestyle='--', linewidth=2, alpha=0.3)
        plt.xticks(fontsize=14, rotation=45)
        plt.yticks(fontsize=14)

        for start, end in window_indices:
            data_to_plot = df.loc[start:end]
            min_PA_data = data_to_plot[instrument.upper() + '_min_pa_angle']
            max_PA_data = data_to_plot[instrument.upper() + '_max_pa_angle']
            plt.fill_between(data_to_plot['MJD'], min_PA_data, max_PA_data, color='grey')

        plt.ylabel('Available Position Angles (Degrees)', fontsize=20)

        if ephemeris.fixed:
            ra, dec = max(df['ra']), max(df['dec'])
            plt.title('RA: {} Dec: {} with {}'.format(ra, dec, instrument.upper()), fontsize=25)
        else:
            plt.title('Target {} with {}'.format(ephemeris.target_name, instrument.upper()), fontsize=25)
        plt.show()

    else:
        instrument_names = ['NIRCAM', 'NIRSPEC', 'NIRISS', 'MIRI', 'FGS', 'V3PA']
        fig, axs = plt.subplots(ncols=3, nrows=2, figsize = (18,8))
        fig.tight_layout()

        if ephemeris.fixed:
            ra, dec = max(df['ra']), max(df['dec'])
            fig.suptitle('RA: {} Dec: {}'.format(ra, dec), fontsize=25)
        else:
            plt.suptitle('Target {}'.format(ephemeris.target_name), fontsize=25)

        for instrument_name, ax in zip(instrument_names, axs.flatten()):
            for start, end in window_indices:
                data_to_plot = df.loc[start:end]
                min_PA_data = data_to_plot[instrument_name + '_min_pa_angle']
                max_PA_data = data_to_plot[instrument_name + '_max_pa_angle']
                ax.fill_between(data_to_plot['MJD'], min_PA_data, max_PA_data, color='grey')
                ax.set_title(instrument_name)
                ax.tick_params('x', labelrotation=45)
                ax.grid(color='k', linestyle='--', linewidth=2, alpha=0.3)
                ax.set_ylabel('Available Position Angles (Degrees)')
        plt.show()

