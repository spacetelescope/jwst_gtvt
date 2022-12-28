import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import numpy as np
from astropy.time import Time
import pandas as pd

from jwst_gtvt.ephemeride_rewrite import Ephemeris

def plot_single_instrument(ax, instrument_name, t, min_pa, max_pa):

    min_pa = np.array(min_pa)
    max_pa = np.array(max_pa)
    t = np.array(t)

    if np.any(min_pa > max_pa):
        minpa_lt_maxpa = min_pa < max_pa
        minpa_gt_maxpa = min_pa > max_pa

        max_pa_upper = np.copy(max_pa)
        min_pa_upper = np.copy(min_pa)
        max_pa_upper[minpa_gt_maxpa] = 360
        max_pa_upper[minpa_lt_maxpa] = np.nan
        min_pa_upper[minpa_lt_maxpa] = np.nan

        max_pa_lower = np.copy(max_pa)
        min_pa_lower = np.copy(min_pa)
        min_pa_lower[minpa_gt_maxpa] = 0
        max_pa_lower[minpa_lt_maxpa] = np.nan
        min_pa_lower[minpa_lt_maxpa] = np.nan


        max_pa[minpa_gt_maxpa] = np.nan
        min_pa[minpa_gt_maxpa] = np.nan

        ax.fill_between(t, min_pa_upper, max_pa_upper, facecolor='.7', edgecolor='.7', lw=2)
        ax.fill_between(t, min_pa_lower, max_pa_lower, facecolor='.7', edgecolor='.7', lw=2)
        ax.fill_between(t, min_pa, max_pa, edgecolor='.7', facecolor='.7', lw=2)
        ax.set_ylabel("Available Position Angle (Degree)")
        ax.set_title(instrument_name)
        ax.fmt_xdata = DateFormatter('%Y-%m-%d')


    else:
        ax.fill_between(t, min_pa, max_pa, edgecolor='none', facecolor='.7')
        ax.set_ylabel("Available Position Angle (Degree)")
        ax.set_title(instrument_name)
        ax.fmt_xdata = DateFormatter('%Y-%m-%d')

eph = Ephemeris()
eph.get_fixed_target_positions(253.2458, 2.4008)

df = eph.dataframe
df = df.loc[df['in_FOR']==True]

search_start = min(df['MJD'])
search_end = max(df['MJD'])

times = Time(df['MJD'], format='mjd').datetime

minV3PA_data = df['V3PA_min_pa_angle']
maxV3PA_data = df['V3PA_max_pa_angle']
minNIRCam_PA_data = df['NIRCAM_min_pa_angle']
maxNIRCam_PA_data = df['NIRCAM_max_pa_angle']
minMIRI_PA_data = df['MIRI_min_pa_angle']
maxMIRI_PA_data = df['MIRI_max_pa_angle']
minNIRSpec_PA_data = df['NIRSPEC_min_pa_angle']
maxNIRSpec_PA_data = df['NIRSPEC_max_pa_angle']
minNIRISS_PA_data = df['NIRISS_min_pa_angle']
maxNIRISS_PA_data = df['NIRISS_max_pa_angle']
minFGS_PA_data = df['FGS_min_pa_angle']
maxFGS_PA_data = df['FGS_max_pa_angle']

years = YearLocator()
months = MonthLocator()
yearsFmt = DateFormatter('%Y')
monthsFmt = DateFormatter('%m')
fig, axes = plt.subplots(2, 3, figsize=(14,8))

axes[0,0].set_title("V3")
plot_single_instrument(axes[0,0], "V3", times, minV3PA_data, maxV3PA_data)
axes[0,0].fmt_xdata = DateFormatter('%Y-%m-%d')
axes[0,0].set_ylabel("Available Position Angle (Degree)")
axes[0,0].set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)
labels = axes[0,0].get_xticklabels()
for label in labels:
    label.set_rotation(30)

plot_single_instrument(axes[0,1], 'NIRCam', times, minNIRCam_PA_data, maxNIRCam_PA_data)
axes[0,1].fmt_xdata = DateFormatter('%Y-%m-%d')
axes[0,1].set_ylabel("Available Position Angle (Degree)")
axes[0,1].set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)
labels = axes[0,1].get_xticklabels()
for label in labels:
    label.set_rotation(30)

axes[0,2].set_title("MIRI")
plot_single_instrument(axes[0,2], 'MIRI', times, minMIRI_PA_data, maxMIRI_PA_data)
axes[0,2].set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)
labels = axes[0,2].get_xticklabels()
for label in labels:
    label.set_rotation(30)

axes[1,0].set_title("NIRSpec")
axes[1,0].fmt_xdata = DateFormatter('%Y-%m-%d')
plot_single_instrument(axes[1,0], 'NIRSpec', times, minNIRSpec_PA_data, maxNIRSpec_PA_data)
axes[1,0].set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)
labels = axes[1,0].get_xticklabels()
for label in labels:
    label.set_rotation(30)

axes[1,1].set_title("NIRISS")
plot_single_instrument(axes[1,1], 'NIRISS', times, minNIRISS_PA_data, maxNIRISS_PA_data)
axes[1,1].set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)
labels = axes[1,1].get_xticklabels()
for label in labels:
    label.set_rotation(30)

axes[1,2].set_title("FGS")
plot_single_instrument(axes[1,2], 'FGS', times, minFGS_PA_data, maxFGS_PA_data)
axes[1,2].set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)
labels = axes[1,2].get_xticklabels()
for label in labels:
    label.set_rotation(30)
plt.show()
plt.close()
