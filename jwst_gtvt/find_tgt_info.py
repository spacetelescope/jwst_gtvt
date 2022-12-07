#! /usr/bin/env python

from __future__ import print_function

import sys
import math


import argparse
from astropy.time import Time
from astropy.table import Table
from astroquery.jplhorizons import Horizons

# Use TkAgg backend by default, but don't change backend if called from a Jupyter notebook with inline plots
# if 'module://ipykernel.pylab.backend_inline' not in matplotlib.rcParams['backend']:
#     matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import numpy as np
from os.path import join, abspath, dirname
import pysiaf
import warnings

from . import ephemeris_old2x as EPH

# ignore astropy warning that Date after 2020-12-30 is "dubious"
warnings.filterwarnings('ignore', category=UserWarning, append=True)
warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)

D2R = math.pi / 180.  #degrees to radians
R2D = 180. / math.pi #radians to degrees
PI2 = 2. * math.pi   # 2 pi
unit_limit = lambda x: min(max(-1.,x),1.) # forces value to be in [-1,1]

def convert_ddmmss_to_float(astring):
    aline = astring.split(':')
    d= float(aline[0])
    m= float(aline[1])
    s= float(aline[2])
    hour_or_deg = (s/60.+m)/60.+d
    return hour_or_deg

def bound_angle(ang):
    while ang<0.:
        ang += 360.
    while ang>360.:
        ang -= 360.
    return ang

def angular_sep(obj1_c1,obj1_c2,obj2_c1,obj2_c2):
    """angular distance betrween two objects, positions specified in spherical coordinates."""
    x = math.cos(obj2_c2)*math.cos(obj1_c2)*math.cos(obj2_c1-obj1_c1) + math.sin(obj2_c2)*math.sin(obj1_c2)
    return math.acos(unit_limit(x))

def calc_ecliptic_lat(ra, dec):
    NEP_ra = 270.000000 * D2R
    NEP_dec = 66.560708 * D2R
    a_sep = angular_sep(ra, dec, NEP_ra, NEP_dec)
    ecl_lat = math.pi/2. - a_sep
    return ecl_lat

def sun_pitch(aV):
    return math.atan2(aV.x,-aV.z)

def sun_roll(aV):
    return math.asin(-aV.y)

def allowed_max_sun_roll(sun_p):
    abs_max_sun_roll = 5.2 *D2R
    if sun_p > 2.5*D2R:
        max_sun_roll = abs_max_sun_roll - 1.7*D2R * (sun_p - 2.5*D2R)/(5.2 - 2.5)/D2R
    else:
        max_sun_roll = abs_max_sun_roll
    max_sun_roll -= 0.1*D2R  #Pad away from the edge
    return max_sun_roll

def allowed_max_vehicle_roll(sun_ra, sun_dec, ra, dec):
    vehicle_pitch = math.pi/2. - angular_sep(sun_ra, sun_dec, ra, dec)
    sun_roll = 5.2 * D2R
    last_sun_roll = 0.
    while abs(sun_roll - last_sun_roll) > 0.0001*D2R:
        last_sun_roll = sun_roll
        sun_pitch = math.asin(unit_limit(math.sin(vehicle_pitch)/math.cos(last_sun_roll)))
        sun_roll = allowed_max_sun_roll(sun_pitch)
        #print sun_roll*R2D,sun_pitch*R2D,vehicle_pitch*R2D
    max_vehicle_roll = math.asin(unit_limit(math.sin(sun_roll)/math.cos(vehicle_pitch)))
    return max_vehicle_roll

def get_target_ephemeris(desg, start_date, end_date, smallbody=False):
    """Ephemeris from JPL/HORIZONS.
    smallbody : bool, optional
      Set to `True` for comets and asteroids, `False` for planets,
      spacecraft, or moons.
    Returns : target name from HORIZONS, RA, and Dec.
    """

    if smallbody:
        bodytype='smallbody'
    else:
        bodytype='majorbody'

    obj = Horizons(id=desg, location='500@-170', id_type=bodytype,
                   epochs={'start':start_date, 'stop':end_date,
                   'step':'1d'})

    eph = obj.ephemerides(cache=False, quantities=(1))

    return eph['targetname'][0], eph['RA'], eph['DEC']


def window_summary_line(fixed, wstart, wend, pa_start, pa_end, ra_start, ra_end, dec_start, dec_end, cvz=False):
    """Formats window summary data for fixed and moving targets."""
    if cvz:
        line = " {0:15} {0:11} {0:11} ".format('CVZ')
    else:
        line = " {:15} {:11} {:11.2f} ".format(Time(wstart, format='mjd').datetime.date().strftime("%Y-%m-%d"),
                                               Time(wend, format='mjd').datetime.date().strftime("%Y-%m-%d"),wend-wstart)
    line += "{:13.5f} {:13.5f} ".format(pa_start*R2D,pa_end*R2D)
    if fixed:
        line += "{:13.5f} {:13.5f} ".format(ra_start*R2D, dec_start*R2D)
    else:
        line += "{:13.5f} {:13.5f} {:13.5f} {:13.5f} ".format(ra_start*R2D, ra_end*R2D, dec_start*R2D, dec_end*R2D)

    return line

def main(args, fixed=True):
    table_output=None
    if args.save_table is not None:
        table_output = open(args.save_table, 'w')

    NRCALL_FULL_V3IdlYang = get_angle('NIRCAM', 'NRCALL_FULL', 'V3IdlYAngle')
    NRS_FULL_MSA_V3IdlYang = get_angle('NIRSPEC', 'NRS_FULL_MSA', 'V3IdlYAngle')
    NIS_V3IdlYang = get_angle('NIRISS', 'NIS_CEN', 'V3IdlYAngle')
    MIRIM_FULL_V3IdlYang = get_angle('MIRI', 'MIRIM_FULL', 'V3IdlYAngle')
    FGS1_FULL_V3IdlYang = get_angle('FGS', 'FGS1_FULL', 'V3IdlYAngle')

    ECL_FLAG = False

    A_eph = EPH.Ephemeris(join(dirname(abspath(__file__)), "horizons_EM_jwst_wrt_sun_2021-2024.txt"),ECL_FLAG, verbose=args.no_verbose)

    search_start = Time(args.start_date, format='iso').mjd if args.start_date is not None else 59574.0  #Dec 26, 2021
    search_end = Time(args.end_date, format='iso').mjd if args.end_date is not None else 60585.0  #Oct 02, 2024

    if not (59574.0 <= search_start <= 60586.0) and args.start_date is not None:
        raise ValueError('Start date {} outside of available ephemeris {} to {}'.format(args.start_date, '2021-12-26', '2024-10-02'))
    if not (59574.0 <= search_end <= 60586.0) and args.end_date is not None:
        raise ValueError('End date {} outside of available ephemeris {} to {}'.format(args.end_date, '2021-12-26', '2024-10-02'))
    if search_start > search_end:
        raise ValueError('Start date {} should be before end date {}'.format(args.start_date, args.end_date))

    if search_start < A_eph.amin:
        print("Warning, search start time is earlier than ephemeris start.", file=table_output)
        search_start = A_eph.amin + 1

    scale = 1  # Channging this value must be reflected in get_target_ephemeris
    span = int(search_end-search_start)


    # if len(sys.argv) < 3:
    #   print "proper usage:"
    #   print "find_tgt_info.py ra dec [pa]"
    #   print "finds full visibility windows over [{}, {}]".format(Time(search_start, format='mjd').isot,
    #     Time(search_start+span/scale, format='mjd').isot)
    #   sys.exit(1)


    pa = 'X'
    if fixed:
        if args.ra.find(':')>-1:  #format is hh:mm:ss.s or  dd:mm:ss.s
          ra  = convert_ddmmss_to_float(args.ra) * 15. * D2R
          dec   = convert_ddmmss_to_float(args.dec) * D2R
        else: #format is decimal
          ra  = float(args.ra) * D2R
          dec   = float(args.dec) * D2R

        # although the coordinates are fixed, we need an array for
        # symmetry with moving target ephemerides
        ra = np.repeat(ra, span * scale + 1)
        dec = np.repeat(dec, span * scale + 1)
    else:
        assert len(args.ra) == span * scale + 1, "{} epochs retrieved for the moving target, but {} expected.".format(len(args.ra), span * scale + 1)
        ra = args.ra * D2R
        dec = args.dec * D2R

    if not args.no_verbose:
        print("", file=table_output)
        print("", file=table_output)
    if fixed:
        if not args.no_verbose:
            print("                ecliptic", file=table_output)
            print("RA      Dec     latitude", file=table_output)
            print("%7.3f %7.3f %7.3f" % (ra[0]*R2D,dec[0]*R2D,calc_ecliptic_lat(ra[0], dec[0])*R2D), file=table_output)

    if not fixed:
        if not args.no_verbose:
            print("TARGET NAME: {}".format(args.name))
            print("=========================================")


    if not args.no_verbose:
        print("", file=table_output)

    if args.v3pa is not None:
        pa     = float(args.v3pa) * D2R
    if not args.no_verbose:
        print("Checked interval [{}, {}]".format(Time(search_start, format='mjd').datetime.date().strftime("%Y-%m-%d"),
            Time(search_start+span, format='mjd').datetime.date().strftime("%Y-%m-%d")), file=table_output)
    if pa == "X":
        iflag_old = A_eph.in_FOR(search_start,ra[0],dec[0])
        if not args.no_verbose:
            print("|           Window [days]                 |    Normal V3 PA [deg]    |", end='', file=table_output)
    else:
        iflag_old = A_eph.is_valid(search_start,ra[0],dec[0],pa)
        if not args.no_verbose:
            print("|           Window [days]                 |   Specified V3 PA [deg]  |", end='', file=table_output)

    if fixed:
        if not args.no_verbose:
            print('\n', end='', file=table_output)
    else:
        if not args.no_verbose:
            print('{:^27s}|{:^27s}|'.format('RA', 'Dec'), file=table_output)

    if not args.no_verbose:
        print("   Start           End         Duration         Start         End    ", end='', file=table_output)
    if fixed:
        if not args.no_verbose:
            print("{:^13s} {:^13s}".format('RA', 'Dec'), file=table_output)
    else:
        if not args.no_verbose:
            print("{:^13s} {:^13s} {:^13s} {:^13s}".format('Start', 'End', 'Start', 'End'), file=table_output)

    if iflag_old:
      twstart = search_start
      ra_start = ra[0]
      dec_start = dec[0]
    else:
      twstart = -1.
    iflip = False

    #Step througth the interval and find where target goes in/out of field of regard.
    for i in range(1,span*scale+1):
        adate = search_start + float(i)/float(scale)
        #iflag = A_eph.in_FOR(adate,ra,dec)
        if pa == "X":
            iflag = A_eph.in_FOR(adate,ra[i],dec[i])
        else:
            iflag = A_eph.is_valid(adate,ra[i],dec[i],pa)
        if iflag != iflag_old:
            iflip = True
            if iflag:
                if pa == "X":
                    twstart = A_eph.bisect_by_FOR(adate,adate-0.1,ra[i],dec[i])
                else:
                    twstart = A_eph.bisect_by_attitude(adate,adate-0.1,ra[i],dec[i],pa)
                ra_start = ra[i]
                dec_start = dec[i]
            else:
                if pa == "X":
                    wend = A_eph.bisect_by_FOR(adate-0.1,adate,ra[i],dec[i])
                else:
                    wend = A_eph.bisect_by_attitude(adate-0.1,adate,ra[i],dec[i],pa)
                if twstart > 0.:
                    wstart = twstart #Only set wstart if wend is valid
                    if pa == "X":
                        pa_start = A_eph.normal_pa(wstart,ra_start,dec_start)
                        pa_end   = A_eph.normal_pa(wend,ra[i],dec[i])
                    else:
                        pa_start = pa
                        pa_end = pa
                    ra_end = ra[i]
                    dec_end = dec[i]
                    if not args.no_verbose:
                        print(window_summary_line(fixed, wstart, wend, pa_start, pa_end, ra_start, ra_end, dec_start, dec_end), file=table_output)
            iflag_old = iflag

    if iflip == True and iflag == True:
        if pa == "X":
            pa_start = A_eph.normal_pa(twstart,ra[i],dec[i])
            pa_end   = A_eph.normal_pa(adate,ra[i],dec[i])
        else:
            pa_start = pa
            pa_end = pa
        if not args.no_verbose:
            print(window_summary_line(fixed, twstart, adate, pa_start, pa_end, ra[i], ra[i], dec[i], dec[i]), file=table_output)

    if iflip == False and iflag == True and pa == "X":
        if dec[i] >0.:
            if not args.no_verbose:
                print(window_summary_line(fixed, 0, 0, 2 * np.pi, 0, ra[0], ra[-1], dec[0], dec[-1], cvz=True), file=table_output)
        else:
            if not args.no_verbose:
                print(window_summary_line(fixed, 0, 0, 0, 2 * np.pi, ra[0], ra[-1], dec[0], dec[-1], cvz=True), file=table_output)

    if 1==1:
        wstart = search_start
        wend = wstart + span
        istart = int(wstart)
        iend = int(wend)
        iflag = A_eph.in_FOR(wstart,ra[0],dec[0])
        tgt_is_in = False
        if iflag:
          tgt_is_in = True

        if not args.no_verbose:
            print("", file=table_output)
            print("", file=table_output)

        if fixed:
            fmt_repeats = 6
            if not args.no_verbose:
                print("                V3PA          NIRCam           NIRSpec         NIRISS           MIRI          FGS", file=table_output)
                print("   Date      min    max      min    max       min    max     min    max      min    max      min    max", file=table_output)
                    #58849.0 264.83 275.18 264.80 264.80  42.32  42.32 264.26 264.26 269.84 269.84 263.58 263.58
        else:
            fmt_repeats = 7
            if not args.no_verbose:
                print("                                V3PA          NIRCam           NIRSpec         NIRISS           MIRI          FGS", file=table_output)
                print("   Date      RA     Dec      min    max      min    max       min    max     min    max      min    max      min    max", file=table_output)
        
        times = []
        minV3PA_data = []
        maxV3PA_data = []
        minNIRCam_PA_data = []
        maxNIRCam_PA_data = []
        minNIRSpec_PA_data = []
        maxNIRSpec_PA_data = []
        minNIRISS_PA_data = []
        maxNIRISS_PA_data = []
        minMIRI_PA_data = []
        maxMIRI_PA_data = []
        minFGS_PA_data = []
        maxFGS_PA_data = []

        for itime in range(istart,iend):
            atime = float(itime)
            i = int((atime - search_start) * float(scale))
            iflag = A_eph.in_FOR(atime,ra[i],dec[i])
            #print atime,A_eph.in_FOR(atime,ra,dec)
            if iflag:
                if not tgt_is_in:
                    if not args.no_verbose:
                        print("", file=table_output)
                tgt_is_in = True

                V3PA = A_eph.normal_pa(atime,ra[i],dec[i])*R2D
                (sun_ra, sun_dec) = A_eph.sun_pos(atime)
                max_boresight_roll = allowed_max_vehicle_roll(sun_ra, sun_dec, ra[i], dec[i]) * R2D
                #sun_ang = angular_sep(sun_ra, sun_dec, ra, dec) * R2D

                minV3PA = bound_angle(V3PA - max_boresight_roll)
                maxV3PA = bound_angle(V3PA + max_boresight_roll)
                minNIRCam_PA = bound_angle(V3PA - max_boresight_roll + NRCALL_FULL_V3IdlYang)
                maxNIRCam_PA = bound_angle(V3PA + max_boresight_roll + NRCALL_FULL_V3IdlYang)
                minNIRSpec_PA = bound_angle(V3PA - max_boresight_roll + NRS_FULL_MSA_V3IdlYang)
                maxNIRSpec_PA = bound_angle(V3PA + max_boresight_roll + NRS_FULL_MSA_V3IdlYang)
                minNIRISS_PA = bound_angle(V3PA - max_boresight_roll + NIS_V3IdlYang)
                maxNIRISS_PA = bound_angle(V3PA + max_boresight_roll + NIS_V3IdlYang)
                minMIRI_PA = bound_angle(V3PA - max_boresight_roll + MIRIM_FULL_V3IdlYang)
                maxMIRI_PA = bound_angle(V3PA + max_boresight_roll + MIRIM_FULL_V3IdlYang)
                minFGS_PA = bound_angle(V3PA - max_boresight_roll + FGS1_FULL_V3IdlYang)
                maxFGS_PA = bound_angle(V3PA + max_boresight_roll + FGS1_FULL_V3IdlYang)

                times.append(Time(atime, format='mjd').datetime)
                minV3PA_data.append(minV3PA)
                maxV3PA_data.append(maxV3PA)
                minNIRCam_PA_data.append(minNIRCam_PA)
                maxNIRCam_PA_data.append(maxNIRCam_PA)
                minNIRSpec_PA_data.append(minNIRSpec_PA)
                maxNIRSpec_PA_data.append(maxNIRSpec_PA)
                minNIRISS_PA_data.append(minNIRISS_PA)
                maxNIRISS_PA_data.append(maxNIRISS_PA)
                minMIRI_PA_data.append(minMIRI_PA)
                maxMIRI_PA_data.append(maxMIRI_PA)
                minFGS_PA_data.append(minFGS_PA)
                maxFGS_PA_data.append(maxFGS_PA)
                #print '%7.1f %6.2f %6.2f %6.2f' % (atime, V3PA, NIRCam_PA, NIRSpec_PA)
                fmt = '{}' + '   {:6.2f} {:6.2f}'*fmt_repeats
                if fixed:
                    if not args.no_verbose:
                        print(fmt.format(
                            Time(atime, format='mjd').datetime.date().strftime("%Y-%m-%d"), minV3PA, maxV3PA,
                            minNIRCam_PA, maxNIRCam_PA, minNIRSpec_PA, maxNIRSpec_PA, minNIRISS_PA,
                            maxNIRISS_PA, minMIRI_PA, maxMIRI_PA, minFGS_PA, maxFGS_PA), file=table_output)#,sun_ang
                else:
                    if not args.no_verbose:
                        print(fmt.format(
                            Time(atime, format='mjd').datetime.date().strftime("%Y-%m-%d"), ra[i]*R2D, dec[i]*R2D, minV3PA, maxV3PA,
                            minNIRCam_PA, maxNIRCam_PA, minNIRSpec_PA, maxNIRSpec_PA, minNIRISS_PA,
                            maxNIRISS_PA, minMIRI_PA, maxMIRI_PA, minFGS_PA, maxFGS_PA), file=table_output)#,sun_ang
            else:
                tgt_is_in = False
                times.append(Time(atime, format='mjd').datetime)
                minV3PA_data.append(np.nan)
                maxV3PA_data.append(np.nan)
                minNIRCam_PA_data.append(np.nan)
                maxNIRCam_PA_data.append(np.nan)
                minNIRSpec_PA_data.append(np.nan)
                maxNIRSpec_PA_data.append(np.nan)
                minNIRISS_PA_data.append(np.nan)
                maxNIRISS_PA_data.append(np.nan)
                minMIRI_PA_data.append(np.nan)
                maxMIRI_PA_data.append(np.nan)
                minFGS_PA_data.append(np.nan)
                maxFGS_PA_data.append(np.nan)

        tab = Table([times, minV3PA_data, maxV3PA_data, minNIRCam_PA_data, maxNIRCam_PA_data,
            minNIRSpec_PA_data, maxNIRSpec_PA_data, minNIRISS_PA_data, maxNIRISS_PA_data,
            minMIRI_PA_data, maxMIRI_PA_data, minFGS_PA_data, maxFGS_PA_data],
            names=('Date', 'V3PA min', 'V3PA max', 'NIRCam min', 'NIRCam max',
                'NIRSpec min', 'NIRSpec max', 'NIRISS min', 'NIRISS max',
                'MIRI min', 'MIRI max', 'FGS min', 'FGS max'))

        # Plot observing windows
        if args.instrument is None:
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

            if fixed:
                axes[0,1].set_title('(R.A. = {}, Dec. = {})\n'.format(args.ra, args.dec)+"NIRCam")
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
            # fig.autofmt_xdate()

        elif args.instrument.lower() not in ['v3', 'nircam', 'miri', 'nirspec', 'niriss', 'fgs']:
            print()
            print(args.instrument+" not recognized. --instrument should be one of: v3, nircam, miri, nirspec, niriss, fgs")
            return

        elif args.instrument.lower() == 'v3':
            fig, ax = plt.subplots(figsize=(14,8))
            plot_single_instrument(ax, 'Observatory V3', times, minV3PA_data, maxV3PA_data)
            ax.set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)

        elif args.instrument.lower() == 'nircam':
            fig, ax = plt.subplots(figsize=(14,8))
            plot_single_instrument(ax, 'NIRCam', times, minNIRCam_PA_data, maxNIRCam_PA_data)
            ax.set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)

        elif args.instrument.lower() == 'miri':
            fig, ax = plt.subplots(figsize=(14,8))
            plot_single_instrument(ax, 'MIRI', times, minMIRI_PA_data, maxMIRI_PA_data)
            ax.set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)

        elif args.instrument.lower() == 'nirspec':
            fig, ax = plt.subplots(figsize=(14,8))
            plot_single_instrument(ax, 'NIRSpec', times, minNIRSpec_PA_data, maxNIRSpec_PA_data)
            ax.set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)

        elif args.instrument.lower() == 'niriss':
            fig, ax = plt.subplots(figsize=(14,8))
            plot_single_instrument(ax, 'NIRISS', times, minNIRISS_PA_data, maxNIRISS_PA_data)
            ax.set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)

        elif args.instrument.lower() == 'fgs':
            fig, ax = plt.subplots(figsize=(14,8))
            plot_single_instrument(ax, 'FGS', times, minFGS_PA_data, maxFGS_PA_data)
            ax.set_xlim(Time(search_start, format='mjd').datetime, Time(search_end, format='mjd').datetime)

        if args.name is not None:
            targname = args.name
        else:
            targname = ''
        if fixed:
            suptitle = '{} (RA = {}, DEC = {})'.format(targname, args.ra, args.dec)
        else:
            suptitle = '{}'.format(targname, args.ra, args.dec)
        fig.suptitle(suptitle, fontsize=18)
        fig.tight_layout()
        fig.subplots_adjust(top=0.88)

        if args.save_plot is None:
            plt.show()
        elif args.save_plot == 'test':
            plt.close()
        else:
            plt.savefig(args.save_plot)


def get_table(ra, dec, instrument=None, start_date=None, end_date=None, save_table=None, v3pa=None, fixed=True, verbose=True):
    """ Returns a table object with the PAs where the target is visible.

    parameters
    ----------
    ra : str
        Right Ascension
    dec : str
        Declination
    instrument : str
        Any of the JWST instruments
    start_date : float
    end_date : float
    save_table : str
        path of file to save plot output
    v3pa : str
        The position angle of the V3 axis.
    fixed : bool
        Whether or not the target is fixed. default = True


    returns
    -------
    astropy.table Table object
    """
    table_output=None
    if save_table is not None:
        table_output = open(save_table, 'w')

    NRCALL_FULL_V3IdlYang = get_angle('NIRCAM', 'NRCALL_FULL', 'V3IdlYAngle')
    NRS_FULL_MSA_V3IdlYang = get_angle('NIRSPEC', 'NRS_FULL_MSA', 'V3IdlYAngle')
    NIS_V3IdlYang = get_angle('NIRISS', 'NIS_CEN', 'V3IdlYAngle')
    MIRIM_FULL_V3IdlYang = get_angle('MIRI', 'MIRIM_FULL', 'V3IdlYAngle')
    FGS1_FULL_V3IdlYang = get_angle('FGS', 'FGS1_FULL', 'V3IdlYAngle')

    ECL_FLAG = False

    A_eph = EPH.Ephemeris(join(dirname(abspath(__file__)), "horizons_EM_jwst_wrt_sun_2021-2024.txt"),ECL_FLAG, verbose=verbose)

    search_start = Time(start_date, format='iso').mjd if start_date is not None else 59574.0  #Dec 26, 2021
    search_end = Time(end_date, format='iso').mjd if end_date is not None else 60585.0  #Oct 02, 2024

    if not (59574.0 <= search_start <= 60586.0) and start_date is not None:
        raise ValueError('Start date {} outside of available ephemeris {} to {}'.format(start_date, '2021-12-26', '2024-10-02'))
    if not (59574.0 <= search_end <= 60586.0) and end_date is not None:
        raise ValueError('End date {} outside of available ephemeris {} to {}'.format(end_date, '2021-12-26', '2024-10-02'))
    if search_start > search_end:
        raise ValueError('Start date {} should be before end date {}'.format(start_date, end_date))

    if search_start < A_eph.amin:
        print("Warning, search start time is earlier than ephemeris start.", file=table_output)
        search_start = A_eph.amin + 1

    scale = 1  # Channging this value must be reflected in get_target_ephemeris
    span = int(search_end-search_start)

    # if len(sys.argv) < 3:
    #   print "proper usage:"
    #   print "find_tgt_info.py ra dec [pa]"
    #   print "finds full visibility windows over [{}, {}]".format(Time(search_start, format='mjd').isot,
    #     Time(search_start+span/scale, format='mjd').isot)
    #   sys.exit(1)

    pa = 'X'

    if ra.find(':')>-1:  #format is hh:mm:ss.s or  dd:mm:ss.s
      ra  = convert_ddmmss_to_float(ra) * 15. * D2R
      dec   = convert_ddmmss_to_float(dec) * D2R
    else: #format is decimal
      ra  = float(ra) * D2R
      dec   = float(dec) * D2R

    # although the coordinates are fixed, we need an array for
    # symmetry with moving target ephemerides
    ra = np.repeat(ra, span * scale + 1)
    dec = np.repeat(dec, span * scale + 1)

    if verbose:
        print("", file=table_output)
        print("       Target", file=table_output)
    
    if verbose:
        if fixed:
            print("                ecliptic", file=table_output)
            print("RA      Dec     latitude", file=table_output)
            print("%7.3f %7.3f %7.3f" % (ra[0]*R2D,dec[0]*R2D,calc_ecliptic_lat(ra[0], dec[0])*R2D), file=table_output)
        print("", file=table_output)


    if v3pa is not None:
        pa     = float(v3pa) * D2R
    if verbose:
        print("Checked interval [{}, {}]".format(Time(search_start, format='mjd').datetime.date().strftime("%Y-%m-%d"),
            Time(search_start+span, format='mjd').datetime.date().strftime("%Y-%m-%d")), file=table_output)
    if pa == "X":
        iflag_old = A_eph.in_FOR(search_start,ra[0],dec[0])
        if verbose:
            print("|           Window [days]                 |    Normal V3 PA [deg]    |", end='', file=table_output)
    else:
        iflag_old = A_eph.is_valid(search_start,ra[0],dec[0],pa)
        if verbose:
            print("|           Window [days]                 |   Specified V3 PA [deg]  |", end='', file=table_output)

    if fixed:
        if verbose:
            print('\n', end='', file=table_output)
    else:
        if verbose:
            print('{:^27s}|{:^27s}|'.format('RA', 'Dec'), file=table_output)

    if verbose:
        print("   Start           End         Duration         Start         End    ", end='', file=table_output)
    if fixed:
        if verbose:
            print("{:^13s} {:^13s}".format('RA', 'Dec'), file=table_output)
    else:
        if verbose:
            print("{:^13s} {:^13s} {:^13s} {:^13s}".format('Start', 'End', 'Start', 'End'), file=table_output)

    if iflag_old:
      twstart = search_start
      ra_start = ra[0]
      dec_start = dec[0]
    else:
      twstart = -1.
    iflip = False

    #Step througth the interval and find where target goes in/out of field of regard.
    for i in range(1,span*scale+1):
        adate = search_start + float(i)/float(scale)
        #iflag = A_eph.in_FOR(adate,ra,dec)
        if pa == "X":
            iflag = A_eph.in_FOR(adate,ra[i],dec[i])
        else:
            iflag = A_eph.is_valid(adate,ra[i],dec[i],pa)
        if iflag != iflag_old:
            iflip = True
            if iflag:
                if pa == "X":
                    twstart = A_eph.bisect_by_FOR(adate,adate-0.1,ra[i],dec[i])
                else:
                    twstart = A_eph.bisect_by_attitude(adate,adate-0.1,ra[i],dec[i],pa)
                ra_start = ra[i]
                dec_start = dec[i]
            else:
                if pa == "X":
                    wend = A_eph.bisect_by_FOR(adate-0.1,adate,ra[i],dec[i])
                else:
                    wend = A_eph.bisect_by_attitude(adate-0.1,adate,ra[i],dec[i],pa)
                if twstart > 0.:
                    wstart = twstart #Only set wstart if wend is valid
                    if pa == "X":
                        pa_start = A_eph.normal_pa(wstart,ra_start,dec_start)
                        pa_end   = A_eph.normal_pa(wend,ra[i],dec[i])
                    else:
                        pa_start = pa
                        pa_end = pa
                    ra_end = ra[i]
                    dec_end = dec[i]
                    if verbose:
                        print(window_summary_line(fixed, wstart, wend, pa_start, pa_end, ra_start, ra_end, dec_start, dec_end), file=table_output)
            iflag_old = iflag

    if iflip == True and iflag == True:
        if pa == "X":
            pa_start = A_eph.normal_pa(twstart,ra[i],dec[i])
            pa_end   = A_eph.normal_pa(adate,ra[i],dec[i])
        else:
            pa_start = pa
            pa_end = pa
        if verbose:
            print(window_summary_line(fixed, twstart, adate, pa_start, pa_end, ra[i], ra[i], dec[i], dec[i]), file=table_output)

    if iflip == False and iflag == True and pa == "X":
        if dec[i] >0.:
            if verbose:
                print(window_summary_line(fixed, 0, 0, 2 * np.pi, 0, ra[0], ra[-1], dec[0], dec[-1], cvz=True), file=table_output)
        else:
            if verbose:
                print(window_summary_line(fixed, 0, 0, 0, 2 * np.pi, ra[0], ra[-1], dec[0], dec[-1], cvz=True), file=table_output)

    if 1==1:
        wstart = search_start
        wend = wstart + span
        istart = int(wstart)
        iend = int(wend)
        iflag = A_eph.in_FOR(wstart,ra[0],dec[0])
        tgt_is_in = False
        if iflag:
          tgt_is_in = True

        if verbose:
            print("", file=table_output)
            print("", file=table_output)
        if fixed:
            fmt_repeats = 6
            if verbose:
                print("                V3PA          NIRCam           NIRSpec         NIRISS           MIRI          FGS", file=table_output)
                print("   Date      min    max      min    max       min    max     min    max      min    max      min    max", file=table_output)
                    #58849.0 264.83 275.18 264.80 264.80  42.32  42.32 264.26 264.26 269.84 269.84 263.58 263.58
        else:
            fmt_repeats = 7
            if verbose:
                print("                                V3PA          NIRCam           NIRSpec         NIRISS           MIRI          FGS", file=table_output)
                print("   Date      RA     Dec      min    max      min    max       min    max     min    max      min    max      min    max", file=table_output)
        
        times = []
        V3PA_data = []
        minV3PA_data = []
        maxV3PA_data = []
        minNIRCam_PA_data = []
        maxNIRCam_PA_data = []
        minNIRSpec_PA_data = []
        maxNIRSpec_PA_data = []
        minNIRISS_PA_data = []
        maxNIRISS_PA_data = []
        minMIRI_PA_data = []
        maxMIRI_PA_data = []
        minFGS_PA_data = []
        maxFGS_PA_data = []

        nomNIRCam_PA_data = []
        nomNIRSpec_PA_data = []
        nomNIRISS_PA_data = []
        nomMIRI_PA_data = []
        nomFGS_PA_data = []

        for itime in range(istart,iend):
            atime = float(itime)
            i = int((atime - search_start) * float(scale))
            iflag = A_eph.in_FOR(atime,ra[i],dec[i])
            #print atime,A_eph.in_FOR(atime,ra,dec)
            if iflag:
                if not tgt_is_in:
                    if verbose:
                        print("", file=table_output)
                tgt_is_in = True

                V3PA = A_eph.normal_pa(atime,ra[i],dec[i])*R2D
                (sun_ra, sun_dec) = A_eph.sun_pos(atime)
                max_boresight_roll = allowed_max_vehicle_roll(sun_ra, sun_dec, ra[i], dec[i]) * R2D
                #sun_ang = angular_sep(sun_ra, sun_dec, ra, dec) * R2D

                minV3PA = bound_angle(V3PA - max_boresight_roll)
                maxV3PA = bound_angle(V3PA + max_boresight_roll)
                minNIRCam_PA = bound_angle(V3PA - max_boresight_roll + NRCALL_FULL_V3IdlYang)
                maxNIRCam_PA = bound_angle(V3PA + max_boresight_roll + NRCALL_FULL_V3IdlYang)
                minNIRSpec_PA = bound_angle(V3PA - max_boresight_roll + NRS_FULL_MSA_V3IdlYang)
                maxNIRSpec_PA = bound_angle(V3PA + max_boresight_roll + NRS_FULL_MSA_V3IdlYang)
                minNIRISS_PA = bound_angle(V3PA - max_boresight_roll + NIS_V3IdlYang)
                maxNIRISS_PA = bound_angle(V3PA + max_boresight_roll + NIS_V3IdlYang)
                minMIRI_PA = bound_angle(V3PA - max_boresight_roll + MIRIM_FULL_V3IdlYang)
                maxMIRI_PA = bound_angle(V3PA + max_boresight_roll + MIRIM_FULL_V3IdlYang)
                minFGS_PA = bound_angle(V3PA - max_boresight_roll + FGS1_FULL_V3IdlYang)
                maxFGS_PA = bound_angle(V3PA + max_boresight_roll + FGS1_FULL_V3IdlYang)

                nomNIRCam_PA = bound_angle(V3PA + NRCALL_FULL_V3IdlYang)
                nomNIRSpec_PA = bound_angle(V3PA + NRS_FULL_MSA_V3IdlYang)
                nomNIRISS_PA = bound_angle(V3PA + NIS_V3IdlYang)
                nomMIRI_PA = bound_angle(V3PA + MIRIM_FULL_V3IdlYang)
                nomFGS_PA = bound_angle(V3PA + FGS1_FULL_V3IdlYang)

                times.append(Time(atime, format='mjd').datetime)
                V3PA_data.append(V3PA)
                minV3PA_data.append(minV3PA)
                maxV3PA_data.append(maxV3PA)
                minNIRCam_PA_data.append(minNIRCam_PA)
                maxNIRCam_PA_data.append(maxNIRCam_PA)
                minNIRSpec_PA_data.append(minNIRSpec_PA)
                maxNIRSpec_PA_data.append(maxNIRSpec_PA)
                minNIRISS_PA_data.append(minNIRISS_PA)
                maxNIRISS_PA_data.append(maxNIRISS_PA)
                minMIRI_PA_data.append(minMIRI_PA)
                maxMIRI_PA_data.append(maxMIRI_PA)
                minFGS_PA_data.append(minFGS_PA)
                maxFGS_PA_data.append(maxFGS_PA)

                nomNIRCam_PA_data.append(nomNIRCam_PA)
                nomNIRSpec_PA_data.append(nomNIRSpec_PA)
                nomNIRISS_PA_data.append(nomNIRISS_PA)
                nomMIRI_PA_data.append(nomMIRI_PA)
                nomFGS_PA_data.append(nomFGS_PA)
                #print '%7.1f %6.2f %6.2f %6.2f' % (atime, V3PA, NIRCam_PA, NIRSpec_PA)
                fmt = '{}' + '   {:6.2f} {:6.2f}'*fmt_repeats
                if fixed:
                    if verbose:
                        print(fmt.format(
                            Time(atime, format='mjd').datetime.date().strftime("%Y-%m-%d"), V3PA, minV3PA, maxV3PA,
                            nomNIRCam_PA, minNIRCam_PA, maxNIRCam_PA, nomNIRSpec_PA, minNIRSpec_PA, maxNIRSpec_PA,
                            nomNIRISS_PA, minNIRISS_PA, maxNIRISS_PA, nomMIRI_PA, minMIRI_PA, maxMIRI_PA,
                            nomFGS_PA, minFGS_PA, maxFGS_PA), file=table_output)#,sun_ang
                else:
                    if verbose:
                        print(fmt.format(
                            Time(atime, format='mjd').datetime.date().strftime("%Y-%m-%d"), ra[i]*R2D, dec[i]*R2D, V3PA, minV3PA, maxV3PA,
                            nomNIRCam_PA, minNIRCam_PA, maxNIRCam_PA, nomNIRSpec_PA, minNIRSpec_PA, maxNIRSpec_PA,
                            nomNIRISS_PA, minNIRISS_PA, maxNIRISS_PA, nomMIRI_PA, minMIRI_PA, maxMIRI_PA,
                            nomFGS_PA, minFGS_PA, maxFGS_PA), file=table_output)#,sun_ang
            else:
                tgt_is_in = False
                times.append(Time(atime, format='mjd').datetime)
                V3PA_data.append(np.nan)
                minV3PA_data.append(np.nan)
                maxV3PA_data.append(np.nan)
                minNIRCam_PA_data.append(np.nan)
                maxNIRCam_PA_data.append(np.nan)
                minNIRSpec_PA_data.append(np.nan)
                maxNIRSpec_PA_data.append(np.nan)
                minNIRISS_PA_data.append(np.nan)
                maxNIRISS_PA_data.append(np.nan)
                minMIRI_PA_data.append(np.nan)
                maxMIRI_PA_data.append(np.nan)
                minFGS_PA_data.append(np.nan)
                maxFGS_PA_data.append(np.nan)

                nomNIRCam_PA_data.append(np.nan)
                nomNIRSpec_PA_data.append(np.nan)
                nomNIRISS_PA_data.append(np.nan)
                nomMIRI_PA_data.append(np.nan)
                nomFGS_PA_data.append(np.nan)


        tab = Table([times, V3PA_data, minV3PA_data, maxV3PA_data, nomNIRCam_PA_data, minNIRCam_PA_data, maxNIRCam_PA_data,
            nomNIRSpec_PA_data, minNIRSpec_PA_data, maxNIRSpec_PA_data, nomNIRISS_PA_data, minNIRISS_PA_data, maxNIRISS_PA_data,
            nomMIRI_PA_data, minMIRI_PA_data, maxMIRI_PA_data, nomFGS_PA_data, minFGS_PA_data, maxFGS_PA_data],
            names=('Date', 'V3PA', 'V3PA min', 'V3PA max', 'NIRCam nom', 'NIRCam min', 'NIRCam max',
                'NIRSpec nom', 'NIRSpec min', 'NIRSpec max', 'NIRISS nom', 'NIRISS min', 'NIRISS max',
                'MIRI nom', 'MIRI min', 'MIRI max', 'FGS nom', 'FGS min', 'FGS max'))

    return tab





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


def get_angle(instrument, aperture, angle_name):
    """Get angle requested by user

    Parameters
    ----------
    instrument : JWST instrument of interest
        type : str
    
    aperture : instrument observing aperture
        type : str

    angle_name : angle of interest
        type : str

    Returns
    -------
    angle : the angle obtained from the SIAF.
        type : float
    """
    
    siaf = pysiaf.Siaf(instrument)
    meta = siaf[aperture]
    angle = getattr(meta, angle_name)

    return angle


if __name__ == '__main__':
    try:
        # see if there is a negative dec in sexagesimal coordinates
        dec_index = [':' in arg and arg.startswith('-') for arg in sys.argv].index(True)
        arg_list = sys.argv[1:]
        dec = arg_list.pop(dec_index-1)
        arg_list.append('--')
        arg_list.append(dec)

    except ValueError:
        arg_list = sys.argv[1:]

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('ra', help='Right Ascension of target in either sexagesimal (hh:mm:ss.s) or degrees')
    parser.add_argument('dec', help='Declination of target in either sexagesimal (dd:mm:ss.s) or degrees')
    parser.add_argument('--pa', help='Specify a desired Position Angle')
    parser.add_argument('--save_plot', help='Path of file to save plot output')
    parser.add_argument('--save_table', help='Path of file to save table output')
    parser.add_argument('--instrument', help='If specified plot shows only windows for this instrument')
    parser.add_argument('--name', help='Target Name to appear on plots')
    parser.add_argument('--start_date', help='Start date for visibility search in yyyy-mm-dd format')
    parser.add_argument('--end_date', help='End date for visibility search in yyyy-mm-dd format')
    parser.add_argument('--no_verbose', help='Suppress table output to screen', type=bool)
    args = parser.parse_args(arg_list)

    main(args)
