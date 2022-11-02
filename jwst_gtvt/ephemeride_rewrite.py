from astropy.time import Time
from astroquery.jplhorizons import Horizons
import csv
import pandas as pd
import pysiaf
import numpy as np
import os
import sys
import requests

from .constants import UNIT_LIMIT, URL

D2R = np.pi/180.  #degrees to radians
R2D = 180. / np.pi #radians to degrees 
unit_limit = lambda x: min(max(-1.,x),1.) # forces value to be in [-1,1]
MIN_SUN_ANGLE = 84.8 * D2R  #minimum Sun angle, in radians
MAX_SUN_ANGLE = 135.0 * D2R #maximum Sun angle, in radians
SUN_ANGLE_PAD = 0.5 * D2R   #pad away from Sun angle limits when constructing safe attitude

obliquity_of_the_ecliptic = -23.439291  # At J2000 equinox
obliquity_of_the_ecliptic *=  D2R
# Qecl2eci = QX(obliquity_of_the_ecliptic)


class Ephemeris:
    def __init__(self, start_date=Time('2021-12-26'), end_date=Time('2024-10-03')):
        """
        ephermeride_filename : str
            path to ephemeris file
        start_date : astropy.Time.time
            Start time of observering, default is current date
        end_date : astropy.Time.time
            End time of observing, default is end of mission lifetime
        verbose : bool
            Print jwst_gtvt results to screen
        """
        self.start_date = start_date
        self.end_date = end_date
        self.ephemeris_filename =  os.path.join(os.path.dirname(__file__), 'ephemeris_2021-12-26_2024-10-03.csv')  # set filename to local copy by default
        self.ephemeris = self.get_ephemeris_data(start_date.strftime('%Y-%m-%d'), end_date=end_date.strftime('%Y-%m-%d'))
        self.dataframe = self.convert_ephemeris_to_df(self.ephemeris)
        self.dataframe = self.dataframe.drop(columns=['VX', 'VY', 'VZ'])  # We don't use the velocities in the GTVT/MTVT

        # Create MJD Column
        self.dataframe['MJD'] = self.dataframe['JDTDB'].values - 2400000.5

        self.start_date_mjd = int(start_date.mjd)
        self.end_date_mjd = int(end_date.mjd)

        # only build dataframe based on start and end date and reset the index
        self.dataframe = self.dataframe[(self.dataframe['MJD'] >= self.start_date_mjd) &
                                        (self.dataframe['MJD'] <= self.end_date_mjd)].reset_index(drop=True)

        # Update positions based on pos() function
        for coordinate in ['X', 'Y', 'Z']:
            self.dataframe[coordinate] = (self.dataframe[coordinate].shift() - self.dataframe[coordinate]) * 1 + self.dataframe[coordinate]

    def get_allowed_max_boresight(self, dataframe):
        dataframe['max_boresight'] = dataframe.apply(lambda x: self.allowed_max_vehicle_roll(x['coord1'], x['coord2'], x['ra_radians'], x['dec_radians']), axis=1)
        dataframe['max_boresight'] *= R2D

        return dataframe

    def allowed_max_sun_roll(self, sun_p):

        abs_max_sun_roll = 5.2 * D2R

        if sun_p > 2.5 * D2R:
            max_sun_roll = abs_max_sun_roll - 1.7 * D2R * (sun_p - 2.5*D2R) / (5.2 - 2.5) / D2R
        else:
            max_sun_roll = abs_max_sun_roll

        max_sun_roll -= 0.1 * D2R  # Pad away from the edge

        return max_sun_roll

    def allowed_max_vehicle_roll(self, sun_ra, sun_dec, ra, dec):
        vehicle_pitch = np.pi/2. - self.angular_sep(sun_ra, sun_dec, ra, dec)
        sun_roll = 5.2 * D2R
        last_sun_roll = 0.
        while abs(sun_roll - last_sun_roll) > 0.0001*D2R:
            last_sun_roll = sun_roll
            sun_pitch = np.arcsin(UNIT_LIMIT(np.sin(vehicle_pitch)/np.cos(last_sun_roll)))
            sun_roll = self.allowed_max_sun_roll(sun_pitch)
            # print(last_sun_roll, sun_pitch, sun_roll)

        max_vehicle_roll = np.arcsin(UNIT_LIMIT(np.sin(sun_roll)/np.cos(vehicle_pitch)))

        return max_vehicle_roll

    def angular_sep(self, obj1_c1, obj1_c2, obj2_c1, obj2_c2):
        """angular distance between two objects, positions specified in spherical coordinates."""
        x = np.cos(obj2_c2) * np.cos(obj1_c2) * np.cos(obj2_c1 - obj1_c1) + np.sin(obj2_c2) * np.sin(obj1_c2)
        
        return np.arccos(UNIT_LIMIT(x))
    
    def calculate_min_max_pa_angles(self, dataframe, instrument, aperture=None, angle_name='V3IdlYAngle'):

        V3PA = dataframe['V3PA']
        max_boresight_roll = dataframe['max_boresight']

        if instrument == 'V3PA':
            minimum_angle = V3PA - max_boresight_roll
            maximum_angle = V3PA + max_boresight_roll
        else:
            instrument_ideal_v3_angle = self.get_angle(instrument, aperture, angle_name)
            minimum_angle = V3PA - max_boresight_roll + instrument_ideal_v3_angle
            maximum_angle = V3PA + max_boresight_roll + instrument_ideal_v3_angle

        dataframe[instrument + '_max_pa_angle'] = maximum_angle
        dataframe[instrument + '_min_pa_angle'] = minimum_angle

        dataframe.loc[dataframe[instrument + '_max_pa_angle'] < 0., instrument + '_max_pa_angle'] += 360
        dataframe.loc[dataframe[instrument + '_max_pa_angle'] > 360., instrument + '_max_pa_angle'] -= 360

        dataframe.loc[dataframe[instrument + '_min_pa_angle'] < 0., instrument + '_min_pa_angle'] += 360
        dataframe.loc[dataframe[instrument + '_min_pa_angle'] > 360., instrument + '_min_pa_angle'] -= 360

        return dataframe

    def calculate_sun_pa(self, dataframe, tgt_coord1, tgt_coord2, sun_coord1, sun_coord2):
        """calculates position angle of object at tgt position."""

        sun_pa_y = np.cos(sun_coord2)*np.sin(sun_coord1 - tgt_coord1)
        sun_pa_x = (np.sin(sun_coord2)*np.cos(tgt_coord2)-np.cos(sun_coord2)*np.sin(tgt_coord2)*np.cos(sun_coord1 - tgt_coord1))

        dataframe['sun_pa'] = np.arctan2(sun_pa_y, sun_pa_x)

        dataframe.loc[dataframe["sun_pa"] < 0, 'sun_pa'] += 2*np.pi
        dataframe.loc[dataframe["sun_pa"] > 2*np.pi, 'sun_pa'] -= 2*np.pi

        return dataframe

    def convert_ephemeris_to_df(self, ephemeris):

        start_index = np.where(ephemeris == '$$SOE')[0][0] + 1
        end_index = np.where(ephemeris == '$$EOE')[0][0]

        row_data = [row_data.split(',') for row_data in ephemeris[start_index:end_index]]
        result = [filter(None, row) for row in row_data]

        df = pd.DataFrame(result, columns=['JDTDB','Calendar Date (TDB)','X','Y','Z','VX', 'VY', 'VZ'])
        convert_dict = {'JDTDB': float,
                        'X': float,
                        'Y': float,
                        'Z': float}
        df = df.astype(convert_dict)

        return df

    def display_ephemeris_footer(self):
        start_index = np.where(self.ephemeris== '$$EOE')[0][0] + 1
        print(self.ephemeris[start_index:])

    def display_ephemeris_header(self):
        end_index = np.where(self.ephemeris == '$$SOE')[0][0] - 2
        print(self.ephemeris[0:end_index])

    def dist(self, dataframe):
        """angular distance betrween two objects, positions specified in spherical coordinates."""

        # obj1 = target
        # obj2 = sun

        obj1_c1, obj1_c2 = dataframe['ra_radians'], dataframe['dec_radians']
        obj2_c1, obj2_c2 = dataframe['coord1'], dataframe['coord2']
    
        data = np.cos(obj2_c2) * np.cos(obj1_c2) * np.cos(obj2_c1 - obj1_c1) + np.sin(obj2_c2) * np.sin(obj1_c2)

        anglar_distance = np.arccos([UNIT_LIMIT(value) for value in data])
        dataframe['dist'] = anglar_distance

        return dataframe

    def get_angle(self, instrument, aperture, angle_name):
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

    def get_ephemeris_data(self, start_date, end_date):
        """dates must be in format YYYY-MM-DD, returns text of ephemeris
        """ 
        try:
            url = URL.format(start_date, end_date)  # Get Horizons url for JWST ephemeris and add user specified dates
            self.eph_request = requests.get(url)
            ephemeris = np.array(self.eph_request.text.splitlines())
        except requests.exceptions.ConnectionError as err:
            print('No internet connection, using local file: {}'.format(self.ephemeris_filename))
            with open(self.ephemeris_filename, newline='') as csvfile:
                ephemeris = csv.reader(csvfile, delimiter=',')
            
        return ephemeris

    def get_moving_target_positions(self, desg, smallbody=False):
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
                       epochs={'start':self.start_date.to_value('iso', subfmt='date'), 
                               'stop':self.end_date.to_value('iso', subfmt='date'),
                               'step':'1d'})

        eph = obj.ephemerides(cache=False, quantities=(1))
        self.target_name = eph['targetname'][0]

        return eph['RA'].data.data, eph['DEC'].data.data

    def in_FOR(self, dataframe):
        dataframe['in_FOR'] = np.where((dataframe["dist"] < MAX_SUN_ANGLE) & (dataframe["dist"] > MIN_SUN_ANGLE), True, False)
        return dataframe

    def normal_pa(self, dataframe, tgt_ra, tgt_dec):
        """Calculate the normal position angle (V3)
        """

        if 'sun_pa' not in dataframe.columns:
            dataframe = self.sun_position_coordinates(dataframe)
            sun_coord1, sun_coord2 = dataframe['coord1'].values, dataframe['coord2'].values
            dataframe['sun_pa'] = self.calculate_sun_pa(dataframe, tgt_ra, tgt_dec, sun_coord1, sun_coord2)

        dataframe['V3PA'] = dataframe['sun_pa'] + np.pi  #  We want -V3 pointed towards sun.

        dataframe.loc[dataframe["V3PA"] < 0, 'V3PA'] += 2 * np.pi
        dataframe.loc[dataframe["V3PA"] > 2 * np.pi, 'V3PA'] -= 2 * np.pi

        dataframe['V3PA'] *= R2D

        return dataframe

    def retrieve_target_positions(self, ra, dec, fixed):
        """start with text file, end with data frame to generate the gtvt figures.
        """
        self.fixed = fixed
        self.dataframe = self.sun_position_vectors(self.dataframe)

        self.dataframe = self.sun_position_coordinates(self.dataframe)
     
        self.dataframe['ra'] = ra
        self.dataframe['dec'] = dec
     
        self.dataframe['ra_radians'] = self.dataframe['ra'] * D2R
        self.dataframe['dec_radians'] = self.dataframe['dec'] * D2R
     
        self.dataframe = self.calculate_sun_pa(self.dataframe,
                                               self.dataframe['ra_radians'],
                                               self.dataframe['dec_radians'],
                                               self.dataframe['coord1'],
                                               self.dataframe['coord2'])
     
        self.dataframe = self.normal_pa(self.dataframe,
                                        self.dataframe['ra_radians'],
                                        self.dataframe['dec_radians'])
     
        self.dataframe = self.dist(self.dataframe)
        self.dataframe = self.in_FOR(self.dataframe)
        self.dataframe = self.get_allowed_max_boresight(self.dataframe)

        instrument_aperture_pairs = [('NIRCAM', 'NRCALL_FULL'), ('NIRSPEC', 'NRS_FULL_MSA'),
                                     ('NIRISS', 'NIS_CEN'), ('MIRI', 'MIRIM_FULL'), 
                                     ('FGS', 'FGS1_FULL'), ('V3PA', None)]

        for instrument, aperture in instrument_aperture_pairs:
            self.dataframe = self.calculate_min_max_pa_angles(self.dataframe, instrument, aperture)

        return self.dataframe

    def sun_position_vectors(self, dataframe):
        """Add the sun vector to ephemeris data frame
        """

        # Create dataframe of only positional vectors.
        sun_vector_df = dataframe[['X', 'Y', 'Z']] * -1

        # Calculate vector magnitude for each row.
        dataframe['Vsun_mag'] = sun_vector_df.apply(np.linalg.norm, axis = 1)

        # Generate normal vectors.
        dataframe['Vsun_X'] = sun_vector_df['X']/dataframe['Vsun_mag']
        dataframe['Vsun_Y'] = sun_vector_df['Y']/dataframe['Vsun_mag']
        dataframe['Vsun_Z'] = sun_vector_df['Z']/dataframe['Vsun_mag']

        return dataframe

    def sun_position_coordinates(self, dataframe):
        """Get sun's position
        """

        dataframe['coord2'] = np.arcsin(dataframe['Vsun_Z'].values)
        dataframe['coord1'] = np.arctan2(dataframe['Vsun_Y'].values, dataframe['Vsun_X'].values)

        dataframe.loc[dataframe["coord1"] < 0, 'coord1'] += 2*np.pi

        return dataframe

    def write_ephemeris(self, write_path):
        outfname = write_path
        with open(outfname, "w") as of:
            for line in self.ephemeris:
                of.write(line)