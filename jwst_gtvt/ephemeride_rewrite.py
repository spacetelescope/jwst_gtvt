from astropy.time import Time
import pandas as pd
import pysiaf
import numpy as np

from .constants import UNIT_LIMIT

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
    def __init__(self, ephmeride_filename, start_date=Time.now(), end_date=None, verbose=True):
        """
        ephermeride_filename : str
            path to ephemeris file
        start_date : astropy.Time.time
            Start time of observering, default is current date
        end_date : astropy.Time.time
            End time of observing, default is end of mission lifetime
        verbose : bool
            Share table data
        """
        self.column_names = ['JDTDB','Calendar Date (TDB)','X','Y','Z','VX', 'VY', 'VZ']
        self.dataframe = pd.read_csv(ephmeride_filename, sep=",", names=self.column_names)
        self.dataframe = self.dataframe.drop(columns=['VX', 'VY', 'VZ'])  # We don't use the velocities in the GTVT/MTVT

        # Create MJD Column
        self.dataframe['MJD'] = self.dataframe['JDTDB'].values - 2400000.5

        self.start_date = int(start_date.mjd)

        if end_date:
            self.end_date = int(end_date.mjd)
        else:
            self.end_date = int(self.dataframe['MJD'].max())

        # only build dataframe based on start and end date and reset the index
        self.dataframe = self.dataframe[(self.dataframe['MJD'] >= self.start_date) &
                                        (self.dataframe['MJD'] <= self.end_date)].reset_index(drop=True)

        # Update positions based on pos() function
        for coordinate in ['X', 'Y', 'Z']:
            self.dataframe[coordinate] = (self.dataframe[coordinate].shift() - self.dataframe[coordinate]) * 1 + self.dataframe[coordinate]

    def build_entire_dataframe(self, ra, dec):
        """start with text file, end with data frame to generate the gtvt figures.
        """

        eph = Ephemeris('jwst_gtvt/data_only_ephemerides.txt')

        eph.dataframe = eph.sun_position_vectors(eph.dataframe)

        eph.dataframe = eph.sun_position_coordinates(eph.dataframe)

        eph.dataframe['ra'] = ra
        eph.dataframe['dec'] = dec

        eph.dataframe['ra_radians'] = eph.dataframe['ra'] * D2R
        eph.dataframe['dec_radians'] = eph.dataframe['dec'] * D2R

        eph.dataframe = eph.calculate_sun_pa(eph.dataframe,
                                             eph.dataframe['ra_radians'],
                                             eph.dataframe['dec_radians'],
                                             eph.dataframe['coord1'],
                                             eph.dataframe['coord2'])

        eph.dataframe = eph.normal_pa(eph.dataframe,
                                      eph.dataframe['ra_radians'],
                                      eph.dataframe['dec_radians'])

        eph.dataframe = eph.dist(eph.dataframe)
        eph.dataframe = eph.in_FOR(eph.dataframe)
        eph.dataframe = eph.get_allowed_max_boresight(eph.dataframe)

        instrument_aperture_pairs = [('NIRCAM', 'NRCALL_FULL'), ('NIRSPEC', 'NRS_FULL_MSA'),
                                     ('NIRISS', 'NIS_CEN'), ('MIRI', 'MIRIM_FULL'), 
                                     ('FGS', 'FGS1_FULL'), ('V3PA', None)]

        for instrument, aperture in instrument_aperture_pairs:
            eph.dataframe = self.calculate_min_max_pa_angles(eph.dataframe, instrument, aperture)

        return eph.dataframe

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

    def calculate_sun_pa(self, dataframe, tgt_coord1, tgt_coord2, sun_coord1, sun_coord2):
        """calculates position angle of object at tgt position."""

        sun_pa_y = np.cos(sun_coord2)*np.sin(sun_coord1 - tgt_coord1)
        sun_pa_x = (np.sin(sun_coord2)*np.cos(tgt_coord2)-np.cos(sun_coord2)*np.sin(tgt_coord2)*np.cos(sun_coord1 - tgt_coord1))

        dataframe['sun_pa'] = np.arctan2(sun_pa_y, sun_pa_x)

        dataframe.loc[dataframe["sun_pa"] < 0, 'sun_pa'] += 2*np.pi
        dataframe.loc[dataframe["sun_pa"] > 2*np.pi, 'sun_pa'] -= 2*np.pi

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

    def in_FOR(self, dataframe):

        dataframe['in_FOR'] = np.where((dataframe["dist"] < MAX_SUN_ANGLE) & (dataframe["dist"] > MIN_SUN_ANGLE), True, False)

        return dataframe

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