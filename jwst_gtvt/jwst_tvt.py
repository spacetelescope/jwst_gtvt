#! /usr/bin/env python

"""Generate visibility plots for GTVT/MTVT

This script contains the main components that generate the dataframe
containing information displayed in the figured created by the jwst_gtvt.


Authors
-------
    - Mees Fix
    - Bryan Holler
    - Melanie Clarke
    - Matt Hill

Use
---
    This script is intended to be executed as such:
    ::
        >>> from jwst_gtvt.ephemeride_rewrite import Ephemeris
        >>> eph = Ephemeris()
"""

from astropy.time import Time
import astropy.units as u
from astroquery.jplhorizons import Horizons
import datetime
import glob
import pandas as pd
import pysiaf
import numpy as np
import os
import re
import requests

from jwst_gtvt.constants import UNIT_LIMIT, URL

D2R = np.pi / 180.0  # degrees to radians
R2D = 180.0 / np.pi  # radians to degrees
unit_limit = lambda x: min(max(-1.0, x), 1.0)  # forces value to be in [-1,1]
LAUNCH_DATE = "2021-12-26"
MIN_SUN_ANGLE = 84.8 * D2R  # minimum Sun angle, in radians
MAX_SUN_ANGLE = 135.0 * D2R  # maximum Sun angle, in radians
SUN_ANGLE_PAD = (
    0.5 * D2R
)  # pad away from Sun angle limits when constructing safe attitude

obliquity_of_the_ecliptic = -23.439291  # At J2000 equinox
obliquity_of_the_ecliptic *= D2R
# Qecl2eci = QX(obliquity_of_the_ecliptic)

NOW = Time.now()

class Ephemeris:
    def __init__(self, start_date=NOW, end_date=NOW + 2 * u.year):
        """
        ephermeride_filename : str
            path to ephemeris file
        start_date : astropy.Time.time
            Start time of observering, default is current date
        end_date : astropy.Time.time
            End time of observing, default is current date plus two years
        verbose : bool
            Print jwst_gtvt results to screen
        """

        # Using a code snippet writeen by Melanie Clarke to check the max date of the ephemeris.
        # On March 20th, 2025, the furthest projected date is 2030-03-16.
        # This date derived from the intentially failed ephemeris call relies on the first line of the
        # returned call containing a message about what the furthest projected date is.
        # IF there are changes to HORIZONS ephemerides structure, this code could potentially fail.
        self.max_date = self.ephemeris_maximum_date()

        if start_date < Time(LAUNCH_DATE) or end_date > Time(self.max_date):
            date_out_of_bound_msg = (
                "Time frame selected {} ----> {} is out of bounds!".format(
                    start_date, end_date
                ),
                "Please select dates between {} ----> {}".format(
                    LAUNCH_DATE, self.max_date
                ),
            )
            raise SystemExit(date_out_of_bound_msg)
        elif Time(start_date) > Time(end_date):
            raise SystemExit("start_date later than end_date, exiting")

        else:
            self.start_date = start_date
            self.end_date = end_date

            project_dirname = os.path.dirname(__file__)
            filename_str = "ephemeris_????-??-??_????-??-??.txt"
            file_search_result = glob.glob(
                os.path.join(project_dirname, "data", filename_str)
            )[0]

            self.ephemeris_filename = (
                file_search_result  # set filename to local copy by default
            )
            self.ephemeris = self.get_ephemeris_data(
                start_date=start_date.strftime("%Y-%m-%d"),
                end_date=end_date.strftime("%Y-%m-%d"),
            )
            self.dataframe = self.convert_ephemeris_to_df(self.ephemeris)
            self.dataframe = self.dataframe.drop(
                columns=["VX", "VY", "VZ"]
            )  # We don't use the velocities in the GTVT/MTVT

            # Create MJD Column
            self.dataframe["MJD"] = self.dataframe["JDTDB"].values - 2400000.5

            self.start_date_mjd = int(start_date.mjd)
            self.end_date_mjd = int(end_date.mjd)

            # only build dataframe based on start and end date and reset the index
            self.dataframe = self.dataframe[
                (self.dataframe["MJD"] >= self.start_date_mjd)
                & (self.dataframe["MJD"] <= self.end_date_mjd)
            ].reset_index(drop=True)

            # Update positions based on pos() function
            for coordinate in ["X", "Y", "Z"]:
                self.dataframe[coordinate] = (
                    self.dataframe[coordinate].shift() - self.dataframe[coordinate]
                ) * 1 + self.dataframe[coordinate]

    def convert_ddmmss_to_float(self, astring):
        """Convert date ra dec to sexigesimal"""
        aline = astring.split(":")
        d = float(aline[0])
        m = float(aline[1])
        s = float(aline[2])

        sign = np.sign(d)

        if sign != 0.0:
            hour_or_deg = sign * ((s / 60.0 + m) / 60.0 + np.abs(d))
        else:
            hour_or_deg = (s / 60.0 + m) / 60.0 + d

        return hour_or_deg

    def get_allowed_max_boresight(self, dataframe):
        """Need Docstring"""
        dataframe["max_boresight"] = dataframe.apply(
            lambda x: self.allowed_max_vehicle_roll(
                x["coord1"], x["coord2"], x["ra_radians"], x["dec_radians"]
            ),
            axis=1,
        )
        dataframe["max_boresight"] *= R2D

        return dataframe

    def allowed_max_sun_roll(self, sun_p):
        """Need Docstring"""
        abs_max_sun_roll = 5.2 * D2R

        if sun_p > 2.5 * D2R:
            max_sun_roll = (
                abs_max_sun_roll - 1.7 * D2R * (sun_p - 2.5 * D2R) / (5.2 - 2.5) / D2R
            )
        else:
            max_sun_roll = abs_max_sun_roll

        max_sun_roll -= 0.1 * D2R  # Pad away from the edge

        return max_sun_roll

    def allowed_max_vehicle_roll(self, sun_ra, sun_dec, ra, dec):
        """Need Docstring

        Parameters
        ----------
        sun_ra : float
            Sun's right ascension
        sun_dec : float
            Sun's declination
        ra : float
            Target's right ascension
        dec : float
            Target's declination

        Returns
        -------
        max_vehicle_roll : float
            Maximum allowed roll of spacecraft.
        """
        vehicle_pitch = np.pi / 2.0 - self.angular_sep(sun_ra, sun_dec, ra, dec)
        sun_roll = 5.2 * D2R
        last_sun_roll = 0.0
        while abs(sun_roll - last_sun_roll) > 0.0001 * D2R:
            last_sun_roll = sun_roll
            sun_pitch = np.arcsin(
                UNIT_LIMIT(np.sin(vehicle_pitch) / np.cos(last_sun_roll))
            )
            sun_roll = self.allowed_max_sun_roll(sun_pitch)
            # print(last_sun_roll, sun_pitch, sun_roll)

        max_vehicle_roll = np.arcsin(
            UNIT_LIMIT(np.sin(sun_roll) / np.cos(vehicle_pitch))
        )

        return max_vehicle_roll

    def angular_sep(self, obj1_c1, obj1_c2, obj2_c1, obj2_c2):
        """angular distance between two objects, positions specified in spherical coordinates.

        Parameters
        ----------
        obj1_c1 : float
            Object one coordinate 1
        obj2_c2 : float
            Object one coordinate 2
        obj2_c1 : float
            Object two coordinate 1
        obj2_c2 : float
            Object two coordinate 2

        Returns
        -------
        np.arccos(UNIT_LIMIT(x)) : float
            Angular seperation of two targets
        """

        x = np.cos(obj2_c2) * np.cos(obj1_c2) * np.cos(obj2_c1 - obj1_c1) + np.sin(
            obj2_c2
        ) * np.sin(obj1_c2)

        return np.arccos(UNIT_LIMIT(x))

    def calculate_min_max_pa_angles(
        self, dataframe, instrument, aperture=None, angle_name="V3IdlYAngle"
    ):
        """Get the minimum and maximum position angle

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Pandas dataframe containing metadata from tool
        instrument : str
            JWST instrument name
        aperture : str
            JWST instrument aperture name
        angle_name : str
            pysiaf angle name [Default : V3IdlYAngle]

        Return
        ------
        dataframe : pd.DataFrame
            Pandas dataframe with updated/new data
        """

        V3PA = dataframe["V3PA"]
        max_boresight_roll = dataframe["max_boresight"]

        if instrument == "V3PA":
            nominal_angle = V3PA
            minimum_angle = V3PA - max_boresight_roll
            maximum_angle = V3PA + max_boresight_roll
        else:
            instrument_ideal_v3_angle = self.get_angle(instrument, aperture, angle_name)
            nominal_angle = V3PA + instrument_ideal_v3_angle
            minimum_angle = V3PA - max_boresight_roll + instrument_ideal_v3_angle
            maximum_angle = V3PA + max_boresight_roll + instrument_ideal_v3_angle

        dataframe[instrument + "_nominal_angle"] = nominal_angle
        dataframe[instrument + "_max_pa_angle"] = maximum_angle
        dataframe[instrument + "_min_pa_angle"] = minimum_angle

        dataframe.loc[
            dataframe[instrument + "_max_pa_angle"] < 0.0, instrument + "_max_pa_angle"
        ] += 360
        dataframe.loc[
            dataframe[instrument + "_max_pa_angle"] > 360.0,
            instrument + "_max_pa_angle",
        ] -= 360

        dataframe.loc[
            dataframe[instrument + "_min_pa_angle"] < 0.0, instrument + "_min_pa_angle"
        ] += 360
        dataframe.loc[
            dataframe[instrument + "_min_pa_angle"] > 360.0,
            instrument + "_min_pa_angle",
        ] -= 360

        dataframe.loc[
            dataframe[instrument + "_nominal_angle"] < 0.0,
            instrument + "_nominal_angle",
        ] += 360
        dataframe.loc[
            dataframe[instrument + "_nominal_angle"] > 360.0,
            instrument + "_nominal_angle",
        ] -= 360

        return dataframe

    def calculate_sun_pa(
        self, dataframe, tgt_coord1, tgt_coord2, sun_coord1, sun_coord2
    ):
        """calculates position angle of object at tgt position.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Pandas dataframe
        tgt_coord1 : float
            Target coordinate 1
        tgt_coord2 : float
            Target coordinate 2
        sun_coord1 : float
            Sun's coordinate 1
        sun_coord2 : float
            Sun's coordinate 2

        Return
        ------
        dataframe : pandas.dataframe
            Pandas dateframe with updated data
        """

        sun_pa_y = np.cos(sun_coord2) * np.sin(sun_coord1 - tgt_coord1)
        sun_pa_x = np.sin(sun_coord2) * np.cos(tgt_coord2) - np.cos(
            sun_coord2
        ) * np.sin(tgt_coord2) * np.cos(sun_coord1 - tgt_coord1)

        dataframe["sun_pa"] = np.arctan2(sun_pa_y, sun_pa_x)

        dataframe.loc[dataframe["sun_pa"] < 0, "sun_pa"] += 2 * np.pi
        dataframe.loc[dataframe["sun_pa"] > 2 * np.pi, "sun_pa"] -= 2 * np.pi

        return dataframe

    def convert_ephemeris_to_df(self, ephemeris):
        """Convert ephemeris data into a dataframe

        Parameters
        ----------
        ephemeris : str
            Ephemeris as python string

        Returns
        -------
        df : pandas.DataFrame
            JWST ephmeris as pandas dataframe
        """
        try:
            start_index = np.where(ephemeris == "$$SOE")[0][0] + 1
            end_index = np.where(ephemeris == "$$EOE")[0][0]
        except IndexError:
            # No positions returned
            idx_err_msg = (
                "No position angles in field of regard! "
                "Check constraints for your target and if it is observable with JWST. \n"
                "Vist: https://jwst-docs.stsci.edu/jwst-observatory-characteristics/jwst-observatory-coordinate-system-and-field-of-regard for more information"
            )
            raise IndexError(idx_err_msg)

        row_data = [
            row_data.split(",") for row_data in ephemeris[start_index:end_index]
        ]
        result = [filter(None, row) for row in row_data]

        df = pd.DataFrame(
            result,
            columns=["JDTDB", "Calendar Date (TDB)", "X", "Y", "Z", "VX", "VY", "VZ"],
        )
        convert_dict = {"JDTDB": float, "X": float, "Y": float, "Z": float}
        df = df.astype(convert_dict)

        return df

    def display_ephemeris_footer(self):
        """Show metadata of ephemeris footer."""
        start_index = np.where(self.ephemeris == "$$EOE")[0][0] + 1
        print(self.ephemeris[start_index:])

    def display_ephemeris_header(self):
        """Show metadata of ephemeris header."""
        end_index = np.where(self.ephemeris == "$$SOE")[0][0] - 2
        print(self.ephemeris[0:end_index])

    def dist(self, dataframe):
        """angular distance betrween two objects, positions specified in spherical coordinates.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Pandas dataframe with updated metadata
        """

        # obj1 = target
        # obj2 = sun

        obj1_c1, obj1_c2 = dataframe["ra_radians"], dataframe["dec_radians"]
        obj2_c1, obj2_c2 = dataframe["coord1"], dataframe["coord2"]

        data = np.cos(obj2_c2) * np.cos(obj1_c2) * np.cos(obj2_c1 - obj1_c1) + np.sin(
            obj2_c2
        ) * np.sin(obj1_c2)

        anglar_distance = np.arccos([UNIT_LIMIT(value) for value in data])
        dataframe["dist"] = anglar_distance

        return dataframe

    def ephemeris_maximum_date(self):
        """Retrieve the last available date for JWST ephemerides."""

        # attempt to retrieve an ephemeris for a date too far in the future
        future_date = "9999-01-01"
        request_url = URL.format(LAUNCH_DATE, future_date)

        try:
            # this should return an error message containing the last good date
            ephemeris_request = requests.get(request_url)

            # parse the date from the message
            m = re.match(
                r".*after A\.D\. (\d{4}-[a-zA-Z]+-\d{1,2})", ephemeris_request.text
            )
            dt = m.groups()[0]
            end_date = datetime.datetime.strptime(dt, "%Y-%b-%d").strftime("%Y-%m-%d")
        except Exception:
            # if the above fails for any reason, fall back to a known good date
            end_date = "2030-03-16"

        return end_date

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

    def get_ephemeris_data(self, start_date=None, end_date=None):
        """Read JWST data and make python object.

        Parameters
        ----------
        start_date : astropy.Time.strftime
            Start date of ephemeris (YYYY-MM-DD)
        end_date : astropy.Time.strftime
            End date of ephemeris (YYYY-MM-DD)

        Returns
        -------
        ephemeris : np.array
            Ephemeris in an array split by lines.
        """

        try:
            self.url = URL.format(
                start_date, end_date
            )  # Get Horizons url for JWST ephemeris and add user specified dates
            self.eph_request = requests.get(self.url)
            ephemeris = np.array(self.eph_request.text.splitlines())
        except Exception as e:
            print(
                "Issue reading ephemeris from HORIZONS. Using local ephemeris {} Full traceback below:".format(
                    self.ephemeris_filename
                )
            )
            print(e)
            with open(self.ephemeris_filename) as f:
                lines = np.array(f.read().splitlines())
            ephemeris = np.array(lines)

        return ephemeris

    def get_fixed_target_positions(self, ra, dec):
        """start with text file, end with data frame to generate the gtvt figures.

        Parameters
        ----------
        ra : float
            Right ascension of target
        dec : float
            Declination of target
        fixed : bool
            True if target is fixed, False if moving.

        Returns
        -------
        dataframe : pandas.DataFrame
            Pandas dataframe with updated metadata
        """
        self.fixed = True

        if ":" in ra:
            ra = round(self.convert_ddmmss_to_float(ra) * 15.0, 4)
            dec = round(self.convert_ddmmss_to_float(dec), 4)
        else:
            ra, dec = float(ra), float(dec)

        self.dataframe["ra"] = ra
        self.dataframe["dec"] = dec

        self.dataframe = self.build_dataframe()

        return self.dataframe

    def get_moving_target_positions(self, desg, smallbody):
        """Ephemeris from JPL/HORIZONS.
        desg : str
            Name of target
        Returns : target name from HORIZONS, RA, and Dec.
        """

        self.fixed = False

        if smallbody:
            id_type = "smallbody"
        else:
            id_type = None

        obj = Horizons(
            id=desg,
            location="500@-170",
            epochs={
                "start": self.start_date.to_value("iso", subfmt="date"),
                "stop": self.end_date.to_value("iso", subfmt="date"),
                "step": "1d",
            },
            id_type=id_type,
        )

        eph = obj.ephemerides(cache=False, quantities=(1))
        self.target_name = eph["targetname"][0]

        self.dataframe["ra"] = eph["RA"].data.data
        self.dataframe["dec"] = eph["DEC"].data.data

        self.dataframe = self.build_dataframe()

        return self.dataframe

    def in_FOR(self, dataframe):
        """Based on parameters calculated in tool, see if the target is in the Field of Regard

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Pandas dataframe

        Returns
        -------
        dataframe : pandas.DataFrame
            Pandas dataframe with updated metadata
        """
        dataframe["in_FOR"] = np.where(
            (dataframe["dist"] < MAX_SUN_ANGLE) & (dataframe["dist"] > MIN_SUN_ANGLE),
            True,
            False,
        )
        return dataframe

    def normal_pa(self, dataframe, tgt_ra, tgt_dec):
        """Calculate the normal position angle (V3)

        Parameters
        ----------
        tgt_ra : float
            Right ascension of target
        tgt_dec : float
            Declination of target

        Returns
        -------
        dataframe : pandas.DataFrame
            Pandas dataframe with updated metadata
        """

        if "sun_pa" not in dataframe.columns:
            dataframe = self.sun_position_coordinates(dataframe)
            sun_coord1, sun_coord2 = (
                dataframe["coord1"].values,
                dataframe["coord2"].values,
            )
            dataframe["sun_pa"] = self.calculate_sun_pa(
                dataframe, tgt_ra, tgt_dec, sun_coord1, sun_coord2
            )

        dataframe["V3PA"] = (
            dataframe["sun_pa"] + np.pi
        )  #  We want -V3 pointed towards sun.

        dataframe.loc[dataframe["V3PA"] < 0, "V3PA"] += 2 * np.pi
        dataframe.loc[dataframe["V3PA"] > 2 * np.pi, "V3PA"] -= 2 * np.pi

        dataframe["V3PA"] *= R2D

        return dataframe

    def build_dataframe(self):
        self.dataframe = self.sun_position_vectors(self.dataframe)

        self.dataframe = self.sun_position_coordinates(self.dataframe)

        self.dataframe["ra_radians"] = self.dataframe["ra"] * D2R
        self.dataframe["dec_radians"] = self.dataframe["dec"] * D2R

        self.dataframe = self.calculate_sun_pa(
            self.dataframe,
            self.dataframe["ra_radians"],
            self.dataframe["dec_radians"],
            self.dataframe["coord1"],
            self.dataframe["coord2"],
        )

        self.dataframe = self.normal_pa(
            self.dataframe, self.dataframe["ra_radians"], self.dataframe["dec_radians"]
        )

        self.dataframe = self.dist(self.dataframe)
        self.dataframe = self.in_FOR(self.dataframe)
        self.dataframe = self.get_allowed_max_boresight(self.dataframe)

        instrument_aperture_pairs = [
            ("NIRCAM", "NRCALL_FULL"),
            ("NIRSPEC", "NRS_FULL_MSA"),
            ("NIRISS", "NIS_CEN"),
            ("MIRI", "MIRIM_FULL"),
            ("FGS", "FGS1_FULL"),
            ("V3PA", None),
        ]

        for instrument, aperture in instrument_aperture_pairs:
            self.dataframe = self.calculate_min_max_pa_angles(
                self.dataframe, instrument, aperture
            )

        return self.dataframe

    def sun_position_vectors(self, dataframe):
        """Add the sun vector to ephemeris data frame

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Pandas dataframe containing ephemeride information

        Returns
        -------
        dataframe : pandas.DataFrame
            Pandas dataframe containing updated metadata.
        """

        # Create dataframe of only positional vectors.
        sun_vector_df = dataframe[["X", "Y", "Z"]] * -1

        # Calculate vector magnitude for each row.
        dataframe["Vsun_mag"] = sun_vector_df.apply(np.linalg.norm, axis=1)

        # Generate normal vectors.
        dataframe["Vsun_X"] = sun_vector_df["X"] / dataframe["Vsun_mag"]
        dataframe["Vsun_Y"] = sun_vector_df["Y"] / dataframe["Vsun_mag"]
        dataframe["Vsun_Z"] = sun_vector_df["Z"] / dataframe["Vsun_mag"]

        return dataframe

    def sun_position_coordinates(self, dataframe):
        """Get sun's positional coordinates

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Pandas dataframe containing ephemeride information

        Returns
        -------
        dataframe : pandas.DataFrame
            Pandas dataframe containing updated metadata.
        """

        dataframe["coord2"] = np.arcsin(dataframe["Vsun_Z"].values)
        dataframe["coord1"] = np.arctan2(
            dataframe["Vsun_Y"].values, dataframe["Vsun_X"].values
        )

        dataframe.loc[dataframe["coord1"] < 0, "coord1"] += 2 * np.pi

        return dataframe

    def write_ephemeris(self, data_frame, write_path):
        """Write ephemeris to file

        Parameters
        ----------
        write_path : str
            Full file path to write file out to.
        """
        if not write_path.endswith(".csv"):
            raise Exception(
                "Writing out files JWST GTVT uses a csv writer, please provide filename with '.csv' extension."
            )
        else:
            data_frame.to_csv(write_path, index=False)

    def update_ephemeris_data(self):
        """Function to write out new text based ephemeris to GTVT."""

        max_date = self.ephemeris_maximum_date()

        request_url = URL.format(LAUNCH_DATE, max_date)
        self.ephemeris_request = requests.get(request_url)

        ephemeris = np.array(self.eph_request.text.splitlines())
        ephemeris_filename = "ephemeris_{}_{}.txt".format(LAUNCH_DATE, max_date)

        with open(
            os.path.join(
                os.path.dirname(__file__), "data/{}".format(ephemeris_filename)
            ),
            "w",
        ) as fp:
            for entry in ephemeris:
                fp.write("%s\n" % entry)
