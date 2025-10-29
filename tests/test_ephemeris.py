import pytest
from jwst_gtvt.jwst_tvt import Ephemeris
from astropy.time import Time
import pandas as pd 
import numpy as np
from datetime import datetime


"""
A test suite that focuses on having correct output for class methods within the Ephemeris Class

"""

### FIXTURES ###


@pytest.fixture(scope="module")
def ephemeris():

    # Keeping constant dates to ensure tests are idempotent
    arbitrary_start_date = Time("2024-09-10 10:30:00", format="iso", scale="utc")
    arbitrary_end_date = Time("2025-09-10 10:30:00", format="iso", scale="utc")

    return Ephemeris(start_date=arbitrary_start_date, end_date=arbitrary_end_date)


### TESTS ###

# TODO: Refactor to use parameters to create more isolated tests, do in for cols
def test_in_FOR(ephemeris):

    # Copied not imported because they shouldn't change, so we should catch changes
    MIN_SUN_ANGLE = 84.8 * np.pi / 180.0
    MAX_SUN_ANGLE = 135.0 * np.pi / 180.0

    original_df = pd.DataFrame({
        "dist" : [MIN_SUN_ANGLE, MAX_SUN_ANGLE, MIN_SUN_ANGLE + 0.1, MAX_SUN_ANGLE - 0.1, MIN_SUN_ANGLE - 0.1, MAX_SUN_ANGLE + 0.1, 1, 2, 3]
    })

    result_df = ephemeris.in_FOR(original_df)

    expected = pd.Series([False, False, True, True, False, False, False, True, False])
    
    assert "in_FOR" in result_df and result_df["in_FOR"].equals(expected)


@pytest.mark.parametrize(
    "input, expected", 
    [
        ("0:0:0", 0),
        ("-0:0:0", 0),
        ("12:30:15", 12.5041666667),
        ("23:59:59", 23.9997222222),
        ("-1:15:37", -1.2602777778),
        ("-89:59:59", -89.9997222222),
    ]
)
def test_convert_ddmmss_to_float_success(ephemeris, input, expected):

    assert ephemeris.convert_ddmmss_to_float(input) == pytest.approx(expected)


@pytest.mark.xfail(reason="Faulty Inputs")
@pytest.mark.parametrize(
    "input",
    [
        "",
        "3:6:",
        "77",
        "a:b:c",
        "5:5:5:5",
    ]
)
def test_convert_ddmmss_to_float_bad_input(ephemeris, input):

    ephemeris.convert_ddmmss_to_float(input)

@pytest.mark.parametrize(
    "input_ra, input_dec, expected",
    [
        (0, 0, 0),
        (180, 0, 0),
        (0, 90, 66.560719),
        (270, 0, 23.439281),
        (200, 10, 17.055630),
    ]
)
def test_calculate_ecliptic_latitude(ephemeris, input_ra, input_dec, expected):

    result = ephemeris.calculate_ecliptic_latitude(input_ra, input_dec)
    assert result == pytest.approx(expected, abs=1e-5)

@pytest.mark.parametrize(
        "input_sun_ra, input_sun_dec, input_ra, input_dec, expected",
        [
            (0, 0, 0, 0, -np.pi / 2),
            (0, 0, np.pi / 4, 0, -0.6682893715286274),
            (np.pi / 2, -np.pi / 3, np.pi / 6, np.pi / 4, 0.09891704637493678),
        ]
)
def test_allowed_max_vehicle_roll(ephemeris, input_sun_ra, input_sun_dec, input_ra, input_dec, expected):
    
    assert ephemeris.allowed_max_vehicle_roll(input_sun_ra, input_sun_dec, input_ra, input_dec) == pytest.approx(expected)


@pytest.mark.parametrize(
    "input_sun_p, expected",
    [
        (0, 5.1 * np.pi / 180),
        (2.5 * np.pi / 180, 5.1 * np.pi / 180),
        (3.85 * np.pi / 180, 4.25 * np.pi / 180),
        (5.2 * np.pi / 180, 3.4 * np.pi / 180),
    ]
)
def test_allowed_max_sun_roll(ephemeris, input_sun_p, expected):

    assert ephemeris.allowed_max_sun_roll(input_sun_p) == pytest.approx(expected)


@pytest.mark.parametrize(
    "obj1_c1, obj1_c2, obj2_c1, obj2_c2, expected",
    [
        (0, 0, 0, 0, 0),
        (0, 0, np.pi/2, 0, np.pi/2),
        (0, np.pi/4, 0, -np.pi/4, np.pi/2),
    ]
)
def test_angular_sep_equal(ephemeris, obj1_c1, obj1_c2, obj2_c1, obj2_c2, expected):

    assert ephemeris.angular_sep(obj1_c1, obj1_c2, obj2_c1, obj2_c2) == pytest.approx(expected)


@pytest.mark.parametrize(
    "obj1_c1, obj1_c2, obj2_c1, obj2_c2",
    [
        (0, np.pi/4, 2*np.pi, np.pi/4),
        (0, 0, np.pi/2, 0),
        (0, np.pi/4, 0, -np.pi/4),
    ]
)
def test_angular_sep_communative(ephemeris, obj1_c1, obj1_c2, obj2_c1, obj2_c2):

    assert ephemeris.angular_sep(obj1_c1, obj1_c2, obj2_c1, obj2_c2) == ephemeris.angular_sep(obj2_c1, obj2_c2, obj1_c1, obj1_c2)


@pytest.mark.parametrize(
    "tgt_coord1, tgt_coord2, sun_coord1, sun_coord2, expected",
    [
        (0, 0, 0, 0, 0.0),
        (np.deg2rad(180), np.deg2rad(0), np.deg2rad(179), np.deg2rad(0), np.deg2rad(270)),
        (np.deg2rad(-90), np.deg2rad(90), np.deg2rad(-60), np.deg2rad(60), np.deg2rad(150)),
    ]
)
def test_calculate_sun_pa(ephemeris, tgt_coord1, tgt_coord2, sun_coord1, sun_coord2, expected):
    original_df = pd.DataFrame({
        "placeholder" : [0]
    })

    result_df = ephemeris.calculate_sun_pa(original_df, tgt_coord1, tgt_coord2, sun_coord1, sun_coord2)

    result_series = result_df["sun_pa"]
    expected_series = pd.Series([expected])

    assert result_series.equals(expected_series)
    

def test_maximum_date(ephemeris):

    result = datetime.strptime(ephemeris.ephemeris_maximum_date(), "%Y-%m-%d")
    expected = datetime(2030, 3, 16)

    # Assuming maximum date is variable, but should always be at least this fallback date
    assert result > expected

@pytest.mark.parametrize(
    "instrument, aperture, angle_name, expected",
    [
        ("NIRCAM", "NRCALL_FULL", "V3IdlYAngle", -0.07457694),
        ("NIRSPEC", "NRS_FULL_MSA", "V3IdlYAngle", 138.5745697),
        ("NIRISS", "NIS_CEN", "V3IdlYAngle", 0.56126717),
        ("MIRI", "MIRIM_FULL", "V3IdlYAngle", 4.83544897),
        ("FGS", "FGS1_FULL", "V3IdlYAngle", -1.24120427),
    ]
)
def test_get_angle(ephemeris, instrument, aperture, angle_name, expected):

    result = ephemeris.get_angle(instrument, aperture, angle_name)
    assert result == expected

@pytest.mark.parametrize(
    "sun_ra, sun_dec, ra, dec, expected",
    [
        (0, 0, 0, 0, -90.0),
        (0, 0, np.pi / 4, 0, -38.29016),
        (np.pi / 2, -np.pi / 3, np.pi / 6, np.pi / 4, 5.667529),
    ]
)
def test_get_allowed_max_boresight(ephemeris, sun_ra, sun_dec, ra, dec, expected):

    input_df = pd.DataFrame({
        "coord1" : [sun_ra],
        "coord2" : [sun_dec],
        "ra_radians" : [ra],
        "dec_radians" : [dec]
    })

    expected_df_output = pd.DataFrame({
        "max_boresight" : [expected]
    })

    df_output = ephemeris.get_allowed_max_boresight(input_df)

    assert ("max_boresight" in df_output.columns) and np.allclose(df_output["max_boresight"], expected_df_output["max_boresight"], atol=1e-6)

@pytest.mark.parametrize(
    "instrument, aperture, v3pa, max_boresight, expected_min, expected_nom, expected_max",
    [
        ("NIRCAM",  "NRCALL_FULL",  82.55,      3 * np.pi / 180,    82.423063,  82.475423,  82.527783),
        ("NIRCAM",  "NRCALL_FULL",  0.0012,     4.22 * np.pi / 180, 359.85297,  359.926623, 0.000276),
        ("MIRI",    "MIRIM_FULL",   18.192831,  0.98 * np.pi / 180, 23.011176,  23.02828,   23.045384),
        ("MIRI",    "MIRIM_FULL",   359.935,    2.5 * np.pi / 180,  4.726816,   4.770449,   4.814082),
        ("V3PA",    None,           16,         5.2 * np.pi / 180,  15.909243,  16,         16.090757), 
    ]
)
def test_calculate_min_max_pa_angles(ephemeris, instrument, aperture, v3pa, max_boresight, expected_min, expected_nom, expected_max):

    input_df = pd.DataFrame({
        "V3PA" : [v3pa],
        "max_boresight" : [max_boresight]
    })

    min_col_name = instrument + "_min_pa_angle"
    nom_col_name = instrument + "_nominal_angle"
    max_col_name = instrument + "_max_pa_angle"

    expected_df_output = pd.DataFrame({
        min_col_name : [expected_min],
        nom_col_name : [expected_nom],
        max_col_name : [expected_max],
    })

    df_output = ephemeris.calculate_min_max_pa_angles(input_df, instrument, aperture, angle_name="V3IdlYAngle")

    assert (
        np.allclose(df_output[min_col_name], expected_df_output[min_col_name], atol=1e-6)
        and np.allclose(df_output[nom_col_name], expected_df_output[nom_col_name], atol=1e-6)
        and np.allclose(df_output[max_col_name], expected_df_output[max_col_name], atol=1e-6)
    )

@pytest.mark.parametrize(
    "sun_x, sun_y, sun_z, expected_coord1, expected_coord2",
    [
        (1, 0, 0, 0, 0),
        (0, 1, 0, np.pi / 2, 0),
        (0, 0, -1, 0,  - np.pi / 2),
        (1 / np.sqrt(2), 1 / np.sqrt(2), 0, np.pi / 4, 0),
        (1 / np.sqrt(2), - 1 / np.sqrt(2), 0, 7 * np.pi / 4, 0),
        (0.5, 0.5, 0.5, np.pi / 4, np.pi / 6),
        (-0.35, 0.22, 0.3, 2.580429, 0.304693),
        (-0.7214, 0.5627, -0.4068, 2.479159, -0.418948)
    ]
)
def test_sun_position_coordinates(ephemeris, sun_x, sun_y, sun_z, expected_coord1, expected_coord2):

    input_df = pd.DataFrame({
        "Vsun_X" : [sun_x],
        "Vsun_Y" : [sun_y],
        "Vsun_Z" : [sun_z],
    })

    expected_df_output = pd.DataFrame({
        "coord1" : [expected_coord1],
        "coord2" : [expected_coord2],
    })

    df_output = ephemeris.sun_position_coordinates(input_df)

    assert (
        np.allclose(df_output["coord1"], expected_df_output["coord1"], atol=1e-6)
        and np.allclose(df_output["coord2"], expected_df_output["coord2"], atol=1e-6)
    )



@pytest.mark.parametrize(
    "tgt_ra, tgt_dec, sun_pa, expected_v3pa",
    [
        (0, 0, 0, 180),
        (3.008, 0.730, 4 * np.pi / 3, 60),
        (5.360, -0.411, 1.73211, 279.242593),
    ]
)
def test_normal_pa_with_sunpa(ephemeris, tgt_ra, tgt_dec, sun_pa, expected_v3pa):

    input_df = pd.DataFrame({
        "sun_pa" : [sun_pa]
    })

    expected_df_output = pd.DataFrame({
        "V3PA" : [expected_v3pa]
    })

    df_output = ephemeris.normal_pa(input_df, tgt_ra, tgt_dec)

    assert (
        np.allclose(df_output["V3PA"], expected_df_output["V3PA"], atol=1e-6)
    )


@pytest.mark.parametrize(
    "sun_x, sun_y, sun_z, tgt_ra, tgt_dec, expected_v3pa",
    [
        (1, 0, 0, 0, 0, 180)
    ]
)
def test_normal_pa_no_sunpa(ephemeris, sun_x, sun_y, sun_z, tgt_ra, tgt_dec, expected_v3pa):

    input_df = pd.DataFrame({
        "Vsun_X" : [sun_x],
        "Vsun_Y" : [sun_y],
        "Vsun_Z" : [sun_z],
    })

    expected_df_output = pd.DataFrame({
        "V3PA" : [expected_v3pa],
    })

    df_output = ephemeris.normal_pa(input_df, tgt_ra, tgt_dec)

    assert (
        np.allclose(df_output["V3PA"], expected_df_output["V3PA"], atol=1e-6)
    )


@pytest.mark.parametrize(
    "x, y, z, expectedVX, expectedVY, expectedVZ, expectedMag",
    [
        (3, 0, 0, -1, 0, 0, 3),
        (10, 2, -11, -10/15, -2/15, 11/15, 15),
        (0.347106, -0.812547, 0.563441, -0.331227, 0.775375, -0.537665, 1.047941),
        (-843.162774, 1659.234981, -247.398215, 0.449076, -0.883724, 0.131767, 1877.548951),
    ]
)
def test_sun_position_vectors(ephemeris, x, y, z, expectedVX, expectedVY, expectedVZ, expectedMag):
    input_df = pd.DataFrame({
        "X" : [x],
        "Y" : [y],
        "Z" : [z],
    })

    expected_df_output = pd.DataFrame({
        "Vsun_mag" : [expectedMag],
        "Vsun_X" : [expectedVX],
        "Vsun_Y" : [expectedVY],
        "Vsun_Z" : [expectedVZ],
    })

    df_output = ephemeris.sun_position_vectors(input_df)

    assert (
        np.allclose(df_output["Vsun_mag"], expected_df_output["Vsun_mag"], atol=1e-6)
        and np.allclose(df_output["Vsun_X"], expected_df_output["Vsun_X"], atol=1e-6)
        and np.allclose(df_output["Vsun_Y"], expected_df_output["Vsun_Y"], atol=1e-6)
        and np.allclose(df_output["Vsun_Z"], expected_df_output["Vsun_Z"], atol=1e-6)
    )


# Same values as tests for angular_sep, just ensuring they are passed correctly
@pytest.mark.parametrize(
    "obj1_c1, obj1_c2, obj2_c1, obj2_c2, expectedDist",
    [
        (0, 0, 0, 0, 0),
        (0, 0, np.pi/2, 0, np.pi/2),
        (0, np.pi/4, 0, -np.pi/4, np.pi/2),
    ]
)
def test_dist(ephemeris, obj1_c1, obj1_c2, obj2_c1, obj2_c2, expectedDist):

    input_df = pd.DataFrame({
        "ra_radians" : [obj1_c1],
        "dec_radians" : [obj1_c2],
        "coord1" : [obj2_c1],
        "coord2" : [obj2_c2]
    })

    expected_df_output = pd.DataFrame({
        "dist" : [expectedDist]
    })
    
    df_output = ephemeris.dist(input_df)

    assert (
        np.allclose(df_output["dist"], expected_df_output["dist"], atol=1e-6)
    )

