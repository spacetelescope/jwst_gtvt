import pytest
from jwst_gtvt.jwst_tvt import Ephemeris
from astropy.time import Time
import pandas as pd 
import numpy as np


### FIXTURES ###


@pytest.fixture(scope="module")
def ephemeris():

    # Keeping constant dates to ensure tests are idempotent
    arbitrary_start_date = Time("2024-09-10 10:30:00", format="iso", scale="utc")
    arbitrary_end_date = Time("2025-09-10 10:30:00", format="iso", scale="utc")

    return Ephemeris(start_date=arbitrary_start_date, end_date=arbitrary_end_date)


### TESTS ###


def test_ephemeris_in_FOR(ephemeris):

    # Copied not imported because they shouldn't change, so we should catch changes
    MIN_SUN_ANGLE = 84.8 * np.pi / 180.0
    MAX_SUN_ANGLE = 135.0 * np.pi / 180.0

    original_df = pd.DataFrame({
        "dist" : [MIN_SUN_ANGLE, MAX_SUN_ANGLE, MIN_SUN_ANGLE + 0.1, MAX_SUN_ANGLE - 0.1, MIN_SUN_ANGLE - 0.1, MAX_SUN_ANGLE + 0.1, 1, 2, 3]
    })

    result_df = ephemeris.in_FOR(original_df)

    assert "in_FOR" in result_df

    expected = pd.Series([False, False, True, True, False, False, False, True, False])
    assert result_df["in_FOR"].equals(expected)


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
def test_ephemeris_convert_ddmmss_to_float_success(ephemeris, input, expected):

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
def test_ephemeris_convert_ddmmss_to_float_bad_input(ephemeris, input):

    ephemeris.convert_ddmmss_to_float(input)