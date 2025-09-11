import pytest
from jwst_gtvt.jwst_tvt import Ephemeris
from astropy.time import Time

@pytest.fixture(scope="module")
def ephemeris_stable_date():
    arbitrary_start_date = Time("2024-09-10 10:30:00", format="iso", scale="utc")
    arbitrary_end_date = Time("2025-09-10 10:30:00", format="iso", scale="utc")
    return Ephemeris(start_date=arbitrary_start_date, end_date=arbitrary_end_date) # Need to pick a specific date so that tests are consistent


def test_random_test_experiment(ephemeris_stable_date):
    ra, dec = '16:52:58.9', '02:24:03'
    ephemeris_stable_date.get_fixed_target_positions(ra, dec)