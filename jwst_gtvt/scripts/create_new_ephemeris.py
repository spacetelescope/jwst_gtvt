from jwst_gtvt.jwst_tvt import Ephemeris


def driver():
    eph = Ephemeris()
    eph.update_ephemeris_data()
