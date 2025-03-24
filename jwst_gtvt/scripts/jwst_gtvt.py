usage = """

Driver for JWST GTVT fixed target tool.

Usage:
  jwst_gtvt --ra=<ra> --dec=<dec> [--start_date=<obs_start>] [--end_date=<obs_end>] [--instrument=<inst>] [--target_name=<name>] [--write_ephemeris=<write_path>] [--write_plot=<plot_path>] [--silent]

Arguments:
  --ra=<ra>     Right ascension of target in degrees
  --dec=<dec>   Declination of target in degrees

Options:
  [--start_date]         Start date for plot (YYYY-MM-DD)
  [--end_date]           End date for plot (YYYY-MM-DD)
  [--instrument]         Instrument (fgs, miri, nircam, niriss, nirspec, default=all)
  [--target_name]        User provided name for target (name for output, double-quoted if there are spaces)
  [--write_ephemeris]    File name to write ephemeris to
  [--write_plot]         File name to write plot out to
  [--silent]             Boolean to print results to screen [default: False]
  --help                 Show this screen.
  --version              Show version.
"""

from astropy.time import Time

from jwst_gtvt.constants import D2R
from jwst_gtvt.display_results import display_results
from docopt import docopt
from jwst_gtvt.jwst_tvt import Ephemeris
from jwst_gtvt.plotting import plot_visibility


def main(args):
    if args['--start_date'] and args['--end_date']:
        start = Time(args['--start_date'])
        end = Time(args['--end_date'])
        if start > end:
            raise ValueError("Date mismatch, end date is before start date.")
        else:
            eph = Ephemeris(start_date=start, end_date=end)
    else:
        eph = Ephemeris()

    eph.get_fixed_target_positions(args['--ra'], args['--dec'])

    if not eph.dataframe['in_FOR'].any():
        in_FOR_msg = ("No position angles in field of regard! "
                      "Check constraints for your target and if it is observable with JWST. \n"
                      "Vist: https://jwst-docs.stsci.edu/jwst-observatory-characteristics/jwst-observatory-coordinate-system-and-field-of-regard for more information")
        raise IndexError(in_FOR_msg)

    if args['--write_ephemeris']:
        args['--silent'] = True
        eph.write_ephemeris(eph.dataframe, args['--write_ephemeris'])

    if not args['--silent']:
            display_results(eph)

    plot_visibility(eph, args['--instrument'], name=args['--target_name'], write_plot=args['--write_plot'])




def driver():
    args = docopt(usage, version='0.6.2')
    main(args)
