1.1 (2025-06-03)
================
- Updated default start and end dates to cover the period from the current local date to 2 years in the future.
- Updated the help text to be more descriptive.
- Removed the .yml environment files from the repository.
- Updated the README to describe the new default dates and installation instructions.
- Added dates for each GTVT release to CHANGES.
- Removed support for Python 3.9

1.0.1 (2023-10-17)
==================
- Bug fix: Negative declinations resulted in the GTVT throwing an
  error because the argument parser accepts arguments that start with
  '--' and '-', so was unable to interpret the minus sign. Updates to
  the jwst_gtvt call were implemented to get around this issue by
  requiring '--ra=' and '--dec=' in front of the right ascension and
  declination values, respectively. The README was updated with this information.

1.0 (2023-06-29)
================
- Package went through a major refactor, list below describes new features
- Tool now dynamically reads ephemerides from JPL HORIZONS using astroquery. Back-up ephemerides are available in repo along with tools to update it.
- Plotting display and control is updated but still presents results similar to previous releases
- Command line arguments now use `docopt`.
- Conda and pip requirements files are included with Package
- Project now uses `pyproject.toml` to setup and install package 
- The --v3pa flag is now deprecated in favor of --instrument v3pa

0.3.0 (2022-11-10)
==================
- Updated the default JWST ephemeris and removed the old ephemeris.
- Updated the default start date to 2021-12-26 (JWST launch) and the default end date to 2024-10-02.
- Updated README examples and figures.

0.2.0 (2021-07-12)
==================
- Adding support for the python >= 3.7
- Adding pysiaf library. pysiaf pulls the v3 reference angle for the GTVT/MTVT.
- Small print/output statements updates.
- Updating repo to follow INS community software standards.
- Switched the CI from travis to Github Actions

0.1.2 (2020-04-23)
==================
- Verbose argument is available via command line. Users can select to suppress output from GTVT/MTVT.
