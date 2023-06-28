1.0
===

- Package went through a major refactor, took an object oriented approach and rewrote most of the API and modules
- Tool now dynamically reads ephemerides from JPL HORIZONS using astroqury. Back up ephemerides are available in repo along with tools to update it.
- Plotting display and control is different but still displays results to previous releases
- Command line arguments now use `docopt`.
- Conda and pip requirements files are included with Package
- Project now uses `pyproject.toml` to setup and install package 

0.3.0
=====

- Updated the default JWST ephemeris and removed the old ephemeris.
- Updated the default start date to 2021-12-26 (JWST launch) and the default end date to 2024-10-02.
- Updated README examples and figures.

0.2.0
=====

- Adding support for the python >= 3.7
- Adding pysiaf library. pysiaf pulls the v3 reference angle for the GTVT/MTVT.
- Small print/output statements updates.
- Updating repo to follow INS community software standards.
- Switched the CI from travis to Github Actions

0.1.2 
=====

- verbose argument is available via command line. Users can select to suppress output from GTVT/MTVT.