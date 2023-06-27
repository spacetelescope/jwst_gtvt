import numpy as np

D2R = np.pi/180.  # degrees to radians
EPSIOLON = 23.43929 * D2R  # obliquity of the ecliptic J2000
MIN_SUN_ANGLE = 84.8 * D2R  # minimum Sun angle, in radians
MAX_SUN_ANGLE = 135.0 * D2R  # maximum Sun angle, in radians
OBLIQUITY_OF_THE_ECLIPTIC = -23.439291 * D2R  # At J2000 equinox
PI2 = 2. * np.pi # 2 pi
R2D = 180. / np.pi  # radians to degrees 
SUN_ANGLE_PAD = 0.5 * D2R  # pad away from Sun angle limits when constructing safe attitude
UNIT_LIMIT = lambda x: min(max(-1.,x),1.)  # forces value to be in [-1,1]
URL = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND=-170&OBJ_DATA=%27NO%27&EPHEM_TYPE=VECTORS&START_TIME=%27{}%27&STOP_TIME=%27{}%27&CENTER=%27500@10%27&STEP_SIZE=%271%20DAYS%27&CSV_FORMAT=%27YES%27&VEC_TABLE=%272%27&REF_SYSTEM=%27ICRF%27&REF_PLANE=%27FRAME%27&VEC_CORR=%27LT%27&OUT_UNITS=%27KM-S%27&VEC_LABELS=%27YES%27&VEC_DELTA_T=%27NO%27'
