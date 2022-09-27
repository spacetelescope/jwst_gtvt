import numpy as np

D2R = np.pi/180.  #  degrees to radians
R2D = 180. / np.pi #  radians to degrees 
PI2 = 2. * np.pi   #  2 pi

UNIT_LIMIT = lambda x: min(max(-1.,x),1.) # forces value to be in [-1,1]

MIN_SUN_ANGLE = 84.8 * D2R  #  minimum Sun angle, in radians
MAX_SUN_ANGLE = 135.0 * D2R #  maximum Sun angle, in radians
SUN_ANGLE_PAD = 0.5 * D2R   #  pad away from Sun angle limits when constructing safe attitude

OBLIQUITY_OF_THE_ECLIPTIC = -23.439291 * D2R #  At J2000 equinox
EPSIOLON = 23.43929 * D2R #  obliquity of the ecliptic J2000