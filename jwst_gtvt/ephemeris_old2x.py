#! /usr/bin/env python
#Module ephemeris.py
from __future__ import print_function

from astroquery.jplhorizons import Horizons
from math import *
import numpy as np
import urllib
import sys
import os
import urllib
from .rotationsx import *
from . import astro_funcx as astro_func

D2R = pi/180.  # degrees to radians
R2D = 180. / pi # radians to degrees 
PI2 = 2. * pi   # 2 pi
unit_limit = lambda x: min(max(-1.,x),1.) # forces value to be in [-1,1]
MIN_SUN_ANGLE = 84.8 * D2R  # minimum Sun angle, in radians
MAX_SUN_ANGLE = 135.0 * D2R # maximum Sun angle, in radians
SUN_ANGLE_PAD = 0.5 * D2R   # pad away from Sun angle limits when constructing safe attitude

obliquity_of_the_ecliptic = -23.439291  # At J2000 equinox
obliquity_of_the_ecliptic *=  D2R
Qecl2eci = QX(obliquity_of_the_ecliptic)


class Ephemeris:
    def __init__(self,start_date,end_date):
        """Eph constructor using astroquery.jplhorizons. Need to add ecliptic conversion in the future
        
        Parameters
        ----------
        start_date : str
            Observation window min date in format YYYY-MM-DD
        end_date : str
            Observation window max date in format YYYY-MM-DD
        """

        obj = Horizons(id='jwst', id_type='id',  location=None, 
                       epochs={'start':start_date, 'stop':end_date, 'step':'1d'})
        
        vectors = obj.vectors()
        au_to_km = 1.496e8

        self.datelist = vectors['datetime_jd'] - 2400000.5
        self.xlist = vectors['x'] * au_to_km
        self.ylist = vectors['y'] * au_to_km
        self.zlist = vectors['z'] * au_to_km

        self.amin = min(self.datelist)
        self.amax = max(self.datelist)


    def pos(self,adate):
        """Get positions
        """
        cal_days = adate - self.datelist[0]
        indx = int(cal_days)
        frac = cal_days - indx

        x = (self.xlist[indx+1] - self.xlist[indx])*frac + self.xlist[indx]  
        y = (self.ylist[indx+1] - self.ylist[indx])*frac + self.ylist[indx]  
        z = (self.zlist[indx+1] - self.zlist[indx])*frac + self.zlist[indx]  
        
        return np.array([x,y,z])


    def Vsun_pos(self,adate):
        """Calculate the sun's velocity unit vector? (v/|v|).
        """

        # -1? Does this have something to do with the direction of the vector?
        Vsun = np.negative(self.pos(adate))/np.linalg.norm(self.pos(adate)) 

        return Vsun
    
    
    def sun_pos(self,adate):
        """Calculate the sun's position?
        """
        
        Vsun = self.Vsun_pos(adate)
        
        coord2 = np.arcsin(Vsun[2])
        coord1 = np.arctan2(Vsun[1], Vsun[0])
        
        if coord1 < 0.: 
            coord1 += PI2
        
        return (coord1,coord2)


    def normal_pa(self,adate,tgt_c1,tgt_c2):
        (sun_c1, sun_c2) = self.sun_pos(adate)
        sun_pa = astro_func.pa(tgt_c1,tgt_c2,sun_c1,sun_c2)
        
        V3_pa = sun_pa + pi  # We want -V3 pointed towards sun.
        
        if V3_pa < 0.: 
            V3_pa += PI2
        
        if V3_pa >= PI2: 
            V3_pa -= PI2
        
        return V3_pa


    def is_valid(self,date,coord_1,coord_2,V3pa):
        """Indicates whether an attitude is valid at a given date."""
        
        # First check that the date is within the time interval of the ephemeris.
        if ((date < self.amin) or (date > self.amax)):
            return False
            
        (sun_1,sun_2) = self.sun_pos(date)
        d = astro_func.dist(coord_1,coord_2,sun_1,sun_2)
        vehicle_pitch = pi/2 - d   #see JI memo from May 2006
        # sun pitch is always equal or greater than sun angle (V1 to sun)
        if (d<MIN_SUN_ANGLE or d>MAX_SUN_ANGLE):
            return False
        
        pa = astro_func.pa(coord_1, coord_2, sun_1, sun_2) + pi
        roll = acos(cos(V3pa - pa))
        sun_roll = asin(sin(roll) * cos(vehicle_pitch))
        
        if (abs(sun_roll)<=5.2*D2R):
            sun_pitch = atan2(tan(vehicle_pitch), cos(roll))
            if (sun_pitch<=5.0*D2R and sun_pitch>=-44.8*D2R):
                return True
        
        return False


    def in_FOR(self,adate,coord_1,coord_2):
        """ Check if the target is in the field of regard.
        """
        (sun_1,sun_2) = self.sun_pos(adate)
        d = astro_func.dist(coord_1,coord_2,sun_1,sun_2)
        # print d*R2D
        # 90 - sun pitch is always equal or greater than sun angle (V1 to sun)
        if (d<MIN_SUN_ANGLE or d>MAX_SUN_ANGLE):
            return False
        return True


    def bisect_by_FOR(self,in_date,out_date,coord_1,coord_2):
        """Why would one bisect by field of regard?
        """
        # FOR = Field Of Regard
        # in and out of FOR, assumes only one "root" in interval
        delta_days = 200.
        mid_date = (in_date+out_date)/2.
        while delta_days > 0.000001:
            (sun_1,sun_2) = self.sun_pos(mid_date)
            d = astro_func.dist(coord_1,coord_2,sun_1,sun_2)
            if (d>MAX_SUN_ANGLE or d<MIN_SUN_ANGLE):
                out_date = mid_date
            else:
                in_date = mid_date
            mid_date = (in_date+out_date)/2.
            delta_days = abs(in_date-out_date)/2.

        if in_date>out_date:
            # ensure returned date always in field of regard
            mid_date = mid_date + 0.000001
        else:
            mid_date = mid_date - 0.000001
        return mid_date


    def bisect_by_attitude(self,in_date,out_date,coord_1,coord_2,pa):
        """ Why would one bisect by attitude?
        """
        # Attitude is the orientation of the craft with respect 
        # to a set of reference axes
        
        # in and out of FOR, assumes only one "root" in interval
        icount = 0
        delta_days = 200.
        mid_date = (in_date+out_date)/2.
        
        while delta_days > 0.000001:
            if self.is_valid(mid_date,coord_1,coord_2,pa):
                in_date = mid_date
            else:
                out_date = mid_date
            
            mid_date = (in_date+out_date)/2.
            delta_days = abs(in_date-out_date)/2.
            
            icount = icount + 1
        
        return mid_date