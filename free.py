# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:50:32 2020

@author: Wallace
"""

import numpy as np
from math import sin, cos, pi
# from geopy.geocoders import Nominatim

def free(t,x): 
    
    g = 9.80665;
    omega = 7.2929e-05; 
    u = x.item(3)
    v = x.item(4)
    w = x.item(5)
               
    valFF = np.zeros([1,6])
    
    # CORIOLIS EFFECT: 2*omega*[Cx,Cy,Cz]
    '''geoloc = Nominatim() # trazendo o objeto para facilitar o estudo
    loc = geoloc.geocode("Rio de Janeiro, Rio de Janeiro, Brazil") # coordenadas RJ'''
    lat = (pi/180)*(-23) # latitude do RJ em radianos
    az = 0 # azimute do lancamento
    
    Cx = -v*cos(lat)*sin(az) - w*sin(lat) 
    Cy = u*cos(lat)*sin(az) + w*cos(lat)*cos(az)
    Cz = u*sin(lat) - v*cos(lat)*cos(az)
    
    valFF[0,0] = u
    valFF[0,1] = v
    valFF[0,2] = w
    valFF[0,3] = 2*omega*Cx   
    valFF[0,4] = - g + 2*omega*Cy
    valFF[0,5] = 2*omega*Cz
    
    return valFF