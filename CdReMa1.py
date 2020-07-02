# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:50:33 2020

@author: Wallace
"""

from math import pi, sqrt, exp, sin, cos
import numpy as np
# from geopy.geocoders import Nominatim

def CdReMa1(tbb,xbb):
    D = 0.1143; # diametro maximo do projetil [m]
    Db = 0.1107; # diametro da base do projetil [m]
    p0 = 1.01325e5; # Pa
    R = 287.05; # J/[Kg.K]
    ka = 1.4; # constante do ar [relacao entre Cp e Cv]
    T0 = 288.15; # K
    Ts = 110.4; # constante de Sutherland [K]
    mu0= 1.79e-5; # constante kg/[m.s.K**1/2]
    g = 9.80665; # aceleracao da gravidade [m/s^2]
    L = 6.5e-3; # constante de dilatacao adiabatica [K/m]
    omega = 7.292e-05 # [rad/s]     
    
    # extraindo dados do array:    
    H = xbb.item(1)
    u = xbb.item(3)
    v = xbb.item(4)
    w = xbb.item(5)
    m = xbb.item(6)
    vel = sqrt(u**2+v**2+w**2)
        
    if (H <= 11000):
      T = T0 - L*H; 
      p = p0*(1 + L*H/T0)**(-g/L/R);
    else:
      T = 216.65;
      p = (2.2632e4)*exp(-g/R/T*(H-11000));
      
    rho = p/R/T
    
    mu = mu0*(T0+Ts)/(T+Ts)*((T/T0))**(1.5) # viscosidade do ar em funcao da temperatura
    vsom = sqrt(ka*R*T) # velocidade do som
    Ma = vel/vsom # numero de Mach
    Rey = rho*D*vel/mu # numero de reynolds associado
    A = 0.25*pi*D**2; # area de referencia do base bleed
    Ab = 0.25*pi*Db**2; # area da base do base bleed
    
    if Ma < 1.0:
        Cdrag = 0.2
    else: 
        Cdrag = 0.3 - 0.1*(Ma-1) 
    
    valBB = np.zeros((1,7))
    valReMa = np.zeros((1,3))
    valPRT = np.zeros((1,3))
    
    mb = 0.005*rho*vel*Ab
    
    # CORIOLIS EFFECT: 2*omega*[Cx,Cy,Cz] 
    '''geoloc = Nominatim() # trazendo o objeto para facilitar o estudo
    loc = geoloc.geocode("Rio de Janeiro, Rio de Janeiro, Brazil") # coordenadas RJ'''
    lat = (pi/180)*(-23); # latitude do RJ em radianos
    az = 0; # azimute do lancamento
    
    Cx = -v*cos(lat)*sin(az) - w*sin(lat) 
    Cy = u*cos(lat)*sin(az) + w*cos(lat)*cos(az)
    Cz = u*sin(lat) - v*cos(lat)*cos(az)    
    
    valBB[0,0] = u
    valBB[0,1] = v
    valBB[0,2] = w
    valBB[0,3] = -(0.5)*Cdrag*rho*vel*u*A/m + u*mb/(vel*m) + 2*omega*Cx
    valBB[0,4] = -(0.5)*Cdrag*rho*vel*v*A/m + v*mb/(vel*m) - g + 2*omega*Cy
    valBB[0,5] = -(0.5)*Cdrag*rho*vel*w*A/m + 2*omega*Cz
    valBB[0,6] = -mb
    
    valReMa[0,0] = Cdrag
    valReMa[0,1] = Ma
    valReMa[0,2] = Rey
    
    valPRT[0,0] = p
    valPRT[0,1] = rho
    valPRT[0,2] = T
    
    return valBB, valReMa, valPRT