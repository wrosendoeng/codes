# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:50:32 2020

@author: Wallace
"""
from math import sin, cos, pi, sqrt, exp
import numpy as np
import pandas as pd

def CdReMa0(tsb,xsb):
    msb = 20.318;
    D = 0.1143; # diametro do projetil
    p0 = 1.01325e5; # Pa
    R = 287.05; # J/[Kg.K]
    ka = 1.4; # constante do ar [relacao entre Cp e Cv]
    T0 = 288.15; # K
    Ts = 110.4; # constante de Sutherland [K]
    mu0= 1.79e-5; # constante kg/[m.s.K**1/2]   
    g = 9.80665; # aceleração da gravidade [m/s^2]   
    L = 6.5e-3; # constante de dilatacao adiabatica [K/m]     
    omega = 7.292e-05; # velocidade angular da Terra [rad/s]
    
    H = xsb.item(1)
    u = xsb.item(3)
    v = xsb.item(4)
    w = xsb.item(5)
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
    Ma = (vel/vsom) # numero de Mach
    Rey = rho*D*vel/mu # numero de reynolds associado
    A = 0.25*pi*D**2
        
    """if Ma < 1.2:
        Cdrag0 = 0.35*(Ma)**2 - 0.27*Ma + 0.14
    elif Ma >= 1.2 and Ma <= 1.3:
        Cdrag0 = 0.32
    else:
        Cdrag0 = 0.3876 - 0.052*Ma"""
        
    df_mach_cd = pd.read_csv('navy_aero.dat', sep='\s+', header=None, skiprows=1)
    cd_function = np.polyfit(df_mach_cd[0],df_mach_cd[1], 4)
    Cdrag0 = np.polyval(cd_function,Ma)
        
    valFF = np.zeros([1,6])
    valReMa = np.zeros([1,3])
    valPRT = np.zeros([1,3])
    
    # CORIOLIS EFFECT: 2*omega*[Cx,Cy,Cz]
    '''geoloc = Nominatim() # trazendo o objeto para facilitar o estudo
    loc = geoloc.geocode("Rio de Janeiro, Rio de Janeiro, Brazil") # coordenadas RJ'''
    lat = (pi/180)*(-23); # latitude do RJ em radianos
    az = 0; # azimute do lancamento
    
    Cx = -v*cos(lat)*sin(az) - w*sin(lat) 
    Cy = u*cos(lat)*sin(az) + w*cos(lat)*cos(az)
    Cz = u*sin(lat) - v*cos(lat)*cos(az) 
    
    valFF[0,0] = u
    valFF[0,1] = v
    valFF[0,2] = w
    valFF[0,3] = -(0.5)*Cdrag0*rho*vel*u*A/msb + 2*omega*Cx   
    valFF[0,4] = -(0.5)*Cdrag0*rho*vel*v*A/msb - g + 2*omega*Cy
    valFF[0,5] = -(0.5)*Cdrag0*rho*vel*w*A/msb + 2*omega*Cz
    
    valReMa[0,0] = Cdrag0
    valReMa[0,1] = Ma
    valReMa[0,2] = Rey
    
    valPRT[0,0] = p
    valPRT[0,1] = rho
    valPRT[0,2] = T
    
    return [valFF, valReMa, valPRT]