# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 20:59:16 2020

@author: Wallace
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 20:14:05 2020

@author: wrosendoeng
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
from math import pi, sin, cos
from CdReMa0 import CdReMa0
from CdReMa1 import CdReMa1
from free import free

elevation_angle = 45.;
azimuth_angle = 0.;
initial_speed = 878.;
elev_rad = pi*elevation_angle/180;
azim_rad = pi*azimuth_angle/180;
num_rows = int(100000); # estimativa
dt = 0.005; # intervalo de tempo [s]
    
def trajectory():

    #Deslocamento livre, sem considerar o arrasto:
    T = [0];
    t = 0;
    x = np.array([0,0,0,initial_speed*cos(elev_rad)*cos(azim_rad),initial_speed*sin(elev_rad)*cos(azim_rad),initial_speed*sin(azim_rad)]);
    Y = x.reshape((1,6));   

    #SEM base bleed:
    T0 = [0];
    tsb = 0; 
    xsb = np.array([0,0,0,initial_speed*cos(elev_rad)*cos(azim_rad),initial_speed*sin(elev_rad)*cos(azim_rad),initial_speed*sin(azim_rad)]);
    Ysb = xsb.reshape((1,6));
    Psb = CdReMa0(tsb,xsb)[2];
    Rsb = CdReMa0(tsb,xsb)[1];
    PRTsb = Psb.reshape((1,3));
    ReMasb = Rsb.reshape((1,3));
    
    #COM base bleed: 
    T1 = [0];
    tbb = 0;
    xbb = np.array([0,0,0,initial_speed*cos(elev_rad)*cos(azim_rad),initial_speed*sin(elev_rad)*cos(azim_rad),initial_speed*sin(azim_rad),20.859]);
    Ybb = xbb.reshape((1,7));
    Pbb = CdReMa1(tbb,xbb)[2];
    Rbb = CdReMa1(tbb,xbb)[1];
    PRTbb = Pbb.reshape((1,3));
    ReMabb = Rbb.reshape((1,3));
    
    for row in range(num_rows):
    
        #Deslocamento livre, sem considerar o arrasto:
        j1 = dt*free(t,x);
        j2 = dt*free(t+0.5*dt,x+0.5*j1);
        j3 = dt*free(t+0.5*dt,x+0.5*j2);
        j4 = dt*free(t+1.0*dt,x+1.0*j3); 
           
        x = x+(1/6)*(j1+2*(j2+j3)+j4);
        t = t+dt;
        T.append(t)
        Y = np.append(Y,x,axis=0) 
        
        if (Y[row+1,1] <= 0):
            break
    
    for row in range(num_rows):
    
        #SEM Base Bleed:
        k1 = dt*CdReMa0(tsb,xsb)[0];
        k2 = dt*CdReMa0(tsb+0.5*dt,xsb+0.5*k1)[0];
        k3 = dt*CdReMa0(tsb+0.5*dt,xsb+0.5*k2)[0];
        k4 = dt*CdReMa0(tsb+1.0*dt,xsb+1.0*k3)[0]; 
        
        xsb = xsb+(1/6)*(k1+2*(k2+k3)+k4);
        tsb = tsb+dt;
        Psb = CdReMa0(tsb,xsb)[2];
        Rsb = CdReMa0(tsb,xsb)[1];
        T0.append(tsb)
        Ysb = np.append(Ysb,xsb,axis=0)
        PRTsb = np.append(PRTsb,Psb,axis=0)
        ReMasb = np.append(ReMasb,Rsb,axis=0)
        
        if (Ysb[row+1,1] <= 0):
            break        

    for row in range(num_rows):
        
        #COM Base Bleed
        l1 = dt*CdReMa1(tbb,xbb)[0];
        l2 = dt*CdReMa1(tbb+0.5*dt,xbb+0.5*l1)[0];
        l3 = dt*CdReMa1(tbb+0.5*dt,xbb+0.5*l2)[0];
        l4 = dt*CdReMa1(tbb+1.0*dt,xbb+1.0*l3)[0];
        
        xbb = xbb+1/6*(l1+2*(l2+l3)+l4);
        tbb = tbb+dt;
        Pbb = CdReMa1(tbb,xbb)[2];
        Rbb = CdReMa1(tbb,xbb)[1];
        T1.append(tbb)
        Ybb = np.append(Ybb,xbb,axis=0)
        PRTbb = np.append(PRTbb,Pbb,axis=0)
        ReMabb = np.append(ReMabb,Rbb,axis=0)
    
        if (Ybb[row+1,1] <= 0):
            break
                  
    return Ysb, Ybb, Y, T, T0, T1, ReMasb, ReMabb, PRTsb, PRTbb

'''Construir os graficos:
Ysb, Ybb, Y, T, T0, T1, ReMasb, ReMabb, PRTsb, PRTbb = trajectory()
       
mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(Ysb[:,0], Ysb[:,2],Ysb[:,1], label='sem BB')
ax.plot(Ybb[:,0], Ybb[:,2],Ybb[:,1], label='com BB')
ax.set_xlabel('[metros]')
ax.set_ylabel('[metros]')
ax.set_zlabel('[metros]')
ax.legend()
plt.title("Alcance x Apogeu")
plt.show()'''




    



