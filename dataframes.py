# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 17:24:49 2020

@author: Wallace
"""
from tcc_py import trajectory
from pandas import DataFrame as df

#Construir as tabelas:
Ysb, Ybb, Y, T, T0, T1, ReMasb, ReMabb, PRTsb, PRTbb = trajectory() 

dsb = {'t': T0[:],'X': Ysb[:,0],'Y': Ysb[:,1],'Z': Ysb[:,2],'Vx': Ysb[:,3],'Vy': Ysb[:,4],'Vz': Ysb[:,5],'Cd': ReMasb[:,0],'Mach': ReMasb[:,1],'Rey': ReMasb[:,2],'P': PRTsb[:,0],'Rho': PRTsb[:,1],'T': PRTsb[:,2]}
dbb = {'t': T1[:],'X': Ybb[:,0],'Y': Ybb[:,1],'Z': Ybb[:,2],'Vx': Ybb[:,3],'Vy': Ybb[:,4],'Vz': Ybb[:,5],'m': Ybb[:,6],'Cd': ReMabb[:,0],'Mach': ReMabb[:,1],'Rey': ReMabb[:,2],'P': PRTbb[:,0],'Rho': PRTbb[:,1],'T': PRTbb[:,2]}

data_sb = df.from_dict(dsb)
data_bb = df.from_dict(dbb)
 
data_sb.to_csv('out_nobb.csv',index=False, encoding='utf-8')
data_bb.to_csv('out_bb.csv',index=False, encoding='utf-8')
