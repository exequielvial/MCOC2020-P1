# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:15:42 2020

@author: Speedy Gonzàlez
"""

import numpy as np
from scipy.integrate import odeint
import  matplotlib.pylab as plt
import xml
import xml.etree.ElementTree as ET
from numpy import zeros
import datetime as dt

G = 6.67 * (10**(-11))  #nm^2/kg^2
km = 1000 #metros

minutos = 60
hora = minutos* 60
omega = 2*np.pi/(24*3600) #rad/s
time=24*3600
masaTierra = 5.972 * 10 ** 24 
radioAtmosfera = (6371+80)*km
R = np.eye(3,dtype = int)
dR_dt = np.eye(3,dtype = int)
dR2_dt2 = np.eye(3,dtype = int)

#omega = 2*sp.pi/86400
ecMovimiento = np.zeros(3)

def utc2time(utc, ut1, EOF_datetime_format = "%Y-%m-%dT%H:%M:%S.%f"):
	t1 = dt.datetime.strptime(ut1,EOF_datetime_format)
	t2 = dt.datetime.strptime(utc,EOF_datetime_format)
	return (t2 - t1).total_seconds()

def leer_eof(fname):
	tree = ET.parse(fname)
	root = tree.getroot()

	Data_Block = root.find("Data_Block")		
	List_of_OSVs = Data_Block.find("List_of_OSVs")

	count = int(List_of_OSVs.attrib["count"])

	t = zeros(count)
	x = zeros(count)
	y = zeros(count)
	z = zeros(count)
	vx = zeros(count)
	vy = zeros(count)
	vz = zeros(count)

	set_ut1 = False
	for i, osv in enumerate(List_of_OSVs):
		UTC = osv.find("UTC").text[4:]
		
		x[i] = osv.find("X").text   #conversion de string a double es implicita
		y[i] = osv.find("Y").text
		z[i] = osv.find("Z").text
		vx[i] = osv.find("VX").text
		vy[i] = osv.find("VY").text
		vz[i] = osv.find("VZ").text

		if not set_ut1:
			ut1 = UTC
			set_ut1 = True

		t[i] = utc2time(UTC, ut1)

	return t, x, y, z, vx, vy, vz

def satelite(z,t):
    
    R=np.array([[np.cos(omega*t),-np.sin(omega*t),0],[np.sin(omega*t),np.cos(omega*t), 0],[0.,0.,1]])
    dR_dt=np.array([[-np.sin(omega*t),-np.cos(omega*t),0],[np.cos(omega*t),-np.sin(omega*t), 0],[0.,0.,0.]])*omega
    dR2_dt2=np.array([[-np.cos(omega*t),np.sin(omega*t), 0],[-np.sin(omega*t),-np.cos(omega*t),0],[0.,0.,0]])*omega**2
    
    zp = np.zeros(6)
    z1 = z[0:3]
    r2 = np.dot(z1,z1)
    r = np.sqrt(r2)
    zp[0:3] = z[3:6]
    z2p=(-G*masaTierra/(r**3))*z[0:3]-R.T@(dR2_dt2@z[0:3]+2*dR_dt@z[3:6])
    zp[3:6] = z2p
    return zp




t, x, y, z, vx, vy, vz = leer_eof("S1B_OPER_AUX_POEORB_OPOD_20200811T110756_V20200721T225942_20200723T005942.EOF") 
z0 = np.array([x[0],y[0],z[0],vx[0],vy[0],vz[0]])

sol = odeint(satelite,z0,t)
x_sol = sol[:,0]
y_sol = sol[:,1]
z_sol = sol[:,2]

#parámetros gráficos
y1=[-5e6,0,5e6];ey=["-5000","0","5000"];x1=[0,18000,36000,54000,72000,90000];ex=["0","5","10","15","20","25"]


plt.figure()
plt.subplot(3,1,1)
plt.plot(t,x)
plt.plot(t,x_sol)
plt.title("Posición")
plt.ylabel("X (KM)")
plt.yticks(y1,ey)
plt.xticks(x1,ex)


plt.subplot(3,1,2)
plt.plot(t,y)
plt.plot(t,y_sol)
plt.ylabel("Y (KM)")
plt.yticks(y1,ey)
plt.xticks(x1,ex)


plt.subplot(3,1,3)
plt.plot(t,z)
plt.plot(t,z_sol)
plt.ylabel("Z (KM)")
plt.yticks(y1,ey)
plt.xticks(x1,ex)
plt.xlabel('Tiempo, t (horas)')


plt.tight_layout()


