# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 23:03:53 2020

@author: Speedy Gonzàlez
"""


import numpy as np
from scipy.integrate import odeint
import  matplotlib.pylab as plt
from sys import argv
import xml
import xml.etree.ElementTree as ET
from numpy import zeros
import datetime as dt

G = 6.67 * (10**(-11))  #nm^2/kg^2
km = 1000 #metros
J2 = 1.75553*(10**10)*1000**5
J3 = -2.61913*(10**11)*1000**6
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
    #parámetros
    zp = np.zeros(6)
    #sin = np.sin(omega*t)
    #cos = np.cos(omega*t)
    z1 = z[0:3]
    r2 = np.dot(z1,z1)
    r = np.sqrt(r2)
    J2 = 1.75553*(10**10)*1000**5
    J3 = -2.61913*(10**11)*1000**6
    u = z[0]#x
    v = z[1]#y
    w = z[2]#z
    
    
    Fx_J2 = J2*(u/r**7)*(6*w**2-3/2*(u**2+v**2))
    Fy_J2 = J2*(v/r**7)*(6*w**2-3/2*(u**2+v**2))
    Fz_J2 = J2*(w/r**7)*(3*w**2-9/2*(u**2+v**2))

    
    Fx_J3 = J3*u*w/r**9*(10*w**2-15/2*(u**2+v**2))
    Fy_J3 = J3*v*w/r**9*(10*w**2-15/2*(u**2+v**2))
    Fz_J3 = J3/r**9*(4*w**2*(w**2-3*(u**2+v**2))+3/2*(u**2+v**2)**2)

    

    R=np.array([[np.cos(omega*t),-np.sin(omega*t),0],[np.sin(omega*t),np.cos(omega*t), 0],[0.,0.,1]])
    dR_dt=np.array([[-np.sin(omega*t),-np.cos(omega*t),0],[np.cos(omega*t),-np.sin(omega*t), 0],[0.,0.,0.]])*omega
    dR2_dt2=np.array([[-np.cos(omega*t),np.sin(omega*t), 0],[-np.sin(omega*t),-np.cos(omega*t),0],[0.,0.,0]])*omega**2
    
    fg = (-G*masaTierra/r**2)*(R@(z1/r))
    zp[3:6]=R.T@(fg-(2*(dR_dt@z[3:6])+(dR2_dt2@z[0:3])))
    
    zp[0:3] = z[3:6] 

    
    zp[3] = zp[3] + Fx_J2 + Fx_J3
    zp[4] = zp[4] + Fy_J2 + Fy_J3
    zp[5] = zp[5] + Fz_J2 + Fz_J3
    
    return zp

eof_out = argv[1]
# SI NO FUNCIONA DESDE LA CONSOLA DESCOMENTAR LA FILA DE ABAJO Y COMENTAR LA DE ARRIBA, GRACIAS!
#eof_out = "S1B_OPER_AUX_POEORB_OPOD_20200811T110756_V20200721T225942_20200723T005942.EOF"

pred = eof_out.replace('.EOF','.PRED')
t, x, y, z, vx, vy, vz = leer_eof(eof_out) 
z0 = np.array([x[0],y[0],z[0],vx[0],vy[0],vz[0]])
zf= np.array([x[-1],y[-1],z[-1],vx[-1],vy[-1],vz[-1]])

sol = odeint(satelite,z0,t)
x_sol = sol[:,0];y_sol = sol[:,1];z_sol = sol[:,2]
vx_sol = sol[:,3]; vy_sol = sol[:,4]; vz_sol = sol[:,5]

with open(eof_out,'w') as fout:
    fout.write('<?xml version="1.0" ?>\n'
'<Earth_Explorer_File>\n'
'  <Earth_Explorer_Header>\n'
'    <Fixed_Header>\n'
'      <File_Name>S1A_OPER_AUX_POEORB_OPOD_20200816T120754_V20200726T225942_20200728T005942</File_Name>\n'
'      <File_Description>Precise Orbit Ephemerides (POE) Orbit File</File_Description>\n'
'      <Notes></Notes>\n'
'      <Mission>Sentinel-1A</Mission>\n'
'      <File_Class>OPER</File_Class>\n'
'      <File_Type>AUX_POEORB</File_Type>\n'
'      <Validity_Period>\n'
'        <Validity_Start>UTC=2020-07-26T22:59:42</Validity_Start>\n'
'        <Validity_Stop>UTC=2020-07-28T00:59:42</Validity_Stop>\n'
'      </Validity_Period>\n'
'      <File_Version>0001</File_Version>\n'
'      <Source>\n'
'        <System>OPOD</System>\n'
'        <Creator>OPOD</Creator>\n'
'        <Creator_Version>0.0</Creator_Version>\n'
'        <Creation_Date>UTC=2020-08-16T12:07:54</Creation_Date>\n'
'      </Source>\n'
'    </Fixed_Header>\n'
'    <Variable_Header>\n'
'      <Ref_Frame>EARTH_FIXED</Ref_Frame>\n'
'      <Time_Reference>UTC</Time_Reference>\n'
'    </Variable_Header>\n'
'  </Earth_Explorer_Header>\n'
'<Data_Block type="xml">\n'
'  <List_of_OSVs count="9361">\n')
    Nt=len(t) 
    for i in range(Nt):
        Dia1 = dt.datetime(2020,7,26,23,00,29,000000)
        Dia2 = dt.datetime(2020,7,26,22,59,52,000000)
        Dia3 = dt.datetime(2020,7,26,22,59,51,787522)
        dias1 = (Dia1 + dt.timedelta(seconds=t[i])).strftime('%Y-%m-%dT%H:%M:%S.%f')
        dias2 = (Dia2 + dt.timedelta(seconds=t[i])).strftime('%Y-%m-%dT%H:%M:%S.%f')
        dias3 = (Dia3 + dt.timedelta(seconds=t[i])).strftime('%Y-%m-%dT%H:%M:%S.%f')
        fout.write('    <OSV>\n'
f'      <TAI>TAI={dias1}</TAI>\n'
f'      <UTC>UTC={dias2}</UTC>\n'
f'      <UT1>UT1={dias3}</UT1>\n'
f'      <Absolute_Orbit>+226632</Absolute_Orbit>\n'
f'      <X unit="m">{x_sol[i]}</X>\n'
f'      <Y unit="m">{y_sol[i]}</Y>\n'
f'      <Z unit="m">{z_sol[i]}</Z>\n'
f'      <VX unit="m/s">{vx_sol[i]}</VX>\n'
f'      <VY unit="m/s">{vy_sol[i]}</VY>\n'
f'      <VZ unit="m/s">{vz_sol[i]}</VZ>\n'
'      <Quality>NOMINAL</Quality>\n'
'    </OSV>\n')

'</Data_Block>\n'
'</Earth_Explorer_File>\n'
