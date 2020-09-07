from scipy.integrate import odeint
import  matplotlib.pylab as plt
import xml.etree.ElementTree as ET
from numpy import zeros
import datetime as dt
import numpy as np



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
    sin = np.sin(omega*t)
    cos = np.cos(omega*t)
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
    FJ2 = np.array([Fx_J2,Fy_J2,Fz_J2])
    
    Fx_J3 = J3*u*w/r**9*(10*w**2-15/2*(u**2+v**2))
    Fy_J3 = J3*v*w/r**9*(10*w**2-15/2*(u**2+v**2))
    Fz_J3 = J3/r**9*(4*w**2*(w**2-3*(u**2+v**2))+3/2*(u**2+v**2)**2)
    FJ3 = np.array([Fx_J3,Fy_J3,Fz_J3])
    

    R=np.array([[np.cos(omega*t),-np.sin(omega*t),0],[np.sin(omega*t),np.cos(omega*t), 0],[0.,0.,1]])
    dR_dt=np.array([[-np.sin(omega*t),-np.cos(omega*t),0],[np.cos(omega*t),-np.sin(omega*t), 0],[0.,0.,0.]])*omega
    dR2_dt2=np.array([[-np.cos(omega*t),np.sin(omega*t), 0],[-np.sin(omega*t),-np.cos(omega*t),0],[0.,0.,0]])*omega**2
    

    zp[0:3] = z[3:6] 
    zp[3:6]= (-G*masaTierra/r**3)*z[0:3]  - R.T@(dR2_dt2@z[0:3] + 2*dR_dt@z[3:6]) +FJ2[0:3]+FJ3[0:3]
    return zp


def eulerint(zp,z0,t,Nsubdivisiones = 1):
    Nt = len(t)
    Ndim = len(np.array(z0))
    z = np.zeros((Nt,Ndim))
    z[0,:] = z0[:]
    for i in range(1,Nt):
        t_anterior = t[i-1]
        dt = (t[i] - t[i-1])/Nsubdivisiones
        z_temp = z[i-1,:].copy()
        for k in range(Nsubdivisiones):
            z_temp += dt * satelite(z_temp,t_anterior + k*dt)
        z[i,:] = z_temp
    return z


t, x, y, z, vx, vy, vz = leer_eof('S1B_OPER_AUX_POEORB_OPOD_20200811T110756_V20200721T225942_20200723T005942.EOF') 
z0 = np.array([x[0],y[0],z[0],vx[0],vy[0],vz[0]]) 


sol = odeint(satelite,z0,t)
x_sol = sol[:,0]
y_sol = sol[:,1]
z_sol = sol[:,2]

zp = sol[:,:] 
sol2 = eulerint(zp,z0,t,Nsubdivisiones = 1)
x_euler = sol2[:,0]
y_euler = sol2[:,1]
z_euler = sol2[:,2]
delta = np.sqrt((x-x_sol)**2+(y-y_sol)**2+(z-z_sol)**2)
delta2 = np.sqrt((x-x_euler)**2+(y-y_euler)**2+(z-z_euler)**2)
plt.plot(t/3600,delta2/1000)
plt.plot(t/3600,delta/1000)
plt.ylabel('Delta (km)')
plt.xlabel('Tiempo, t (horas)')
plt.show()

