# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:15:42 2020

@author: Speedy Gonzàlez
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pylab as plt

G = 6.67 * (10**(-11))  #nm^2/kg^2
km = 1000 #metros
r = 7071*km
minutos = 60
hora = minutos* 60
rotacionTierra = 2*np.pi/(24*3600) #rad/s
time=24*3600
masaTierra = 5.972 * 10 ** 24 
radioAtmosfera = 6451*km
R = np.eye(3,dtype = int)
Rp = np.eye(3,dtype = int)
Rpp = np.eye(3,dtype = int)
zp = np.zeros(6)
ecMovimiento = np.zeros(3)
def satelite(z,t):
    
    ang = rotacionTierra*t
    R[0][0]= np.cos(ang); R[0][1] = -np.sin(ang); R[1][0] = -R[0][1]; R[1][1] = R[0][0]
    Rp[0][0]= (-np.sin(ang))*rotacionTierra; Rp[0][1] = (-np.cos(ang))*rotacionTierra; Rp[1][0] = -Rp[0][1]; Rp[1][1] = Rp[0][0];Rp[2][2]= 0
    Rpp[0][0]= (-np.cos(ang))*(rotacionTierra**2); Rpp[0][1] = (np.sin(ang))*(rotacionTierra**2); Rpp[1][0] = -Rpp[0][1]; R[1][1] = Rpp[0][0];Rpp[2][2]= 0
 
    
    
    zp[0:3] = z[3:6]
    ecMovimiento = (-G*masaTierra/r**3)*z[0:3] - R.T@(Rpp@z[0:3]+ 2*Rp@z[3:6])
    zp[3:6] = ecMovimiento 
    
    return zp

# posicion inicial
x = r; y = 0;z = 0
# velocidad
vx = 0; vy = 7507.54 ; vz = 0


t = np.linspace(0,time,1001)
z0 = np.array([x,y,z,vx,vy,vz])  #condición inicial


sol = odeint(satelite,z0,t)
x = sol[:,0]
y = sol[:,1]
z = sol[:,2]

#gráfico radio(t)
'''/plt.plot(t, np.sqrt(x**2+y**2+z**2))

plt.hlines(y=6371000,xmin=0,xmax=13000,color="r")
plt.hlines(y=6451000,xmin=0,xmax=13000,color="g")
plt.hlines(y=-6371000,xmin=0,xmax=13000,color="r")
plt.hlines(y=-6451000,xmin=0,xmax=13000,color="g")
plt.grid()
plt.legend(["r (t)","radio tierra","radio orbital"]) 
plt.xlabel('Tiempo (s)') ; plt.ylabel('Radio');plt.ylim(6000000);plt.savefig('Radio_fn_tiempo.png')/'''

#gráfico posición vs tiempo
'''/plt.plot(t,x);plt.plot(t,y);plt.plot(t,z)
plt.title('Posición en el tiempo')
plt.xlabel('Tiempo (s)');plt.ylabel('Posición')
plt.xlim(0)
plt.grid()
plt.savefig('graficoPosición.png')/'''

teta = np.linspace(0,2*np.pi,100) 
a = radioAtmosfera*np.cos(teta);b = radioAtmosfera*np.sin(teta)
figure,circle = plt.subplots(1)
circle.plot(a,b);circle.set_aspect(1)
plt.plot(x,y)    

# la velocidad de orbita está descrita por la siguiente expresión
v = (G*masaTierra/r)**0.5
print ('Velocidad de órbita: ', v)

plt.title('Órbita satélite Sentinel 1-A' )
plt.legend(['Atmósfera','V = 7505.54 m/s'])
plt.grid()
plt.savefig('Satélite_atmósfera')

plt.show

    
