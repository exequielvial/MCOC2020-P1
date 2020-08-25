# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:59:48 2020

@author: Speedy Gonzàlez
"""

import scipy as sp 
from scipy.integrate import odeint
import matplotlib.pylab as plt

r_aire = 1.225  #kg/m**3
cdrag = 0.47
cm = 0.01
inch=2.54*cm
g=9.81
m = 15

D_pelota=8.5*inch
radio = D_pelota/2
A = sp.pi*radio**2
CD = 0.5*r_aire*cdrag*A

#viento
vs = [0,10.0,20.0]

#función a integrar
#z es el vecto de estado
#z = [x,y,vx,vy]
#dz/dt = bala (z,t)
#         [    z2    ]
# dz/dt = [          ] (modelo)
#         [  Fd/m  -g]

#ds/dt = []

#vector de estado

def bala(z,t):
    zp = sp.zeros(4)
    zp[0] = z[2]
    zp[1] = z[3]     
    v = z[2:4]  #saca los últimos 2 componentes
    v[0] = v[0] -V
    vnorm = sp.sqrt(sp.dot(v,v))
    FD = -CD*sp.dot(v,v)*(v/vnorm)
    zp[2] = FD[0]/m
    zp[3] = FD[1]/m -g
    
    return zp

t = sp.linspace(0,30,1001)

#parte en el origen y tiene vx = vy = 2m/s
vi= 100*1000/3600
z0= sp.array([0,0,vi,vi])

plt.figure(1)
xs = []
ys = []
for V in vs:
    sol = odeint(bala,z0,t)
    x = sol[:,0]
    y = sol[:,1]
    plt.ylim(0,50)
    plt.xlim(0,150)
    V = str(V)
    plt.plot(x,y, label = 'V = '+ V +' m/s')
    
plt.grid(True)
un = 'm/s'
plt.title('Trayectoria para distintos vientos')
plt.legend()
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.savefig('Grafico_balistica.png')




