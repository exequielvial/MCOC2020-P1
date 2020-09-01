# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 10:43:07 2020

@author: Speedy Gonzàlez
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pylab import *
from scipy.integrate import odeint
from eulerint import eulerint

# euler int es inestable

#condición inicial:
x0 = 1
x0punto = 1
z0 = [x0,x0punto]  #vector inicial

m = 1
f = 1
eRara = 0.2
w = 2*np.pi*f
k = m*w**2
c = 2*eRara*w*m

def zpunto(z,t):
    zp = np.zeros(2)
    zp[0] = z[1]
    zp[1] = -(c*z[1]+k*z[0])/m

    return zp


t = linspace(0,4.,100)

sol = odeint(zpunto,z0, t)
z_odeint=(sol[:,0])        #odeint
plt.plot(t,z_odeint,'b', label='Odeint')

z_real = np.exp(-c*t/2)*np.cos(w*t)
plt.plot(t, z_real, 'k', linewidth = 2, label = 'Solución analítica')  #solución analítica
sol = eulerint(zpunto, z0,t,Nsubdivisiones= 1)  #1
z_euler = sol[:,0]
plt.plot(t,z_euler,':',color ='g',label='1 subdv')
sol = eulerint(zpunto, z0,t,Nsubdivisiones= 10)  #10
z_euler = sol[:,0]
plt.plot(t,z_euler,':',color ='r',label='10 subdv')
sol = eulerint(zpunto, z0,t,Nsubdivisiones= 100)  # 100
z_euler = sol[:,0]
plt.plot(t,z_euler,':',color ='orange',label='100 subdv')

plt.title('Oscilador armónico')
plt.ylabel('X(t)')
plt.xlabel('Tiempo')
legend()
show()