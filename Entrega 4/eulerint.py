# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 11:06:45 2020

@author: Speedy Gonzàlez
"""
from matplotlib.pylab import *
from scipy.integrate import odeint

def eulerint(zp,z0,t,Nsubdivisiones=1):
    Nt = len(t)
    Ndim = len(np.array(z0))
    
    z = zeros((Nt,Ndim))
    z[0,:] = z0[0]
    
    #z (i+1) = zp_1 * dt + z_i
    for i in range(1,Nt):
        t_anterior = t[i-1]
        dt = (t[i] - t[i-1])/Nsubdivisiones
        
        z_temp = z[i-1, :].copy()
        for k in range(Nsubdivisiones):
            z_temp += dt* zp(z_temp,t_anterior+k*dt)
            
        z[i,:] = z_temp
        
    return z