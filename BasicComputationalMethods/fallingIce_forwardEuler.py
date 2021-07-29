# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:33:49 2021

@author: FanE

The falling ice problem using the forward Euler method.

https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-90-computational-methods-in-aerospace-engineering-spring-2014/numerical-integration-of-ordinary-differential-equations/discretizing-odes/1690r-the-forward-euler-method/
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


# particle radius
a = 0.01 #m
# particle density
rho_p = 917 # kg/m3
# particle velocity
m_p = rho_p*math.pi*(a**3)*4/3

# gas properties
rho_g = 0.9
mu_g = 1.69e-5
g = 9.8

# initial state
u0 = 1e-4
u = u0

# construct the driving force
# Re = 2*rho_g*u*a/mu_g
# Cd = 24/Re + 6/(1+math.sqrt(Re))+0.4
# D = 0.5*rho_g*math.pi*a^2*u^2*Cd
# f = g - D/m_p

# construct the df/du, df/dt

# df/du = - (dD/dt)/m_p
# dRedu = 2*rho_g*a/mu_g
# dDdRe = -24/(Re^2)-6/((1+math.sqrt(Re))^2)*0.5*(Re^(-0.5))
# dDdu = 0.5*rho_g*math.pi*a^2*[2*u*Cd+u^2*dDdRe*dRedu]
# dfdu = -dDdu/m_p

dfdt = 0

dt = 0.25
N = 100
results = np.zeros((4,N),dtype=float)
i = 0
# forward Euler
while (i<N):
    
    Re = 2*rho_g*u*a/mu_g
    Cd = 24/Re + 6/(1+math.sqrt(Re))+0.4
    D = 0.5*rho_g*math.pi*(a**2)*u*u*Cd
    f = g - D/m_p

    results[0,i]=dt*i
    results[1,i]=u
    results[2,i]=Re
    results[3,i]=Cd
    if(i==0): results[3,i]=0
    
    u = u + f*dt
    i = i+1
    
# plt.plot(results[0,:], results[1,:])       # Plot the sine of each x point
# plt.show()                   # Display the plot



x = np.linspace(0.75, 1.25, 201)

with plt.style.context(['science', 'ieee','grid']):
    fig, ax = plt.subplots(3)

    m_markevery=5
    m_markersize=1

    fig_name = ['Velocity','Re','Cd']

    fig_index = 0
    ax[fig_index].plot(results[0,:], results[1,:],marker='None',markersize=m_markersize,markevery=m_markevery)
    # ax[fig_index].legend(title='Forward Euler Method')
    ax[fig_index].autoscale(tight=True)
    pparam = dict(ylabel='$u (m/s)$')
    ax[fig_index].set_xlim([0,25])
    ax[fig_index].set_ylim([0,30])
    ax[fig_index].set(**pparam)
    ax[fig_index].get_xaxis().set_ticklabels([])

    fig_index = 1
    ax[fig_index].plot(results[0,:], results[2,:],marker='None',markersize=m_markersize,markevery=m_markevery)
    # ax[fig_index].legend(title=fig_name[fig_index])
    ax[fig_index].autoscale(tight=True)
    pparam = dict(ylabel='$Re$')
    ax[fig_index].set_xlim([0,25])
    ax[fig_index].set_ylim([0,3e4])
    ax[fig_index].set(**pparam)
    ax[fig_index].get_xaxis().set_ticklabels([])
    ax[fig_index].get_yaxis().set_ticklabels([1e4,2e4,3e4])
    ax[fig_index].yaxis.set_major_formatter(FormatStrFormatter('%.0e'))

    fig_index = 2
    ax[fig_index].plot(results[0,1:], results[3,1:],marker='None',markersize=m_markersize,markevery=m_markevery)
    # ax[fig_index].legend(title=fig_name[fig_index])
    ax[fig_index].autoscale(tight=True)
    pparam = dict(xlabel='$t (s)$', ylabel='$C_D$')
    ax[fig_index].set_xlim([0,25])
    ax[fig_index].set_ylim([0.4,0.6])
    ax[fig_index].set(**pparam)

    # plt.show()
    # fig.savefig('figures/fig1.pdf')
    fig.suptitle('Forward Euler Method')
    fig.savefig('ForwardEulerMethod.jpg', dpi=600)