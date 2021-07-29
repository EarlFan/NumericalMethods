# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:48:49 2021

@author: FanE

The falling ice problem using the middle point method.

https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-90-computational-methods-in-aerospace-engineering-spring-2014/numerical-integration-of-ordinary-differential-equations/discretizing-odes/1690r-the-midpoint-method/
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

T_total=25
N = 1000
dt = T_total/N

results = np.zeros((4,N),dtype=float)
i = 0

def cal_f(u):
    Re = 2*rho_g*u*a/mu_g
    Cd = 24/Re + 6/(1+math.sqrt(Re))+0.4
    D = 0.5*rho_g*math.pi*(a**2)*u*u*Cd
    f = g - D/m_p
    return f

# middle point method 
# dudt(n) = 0.5*(u(n+1)-u(n-1))/dt                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
while (i<N-1):
    print("i = %d.\n"%(i))

    if(i==0):
        Re = 2*rho_g*u*a/mu_g
        Cd = 24/Re + 6/(1+math.sqrt(Re))+0.4
        D = 0.5*rho_g*math.pi*(a**2)*u*u*Cd
        f = g - D/m_p
        # update other variables
        results[0,i]=dt*i
        results[1,i]=u
        results[2,i]=Re
        results[3,i]=Cd

    # update u

    # the first step cannot be estimated by middle point rule
    # a forward euler method is appied instead
    # u(n+1) = u(n)+dt*f(u(n))

    if (i+1<2):
        # the step of forward Euler
        u = u + f*dt

    # NOTES: be careful for the time step of the temporal advancement
    # write it in this format:
    # u(n+1) = u(n-1) + 2*dt*f(u(n))

    else:
        # the middle point rule
        u = results[1,i-1]+2*dt*f
    
    i = i+1
    # update other variables
    # NOTES: to avoid the problem of negative Re, use abs(u) instead
    Re = 2*rho_g*abs(u)*a/mu_g
    Cd = 24/Re + 6/(1+math.sqrt(Re))+0.4
    D = 0.5*rho_g*math.pi*(a**2)*u*u*Cd
    f = g - D/m_p
    results[0,i]=dt*i
    results[1,i]=u
    results[2,i]=Re
    results[3,i]=Cd
    
# plt.plot(results[0,:], results[1,:])       # Plot the sine of each x point
# plt.show()                   # Display the plot


with plt.style.context(['science', 'ieee','grid']):
    fig, ax = plt.subplots(3)

    m_markevery=5
    m_markersize=1

    fig_name = ['Velocity','Re','Cd']
    label_format = '{:,.0e}'

    fig_index = 0
    ax[fig_index].plot(results[0,:], results[1,:],marker='None',markersize=m_markersize,markevery=m_markevery)
    # ax[fig_index].legend(title='Forward Euler Method')
    ax[fig_index].autoscale(tight=True)
    pparam = dict(ylabel='$u (m/s)$')
    ax[fig_index].set_xlim([0,25])
    ax[fig_index].set_ylim([-100,50])
    ax[fig_index].set(**pparam)
    ax[fig_index].get_xaxis().set_ticklabels([])

    fig_index = 1
    ax[fig_index].plot(results[0,:], results[2,:],marker='None',markersize=m_markersize,markevery=m_markevery)
    # ax[fig_index].legend(title=fig_name[fig_index])
    ax[fig_index].autoscale(tight=True)
    pparam = dict(ylabel='$Re$')
    ax[fig_index].set_xlim([0,25])
    ax[fig_index].set_ylim([-1e5,1e5])
    ax[fig_index].set(**pparam)
    ax[fig_index].get_xaxis().set_ticklabels([])
    # ax[fig_index].get_yaxis().set_major_locator(mticker.FixedLocator(ticks_loc))
    ax[fig_index].yaxis.set_major_formatter(FormatStrFormatter('%.0e'))

    fig_index = 2
    ax[fig_index].plot(results[0,1:], results[3,1:],marker='None',markersize=m_markersize,markevery=m_markevery)
    # ax[fig_index].legend(title=fig_name[fig_index])
    ax[fig_index].autoscale(tight=True)
    pparam = dict(xlabel='$t (s)$', ylabel='$C_D$')
    ax[fig_index].set_xlim([0,25])
    ax[fig_index].set_ylim([0.3,0.7])
    ax[fig_index].set(**pparam)

    # plt.show()
    # fig.savefig('figures/fig1.pdf')
    fig.suptitle('Forward Euler Method')
    fig.savefig('1_2_5-MiddlePointMethod.jpg', dpi=600)