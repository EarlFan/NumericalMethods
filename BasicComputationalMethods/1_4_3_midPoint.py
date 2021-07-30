# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 21:50 2021

@author: FanE

u = u(t)
u_t = f = u*u, u(0) = 1

https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-90-computational-methods-in-aerospace-engineering-spring-2014/numerical-integration-of-ordinary-differential-equations/convergence/1690r-rate-of-convergence--global-order-of-accuracy-/
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

dfdt = 0
u = 1
T =10
N1 = 250
dt = T/N1
results1 = np.zeros((3,N1),dtype=float)
i = 0
# forward Euler
while (i<N1-1):
    print("i = %d.\n"%(i))

    if(i==0):
        f = -u*u
        # update other variables
        results1[0,i]=dt*i
        results1[1,i]=u

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
        u = results1[1,i-1]+2*dt*f
    
    i = i+1
    # update other variables
    # NOTES: to avoid the problem of negative Re, use abs(u) instead
    f = -u*u
    results1[0,i]=dt*i
    results1[1,i]=u


u = 1
N1 = 500
dt = T/N1
results2 = np.zeros((3,N1),dtype=float)
i = 0
# forward Euler
while (i+1<N1):
    print("i = %d.\n"%(i))

    if(i==0):
        f = -u*u
        # update other variables
        results2[0,i]=dt*i
        results2[1,i]=u

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
        u = results2[1,i-1]+2*dt*f
    
    i = i+1
    # update other variables
    # NOTES: to avoid the problem of negative Re, use abs(u) instead
    f = -u*u
    results2[0,i]=dt*i
    results2[1,i]=u

u = 1
N1 = 1000
dt = T/N1
results3 = np.zeros((3,N1),dtype=float)
i = 0
# forward Euler
while (i+1<N1):
    print("i = %d.\n"%(i))

    if(i==0):
        f = -u*u
        # update other variables
        results3[0,i]=dt*i
        results3[1,i]=u

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
        u = results3[1,i-1]+2*dt*f
    
    i = i+1
    # update other variables
    # NOTES: to avoid the problem of negative Re, use abs(u) instead
    f = -u*u
    results3[0,i]=dt*i
    results3[1,i]=u

results1[2,:]=abs(results1[1,:]-1/(results1[0,:]+1))
results2[2,:]=abs(results2[1,:]-1/(results2[0,:]+1))
results3[2,:]=abs(results3[1,:]-1/(results3[0,:]+1))

# plt.plot(results1[0,:], results1[1,:])       # Plot the sine of each x point
# plt.show()                   # Display the plot

with plt.style.context(['science', 'ieee','grid']):
    fig, ax = plt.subplots(2)

    m_markevery=5
    m_markersize=1

    fig_name = ['Velocity','Re','Cd']
    label_format = '{:,.0e}'

    fig_index = 0
    ax[0].plot(results1[0,:], results1[1,:],label="dt=0.04",color='b',marker='None',markersize=m_markersize,markevery=m_markevery)
    ax[0].plot(results2[0,:], results2[1,:],label="dt=0.02",color='r',marker='None',markersize=m_markersize,markevery=m_markevery)
    ax[0].plot(results3[0,:], results3[1,:],label="dt=0.01",color='g',marker='None',markersize=m_markersize,markevery=m_markevery)
    ax[0].plot(results3[0,:], 1/(1+results3[0,:]),label="analytical",marker='o',color='k',linestyle="None",markersize=1,markevery=10)
    ax[0].legend(bbox_to_anchor=(1.05, 0.8, 0.3, 0.2), loc='upper left')
    ax[0].autoscale(tight=True)
    pparam = dict(ylabel='$u$',xlabel='$t$')
    ax[0].set_xlim([0,10])
    ax[0].set_ylim([0,1])
    ax[0].set(**pparam)

    ax[1].plot(results1[0,:], results1[2,:],label="dt=0.04",color='b',marker='None',markersize=m_markersize,markevery=m_markevery)
    ax[1].plot(results2[0,:], results2[2,:],label="dt=0.02",color='r',marker='None',markersize=m_markersize,markevery=m_markevery)
    ax[1].plot(results3[0,:], results3[2,:],label="dt=0.01",color='g',marker='None',markersize=m_markersize,markevery=m_markevery)
    # ax[1].plot(results3[0,:], 1/(1+results3[0,:]),label="analytical",marker='*',linestyle="None",markersize=3,markevery=1)
    ax[1].legend(bbox_to_anchor=(1.05, 0.8, 0.3, 0.2), loc='upper left')
    ax[1].autoscale(tight=True)
    pparam = dict(ylabel="error",xlabel='$t$')
    ax[1].set_xlim([0,10])
    # ax[1].set_ylim([0,1])
    ax[1].set(**pparam)


    # plt.show()
    # fig.savefig('figures/fig1.pdf')
    fig.suptitle('Middle Point, $u_t = -u^2$')
    fig.savefig('1_4_3-MidPoint.jpg', dpi=600)