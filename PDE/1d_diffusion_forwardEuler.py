# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:33:49 2021

@author: FanE

The 1-d heat diffusion solved by forward Euler method

http://hplgit.github.io/num-methods-for-PDEs/doc/pub/diffu/sphinx/._main_diffu001.html
"""
import numpy as np
import math

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import scipy

def I(x,type,L):
    """Plug profile as initial condition."""
    if(type == "plug"):
        # print("plug")
        if abs(x-L/2.0) > 0.1:
            return 0
        else:
            return 1
    elif(type == "gaussian"):
        # print("gaussian")
        sigma=0.05
        return math.exp(-0.5*((x-L/2.0)**2)/sigma**2)
    elif(type == "expsin"):
        # print("expsin")
        L = 10.0
        a = 1
        T = 1.2
        pi = math.pi
        return exp(-m**2*pi**2*a/L**2*t)*sin(m*pi/L*x)


def solver_FE_simple(a, L, dt, F, T,type):
    """
    Simplest expression of the computational algorithm
    using the Forward Euler method and explicit Python loops.
    For this method F <= 0.5 for stability.
    """
    import time;  t0 = time.clock()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)
    u_n = np.zeros(Nx+1)

    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u_n[i] = I(x[i],type,L)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_n[i] + F*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) 
        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

        if(n%(50)==0 or n==Nt-1):
            with plt.style.context(['science', 'ieee','grid']):
                fig, ax = plt.subplots()

                m_markevery=5
                m_markersize=1

                fig_index = 0
                if(n==0):
                    ax.plot(x,u,marker='None',markersize=m_markersize,markevery=m_markevery)
                else:
                    ax.plot(x,u_n,marker='None',markersize=m_markersize,markevery=m_markevery)
                # ax.legend(title='Forward Euler Method')
                ax.autoscale(tight=True)
                pparam = dict(ylabel='$T/T_{0,max}$')
                ax.set_xlim([0,1])
                ax.set_ylim([-0.1,1.1])
                ax.set(**pparam)
                # ax.get_xaxis().set_ticklabels([])
                fig.suptitle('t = %.2f, F = %.2f, FE.'%(n*dt,F))
                fileName = "%s_F_%.2f_FE_%.2e.jpg"%(type,F,n*dt)
                fig.savefig(fileName, dpi=600)

    t1 = time.clock()
    return u_n, x, t, t1-t0  # u_n holds latest u

def solver_BE_simple(a, L, dt, F, T,type):
    """
    Simplest expression of the computational algorithm
    for the Backward Euler method, using explicit Python loops
    and a dense matrix format for the coefficient matrix.
    """
    import time;  t0 = time.clock()  # for measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)
    u_n = np.zeros(Nx+1)

    # Data structures for the linear system
    A = np.zeros((Nx+1, Nx+1))
    b = np.zeros(Nx+1)

    for i in range(1, Nx):
        A[i,i-1] = -F
        A[i,i+1] = -F
        A[i,i] = 1 + 2*F
    A[0,0] = A[Nx,Nx] = 1
    A[0,1] = A[Nx,Nx-1] = 0

    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u_n[i] = I(x[i],type,L)

    for n in range(0, Nt):
        # Compute b and solve linear system
        for i in range(1, Nx):
            b[i] = u_n[i] 
        b[0] = b[Nx] = 0
        # Note: solve the Ax = b equation
        u[:] = np.linalg.solve(A, b)

        # Update u_n before next step
        u_n, u = u, u_n

        if(n%(50)==0 or n==Nt-1):
            with plt.style.context(['science', 'ieee','grid']):
                fig, ax = plt.subplots()

                m_markevery=5
                m_markersize=1

                fig_index = 0
                if(n==0):
                    ax.plot(x,u,marker='None',markersize=m_markersize,markevery=m_markevery)
                else:
                    ax.plot(x,u_n,marker='None',markersize=m_markersize,markevery=m_markevery)
                # ax.legend(title='Forward Euler Method')
                ax.autoscale(tight=True)
                pparam = dict(ylabel='$T/T_{0,max}$')
                ax.set_xlim([0,1])
                ax.set_ylim([-0.1,1.1])
                ax.set(**pparam)
                # ax.get_xaxis().set_ticklabels([])
                fig.suptitle('t = %.2f, F = %.2f, BE.'%(n*dt,F))
                fileName = "%s_F_%.2f_BE_%.2e.jpg"%(type,F,n*dt)
                fig.savefig(fileName, dpi=600)

    t1 = time.clock()
    return t1-t0

# cannot be use, scipy is not correctly installed
def solver_BE(a, L, dt, F, T,type):
    """
    Simplest expression of the computational algorithm
    for the Backward Euler method, using explicit Python loops
    and a dense matrix format for the coefficient matrix.
    """
    import time;  t0 = time.clock()  # for measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)
    u_n = np.zeros(Nx+1)

    # Data structures for the sparse matrix
    main = np.zeros(Nx+1)
    lower = np.zeros(Nx)
    upper = np.zeros(Nx)
    b = np.zeros(Nx+1)

    main[:] = 1+2*F
    lower[:] = -F
    upper[:] = -F

    # insert boundary condition
    main[0] = main[Nx] = 1

    A = scipy.sparse.diags(
        diagonals=[main, lower, upper],
        offsets=[0, -1, 1], shape=(Nx+1, Nx+1),
        format='csr')

    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u_n[i] = I(x[i],type,L)

    for n in range(0, Nt):
        # Compute b and solve linear system
        for i in range(1, Nx):
            b[i] = u_n[i] 
        b[0] = b[Nx] = 0
        # Note: solve the Ax = b equation
        u[:] = scipy.sparse.linalg.spsolve(A, b)

        # Update u_n before next step
        u_n, u = u, u_n

        if(n%(50)==0 or n==Nt-1):
            with plt.style.context(['science', 'ieee','grid']):
                fig, ax = plt.subplots()

                m_markevery=5
                m_markersize=1

                fig_index = 0
                if(n==0):
                    ax.plot(x,u,marker='None',markersize=m_markersize,markevery=m_markevery)
                else:
                    ax.plot(x,u_n,marker='None',markersize=m_markersize,markevery=m_markevery)
                # ax.legend(title='Forward Euler Method')
                ax.autoscale(tight=True)
                pparam = dict(ylabel='$T/T_{0,max}$')
                ax.set_xlim([0,1])
                ax.set_ylim([-0.1,1.1])
                ax.set(**pparam)
                # ax.get_xaxis().set_ticklabels([])
                fig.suptitle('t = %.2f, F = %.2f, BE.'%(n*dt,F))
                fileName = "%s_F_%.2f_BE_%.2e.jpg"%(type,F,n*dt)
                fig.savefig(fileName, dpi=600)

    t1 = time.clock()
    return t1-t0

# plug test
L = 1
a = 1
T = 0.1
F = 0.5
Nx = 50

dx = L/Nx
dt = F/a*dx**2


# solver_FE_simple(a, L, dt, 0.51, T,"gaussian")

solver_BE(a, L, dt, 1, T,"plug")

