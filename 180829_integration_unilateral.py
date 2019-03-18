#
#
#

from random import random
from math import pi, sin, cos, e, sqrt, tanh


# EXP = 1 # use exponential response function (1/tau)*e**(-t/tau)
# EXP = 0 # use power law response function  t**tau
# tau     # timescale, according to type (heavy tailed power law tau=-1.5)
# nu      # tropic sensitivity
# gamma   # propioceptive sensitivity

def unilateral_integration(stim, EXP=1, tau=5.0, nu=1.0, gamma=0.5, return_max=1):
    
    # constants
    L = 1.0    # length of the shoot
    N = 100    # number of bins dividing the shoot length
    T = len(stim)
    ds = L/N   # segment size
    dt = 0.005 # timestep
    
    if EXP==1:
        thist = 3.*tau # 5*tau
    if EXP==0:
        thist = 3.*abs(tau)

        
    nthist = int(thist/dt)
    stim = [0]*nthist + stim
        
        
    #
    # Memory parameters
    #
    # if exponential memory kernel
    #if EXP==1:
    #    mem_norm = sum([dt*((1/tau)*e**(-t*dt/tau)) for t in range(nthist)])
    #if EXP==0:
    #    mem_norm = sum([dt*(((t*dt)**(tau))) for t in range(1, nthist)])

        
    # Initial conditions:
    cu = [0.0]*N                  # cu(s, t=0) = d(theta)/dt = 0, all straight
    th = [0]*N                    # initial condition, theta(s, t=0)=0 all vertical
    th_mem = [[0]*N]*nthist        # it was with this same initial condition for a while
    shy = [0.0 for i in range(N)] # y position of each shoot segment
    shx = [s*ds*sin(th[s]) for s in range(N)] # x position of each shoot segment

    # Boundary conditions:
    th[0] = 0.0 # clamped vertically at the bottom theta(s=0, t)=0
    cu[0] = 0.0 # clamped at the bottom cu(s=0, t)=0

    theta_maxs = [] # save the maximal theta at the tip during the run.

    #
    # Euler integration in time
    #
    for t in range(0, T-1):
        if t%100==0:
            print t,
        # initialize new timesteps
        th_n = [0]*N # theta   
        cu_n = [0]*N # curvature
        shy_n = [0]*N # y position of each shoot segment
        shx_n = [0]*N # x position of each shoot segment

        for s in range(1, N):
            # integrate theta back in time
            if EXP==1: # exponential
                cu_n[s] = cu[s] - dt*(gamma*cu[s] + \
                                    sum([dt*(nu*stim[t-tt+nthist])*sin(th_mem[-tt][s] - pi/2)*\
                                            (e**(-(tt-1)*dt/tau))
                                        for tt in range(0, nthist-1)]) )
            if EXP==0: # power law
                cu_n[s] = cu[s] - dt*(gamma*cu[s] + \
                                    sum([dt*(nu*stim[t-tt+nthist])*sin(th_mem[-tt][s] - pi/2)*\
                                            ((dt*tt)**(tau))
                                        for tt in range(1, nthist)]))

            th_n[s] = ds*cu_n[s] + th_n[s-1]
            shy_n[s] = shy_n[s-1] + ds*cos(th_n[s-1])
            shx_n[s] = shx_n[s-1] + ds*sin(th_n[s]) # s*ds*sin(th_n[s])
        cu = cu_n[:]    
        th = th_n[:]    
        shy = shy_n[:]    
        shx = shx_n[:]
        th_mem = th_mem[1:] + [th]
        # update maximal theta
        theta_maxs.append(th[-1]*180/pi)

    theta_max = max(theta_maxs)

    if return_max == 1:
        return theta_max

    if return_max != 1:
        return theta_maxs
