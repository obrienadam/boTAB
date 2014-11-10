#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""

boTAB
=====

This solver uses the popular TAB model to simulate the atomization of droplets

Author: Adam O'Brien

"""

from math import exp, cos, sin
from fluid import Freestream, Droplet, Vector

# Returns the axis distortion y^n+1

def computeY(Wec, dt, td, t, omega, yn, dydtn):
    
    return Wec + exp(-dt/td)*((yn - Wec)*cos(omega*t) + 1/omega*(dydtn + \
            (yn - Wec)/td)*sin(omega*t))
            
# Returns the time rate of change of axis distortion (dy/dt)^n+1     
            
def computeDyDt(Wec, dt, td, t, omega, yn1, yn, dydtn):
    
    return (Wec - yn1)/td + omega*exp(-dt/td)*(1/omega*(dydtn + (yn - Wec)/td)* \
            cos(omega*dt) - (yn - Wec)*sin(omega*dt))

def main():
    
    # Define the empirical constants
    
    Cb = 0.5
    Ck = 8.
    Cd = 5.
    Cf = 1/3.
    K = 10/3.
    Cv = 1.
    
    freestream = Freestream(1.205,             # density
                            18.27e-6,          # viscosity
                            Vector(200., 0.))  # velocity
        
    droplet = Droplet(0.0005,           # radius
                      998.,             # density
                      8.94e-4,          # viscosity
                      0.07262,          # surface tension coefficient
                      Vector(0., 0.),   # position
                      Vector(0., 0.))   # velocity
                      
    # Define simulation parameters
                      
    maxTime = 0.1
    nTimeSteps = 1000
    dt = maxTime/nTimeSteps
    nBreakups = 0
    
    # Begin the simulation
    
    for stepNo in range(1, nTimeSteps + 1):
        
        # Compute relevant droplet oscillation parameters
    
        We = droplet.weberNo(freestream)
        Wec = droplet.weberNoCrit(freestream, Cf, Ck, Cb)       
        td = droplet.td(Cd)
        omega = droplet.omega(Cd, Ck)
        
        # Compute the droplet distortion       
        
        yn = droplet.y
        droplet.y = computeY(Wec, dt, td, droplet.t, omega, droplet.y, \
        droplet.dydt)
        droplet.dydt = computeDyDt(Wec, dt, td, droplet.t, omega, droplet.dydt, \
        yn, droplet.dydt)
        
        #print droplet.y
        
        # Advect the droplet
        
        
        
        
        
        droplet.advect(freestream, Cb, dt)
        
        # Check for breakup. Create a new droplet in the event of breakup
        
        if droplet.checkBreakup():
            
            # Compute the Sauter Mean Radius (SMR) of the new droplet
            # distribution. For now, only one droplet is tracked            
            """
            r32 = droplet.radius/(1. + 8.*K*droplet.y**2/20. + \
            droplet.rho*droplet.radius**3*droplet.dydt**2/droplet.sigma* \
            (6.*K - 5.)/120.)
            """
            r32 = droplet.radius/2
            # Compute the velocity of the new child droplet. Assume it just
            # travels normal to the path of the parent droplet

            v = droplet.velocity


            # Create new droplet
            
            droplet = Droplet(r32,   # radius
                      998.,             # density
                      8.94e-4,          # viscosity
                      0.07262,          # surface tension coefficient
                      droplet.position,   # position
                      v)   # velocity
                      
            nBreakups += 1
                      
            print droplet
        
    
# Execute the main function    
    
if __name__ == "__main__":
    main()