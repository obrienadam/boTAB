#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""

boTAB
=====

This solver uses the popular TAB model to simulate the atomization of droplets

Author: Adam O'Brien

"""

from math import exp, cos, sin, sqrt
from fluid import Freestream, Droplet, Vector

# Returns the axis distortion y^n+1

def computeY(Wec, dt, td, t, omega, yn, dydtn):
    
    return Wec + exp(-dt/td)*((yn - Wec)*cos(omega*t) + 1/omega*(dydtn + \
    (yn - Wec)/td)*sin(omega*t))
            
# Returns the time rate of change of axis distortion (dy/dt)^n+1     
            
def computeDyDt(Wec, dt, td, t, omega, yn1, yn, dydtn):
    
    return (Wec - yn1)/td + omega*exp(-dt/td)*(1/omega*(dydtn + (yn - Wec)/td)* \
    cos(omega*t) - (yn - Wec)*sin(omega*t))

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
        
    droplet = Droplet(0.001,           # radius
                      998.,             # density
                      8.94e-4,          # viscosity
                      0.07262,          # surface tension coefficient
                      Vector(0., 0.),   # position
                      Vector(0., 0.))   # velocity
                      
    # Define simulation parameters
                      
    maxTime = 1.
    nTimeSteps = 1800
    dt = maxTime/nTimeSteps
    nBreakups = 0
    
    # Begin the simulation
    
    for stepNo in range(1, nTimeSteps + 1):
        
        # Advect the droplet
        
        droplet.advect(freestream, Cf, dt)
        
        # Compute relevant droplet oscillation parameters
    
        Wec = droplet.weberNoCrit(freestream, Cf, Ck, Cb)
        td = droplet.td(Cd)
        omega = droplet.omega(Cd, Ck)
        
        # Compute undamped oscillation, check if breakup is possible
        
        A = sqrt((droplet.y - Wec)**2 + ((droplet.dydt)/omega)**2)
        
        if Wec + A > 1:
            
            pass
            
        # Update droplet distortion
            
        yn = droplet.y            
            
        droplet.y = computeY(Wec, dt, td, droplet.t, omega, droplet.y, droplet.dydt)
        
        droplet.dydt = computeDyDt(Wec, dt, td, droplet.t, omega, droplet.y, yn, droplet.dydt)
        
        #print droplet.y
        
        # Break the drop if necessary
        
        if droplet.checkBreakup():
            
            r32 = droplet.radius/(1. + 8.*K*droplet.y**2/20. + \
            droplet.rho*droplet.radius**3*droplet.dydt**2/droplet.sigma*(6.*K - \
            5.)/120.)
            
            vn = droplet.velocity.unitVector().unitVector().scale(Cv*Cb*droplet.radius*droplet.dydt)
            
            droplet = Droplet(r32,             # radius
                      998.,                    # density
                      8.94e-4,                 # viscosity
                      0.07262,                 # surface tension coefficient
                      droplet.position,        # position
                      droplet.velocity + vn)   # velocity
            
            print "The droplet has broken!"
            print droplet
    
    print "\nThe final droplet:"        
    print droplet

# Execute the main function   
    
if __name__ == "__main__":
    main()