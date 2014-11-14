#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""

boTAB
=====

This solver uses the popular TAB model to simulate the atomization of droplets

Author: Adam O'Brien

"""
from input import readInputFile
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

    print ""
    print "boTAB |"
    print "-------"
    print "        Compute the break-up of a drop in a uniform flow", "\n"
    
    # Open up a configuration file
    
    userInput = readInputFile()
    
    # Set-up the constants in accordance with the input
    
    Cb = userInput["Cb"]
    Ck = userInput["Ck"]
    Cd = userInput["Ck"]
    Cf = userInput["Cf"]
    K = userInput["K"]
    Cv = userInput["Cv"]
    
    # Set-up the freestream and droplet in accordance with the input
                      
    freestream = Freestream(userInput["freestreamRho"],       # density
                            userInput["freestreamMu"],        # viscosity
                            userInput["freestreamVelocity"])  # velocity
        
    droplet = Droplet(userInput["radius"],           # radius
                      userInput["dropletRho"],       # density
                      userInput["dropletMu"],        # viscosity
                      userInput["sigma"],            # surface tension coefficient
                      userInput["dropletPosition"],  # position
                      userInput["dropletVelocity"])  # velocity
                      
    # Set-up the simulation parameters in accordance with the input
                      
    maxTime = userInput["maxTime"]
    nTimeSteps = userInput["nTimeSteps"]
    
    # Initialize misc parameters    
    
    dt = maxTime/nTimeSteps
    t = 0.
    nBreakups = 0
    
    # Open a file
    
    outFile = open("dropBreakup.txt", "w")
    outFile.write("time, droplet_radius, position, velocity\n")
    
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

        # Check the undamped oscillation amplitude. If this condition
        # is not met, break-up is not possible
        
        if Wec + A > 1.:
            
            breakupPossible = True

        else:

            breakupPossible = False
            
        # Update droplet distortion
            
        yn = droplet.y            
            
        droplet.y = computeY(Wec, dt, td, droplet.t, omega, droplet.y, droplet.dydt)
        
        droplet.dydt = computeDyDt(Wec, dt, td, droplet.t, omega, droplet.y, yn, droplet.dydt)
        
        # Break the drop if necessary
        
        if droplet.checkBreakup() and breakupPossible:
            
            # Compute the Sauter Mean Radius (SMR) and use this as the new droplet radius.
            # Since only one droplet needs to be tracked, the others are just ignored
            
            r32 = droplet.radius/(1. + 8.*K*droplet.y**2/20. + \
            droplet.rho*droplet.radius**3*droplet.dydt**2/droplet.sigma*(6.*K - \
            5.)/120.)

            # Compute the normal component of the velocity
            
            vn = droplet.velocity.unitVector().normalVector().scale(Cv*Cb*droplet.radius*droplet.dydt)
            
            droplet = Droplet(r32,             # new radius
                      998.,                    # density
                      8.94e-4,                 # viscosity
                      0.07262,                 # surface tension coefficient
                      droplet.position,        # position
                      droplet.velocity + vn)   # new velocity

            nBreakups += 1
            
            print "\nThe droplet has broken!"
            print "time =", str(t), "s"
            print droplet
            
        t += dt
            
        outFile.write("%s, %s, %s, %s\n"%(str(t), str(droplet.radius), \
        droplet.position, droplet.velocity))
    
    print "\nThe final droplet:"
    print "time =", str(t), "s"
    print "number of breakups:", str(nBreakups)
    print droplet

# Execute the main function   
    
if __name__ == "__main__":
    main()
