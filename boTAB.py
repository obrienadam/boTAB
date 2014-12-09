#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""

boTAB
=====

This solver uses the popular TAB model to simulate the atomization of droplets

Author: Adam O'Brien

"""
from input import *
from math import exp, cos, sin, sqrt
from fluid import *
from evaporation import *
from matplotlib.pyplot import *

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

    # Set-up the freestream

    freestream = Freestream(userInput["freestreamRho"],       # density
                            userInput["freestreamMu"],        # viscosity
                            userInput["temperature"],         # temperature
                            userInput["freestreamK"],         # thermal conductivity
                            userInput["freestreamVelocity"])  # velocity
                            
    # Set-up droplet initial conditions

    droplet = Droplet(userInput["radius"],           # radius
                      userInput["dropletRho"],       # density
                      userInput["dropletMu"],        # viscosity
                      userInput["sigma"],            # surface tension coefficient
                      userInput["boilingTemp"],      # boiling temperature
                      userInput["latentHeat"],       # latent heat of evaporation
                      userInput["specificHeat"],     # specific heat
                      userInput["dropletK"],         # thermal conductivity
                      userInput["dropletPosition"],  # position
                      userInput["dropletVelocity"])  # velocity

    # Set-up the simulation parameters in accordance with the input

    maxTime = userInput["maxTime"]
    nTimeSteps = userInput["nTimeSteps"]
    nDroplets = userInput["nDroplets"]
    
    # Initialize a droplet list
    
    droplets = [droplet]*nDroplets

    # Initialize misc parameters

    dt = maxTime/nTimeSteps
    t = [0.]
    smr = [droplet.radius]
    nBreakups = 0

    # Open a file

    outFile = open("dropBreakup.txt", "w")
    outFile.write("time, droplet_radius, position, velocity\n")

    # Begin the simulation

    print "\nBeginning time-stepping..."

    for stepNo in range(1, nTimeSteps + 1):

        pass

    outFile.close()

    print "\nTime-stepping complete."

    print "\nThe final droplet:"
    print "time =", str(t[-1]), "s"
    print "number of breakups:", str(nBreakups)
    print droplet

    plot(t, smr, linewidth=2.0)
    axis([0., t[-1], 0., 1.1*userInput["radius"]])
    title('Droplet Sauter Mean Radius History')
    xlabel('Time (s)')
    ylabel('Sauter Mean Radius (SMR) (m)')
    grid(True)
    show()

# Execute the main function

if __name__ == "__main__":
    main()
