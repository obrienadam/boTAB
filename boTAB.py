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
from output import *
import copy as cp

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
                            userInput["freestreamVelocity"],  # velocity
                            userInput["freestreamGravity"])   # gravity
    # Set-up droplet initial conditions

    initialDroplet = Droplet(userInput["radius"],    # radius
                      userInput["dropletRho"],       # density
                      userInput["dropletMu"],        # viscosity
                      userInput["sigma"],            # surface tension coefficient
                      userInput["boilingTemp"],      # boiling temperature
                      userInput["latentHeat"],       # latent heat of evaporation
                      userInput["specificHeat"],     # specific heat
                      userInput["dropletK"],         # thermal conductivity
                      userInput["dropletPosition"],  # position
                      userInput["dropletVelocity"])  # velocity

    # Set-up the droplet inlet

    dropletInlet = DropletInlet(userInput["newDropletFrequency"],
                                userInput["inletWidth"],
                                userInput["velocityDeviation"])

    # Set-up the simulation parameters in accordance with the input

    maxTime = userInput["maxTime"]
    nTimeSteps = userInput["nTimeSteps"]

    # Initialize a droplet list, with one copy of the initial droplet

    droplets = [cp.deepcopy(initialDroplet)]

    # Initialize misc parameters

    dt = maxTime/nTimeSteps
    t = [0.]
    smr = [initialDroplet.radius]
    nBreakups = 0

    # Open a file

    outFile = open("dropBreakup.txt", "w")
    outFile.write("time, droplet_radius, position, velocity\n")

    # Begin the simulation

    print "\nBeginning time-stepping..."

    for stepNo in range(1, nTimeSteps + 1):

        for droplet in droplets:

            droplet.advectPredictorCorrector(freestream, dt)

        dropletInlet.addDrops(initialDroplet, droplets, dt)

    outFile.close()

    plotDroplets(droplets)

# Execute the main function

if __name__ == "__main__":
    main()
