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
from TAB import *
from output import *
import copy as cp

def main():

    print ""
    print "boTAB |"
    print "-------"
    print "        Compute the break-up of a drop in a uniform flow", "\n"

    # Open up a configuration file

    userInput = readInputFile()

    # Set-up the freestream

    freestream = Freestream(userInput["freestreamVelocity"],
                            userInput["freestreamGravity"])

    # Set-up droplet initial conditions

    initialDroplet = Droplet(userInput["dropletRadius"],
                             userInput["dropletPosition"],
                             userInput["dropletVelocity"])

    # Set-up the droplet inlet

    dropletInlet = DropletInlet(userInput["inletDropletCreationFrequency"],
                                userInput["inletWidth"],
                                userInput["inletVelocityDeviation"])

    # Set-up the simulation parameters in accordance with the input

    maxTime = userInput["maxTime"]
    nTimeSteps = userInput["nTimeSteps"]

    # Initialize a droplet list, with one copy of the initial droplet

    droplets = [cp.deepcopy(initialDroplet)]

    # Initialize misc parameters

    dt = maxTime/nTimeSteps
    t = [0.]

    # Open a file

    outFile = open("dropBreakup.txt", "w")
    outFile.write("time, droplet_radius, position, velocity\n")

    # Begin the simulation

    print "\nBeginning time-stepping..."

    ###########################################################################
    #                                                                         #
    #                        Main Iteration Loop                              #
    #                                                                         #
    ###########################################################################

    for stepNo in range(1, nTimeSteps + 1):

        for droplet in droplets:

            droplet.advectPredictorCorrector(freestream, dt)

        evaporate(freestream, droplets, dt)
        breakupTab(freestream, droplets, dt)

        dropletInlet.addDrops(initialDroplet, droplets, dt)
        t.append(t[-1] + dt)

    outFile.close()

    print "\nTime-stepping complete. Finalizing output..."

    plotDroplets(droplets)

# Execute the main function

if __name__ == "__main__":
    main()
