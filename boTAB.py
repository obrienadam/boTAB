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
    print "        Compute the break-up of a drop in a uniform cross-flow", "\n"

    # Open up a configuration file

    userInput = readInputFile()

    freestream = Freestream()
    initialDroplet = Droplet()
    dropletInlet = DropletInlet()

    # Set object parameters from the input file

    setObjectParametersFromInput(userInput, freestream, initialDroplet, dropletInlet)

    # Set-up the simulation parameters in accordance with the input

    maxTime = userInput["maxTime"]
    nTimeSteps = userInput["nTimeSteps"]

    # Initialize a droplet list, with one copy of the initial droplet

    droplets = [cp.deepcopy(initialDroplet)]

    # Initialize misc parameters

    dt = maxTime/nTimeSteps
    t = [0.]
    nChildDroplets = 0

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
        nChildDroplets += breakupTab(freestream, droplets, dt)

        dropletInlet.addDrops(initialDroplet, droplets, dt)
        t.append(t[-1] + dt)

        if stepNo%(nTimeSteps/20) == 0:

            completionPercentage = float(stepNo)/float(nTimeSteps)*100.

            print "-----------------------------------------------------------"
            print "Time-stepping completion     : %s%%"%(completionPercentage)
            print "Number of droplets in domain :", len(droplets)
            print "Simulation time elapsed      : %s seconds"%(t[-1])
            print "Simulation time remaining    : %s seconds"%(maxTime - t[-1])
            print "Number of child drops        :", nChildDroplets

    print "\nTime-stepping complete. Finalizing output..."

    plotDroplets(droplets)

# Execute the main function

if __name__ == "__main__":
    main()
