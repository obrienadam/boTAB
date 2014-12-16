# -*- coding: utf-8 -*-

"""

boTAB
=====

This module is for simple file input and configuration

Author: Adam O'Brien

"""

from fluid import Vector

def str2num(string):

    if string.partition("/")[1] == "/":

        string = string.partition("/")

        return float(string[0])/float(string[2])

    elif string.partition(",")[1] == ",":

        string = string.partition(",")

        return Vector(float(string[0]), float(string[2]))

    else:

        try:

            return int(string)

        except ValueError:

            return float(string)


def process(inputDict, line):

    line = line.replace("\n", "")
    line = line.replace(" ", "")
    line = line.partition("#")
    line = line[0]

    if line == "":
        return

    line = line.partition("=")

    if not line[1] == "=":
        print "Warnining, potentially bad input:", line

    inputDict[line[0]] = str2num(line[2])

    return

def readInputFile(filename = "config.in"):

    userInput = {}

    print "Reading from input file \"config.in\"..."

    with open(filename) as inFile:

        for line in inFile:

            process(userInput, line)

    print "The following input parameters have been loaded from \"{}\"".format(filename), "\n"

    for parameter in userInput:

        print parameter, "=", userInput[parameter]

    print "Finished reading from input file \"config.in\"."

    return userInput

def setObjectParametersFromInput(userInput, freestream, droplet, inlet):

    # Freestrean properties

    freestream.rho = userInput["freestreamRho"]
    freestream.mu = userInput["freestreamMu"]
    freestream.Tambient = userInput["freestreamTambient"]
    freestream.Cp = userInput["freestreamCp"]
    freestream.K = userInput["freestreamK"]
    freestream.velocity = userInput["freestreamVelocity"]
    freestream.gravity = userInput["freestreamGravity"]
    freestream.M = userInput["freestreamM"]

    # Initial droplet properties

    droplet.radius = userInput["dropletRadius"]
    droplet.rho = userInput["dropletRho"]
    droplet.mu = userInput["dropletMu"]
    droplet.sigma = userInput["dropletSigma"]
    droplet.Tboil = userInput["dropletTboil"]
    droplet.L = userInput["dropletL"]
    droplet.Cp = userInput["dropletCp"]
    droplet.K = userInput["dropletK"]
    droplet.position = userInput["dropletPosition"]
    droplet.velocity = userInput["dropletVelocity"]
    droplet.M = userInput["dropletM"]

    # Droplet inlet properties

    inlet.newDropletFrequency = userInput["inletDropletCreationFrequency"]
    inlet.inletWidth = userInput["inletWidth"]
    inlet.velocityDeviation = userInput["inletVelocityDeviation"]
