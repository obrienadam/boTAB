# -*- coding: utf-8 -*-

"""

boTAB
=====

This module contains classes and function for the modelling of droplet
evaporation in freestream flows

Author: Adam O'Brien

"""

from math import log, sqrt, fabs

def evaporate(freestream, droplet, dt):

    B = droplet.specificHeat*(freestream.temperature - droplet.boilingTemp)/\
    droplet.latentHeat

    beta = 8.*freestream.k/(droplet.rho*droplet.specificHeat)*log(1. + B)

    D02 = (2.*droplet.radius)**2

    D2 = D02 + beta*dt

    if D2 <= 1e-14:

        droplet.radius = 1e-7;

    else:

        droplet.radius = sqrt(D2)*0.5