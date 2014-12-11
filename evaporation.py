# -*- coding: utf-8 -*-

"""

boTAB
=====

This module contains classes and function for the modelling of droplet
evaporation in freestream flows

Author: Adam O'Brien

"""

from math import log, sqrt, fabs, pi

def evaporate(freestream, droplets, dt):

    for droplet in droplets:

        Nu = 2. + droplet.Re(freestream)**0.5*freestream.Pr