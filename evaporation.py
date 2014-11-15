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
    