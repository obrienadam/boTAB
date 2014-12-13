# -*- coding: utf-8 -*-

"""

boTAB
=====

This module contains classes and function for the modelling of droplet
evaporation in freestream flows

Author: Adam O'Brien

"""

from math import log, sqrt, fabs, pi

def clausiusClapeyron(freestream, droplet):

    # This function determines the water vapour mass fraction at the surface
    # of a droplet

    return 1./(1. + freestream.Pambient*freestream.M/((droplet.Pvap(freestream) - 1.)*droplet.M))

def evaporate(freestream, droplets, dt):

    for i in range(0, len(droplets)):

        Yls = clausiusClapeyron(freestream, droplets[i])

        BM = Yls/(1. - Yls)

        Re = droplets[i].Re(freestream)
        Pr = freestream.Pr

        gamma = 8.*freestream.k*log(1. + BM)/(freestream.Cp*droplets[i].rho)*(1. + 0.3*Re**0.5*Pr**(1./3.))

        D2 = droplets[i].diameter()**2 - gamma*dt

        if D2 > 0.:

            droplets[i].radius = 0.5*D2**0.5

        else:

            droplets[i].radius = 0.

    # Discard any droplets that have a radius less than the tolerance, ie they
    # are completely evaporated

    droplets[:] = [droplet for droplet in droplets if not droplet.radius <= 1e-10]