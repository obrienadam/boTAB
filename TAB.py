# -*- coding: utf-8 -*-

"""

boTAB
=====

This module contains the functions necessary for TAB breakup

Author: Adam O'Brien

"""

from math import exp, sin, cos

def breakupTab(freestream, droplets, dt):

    Cb = 0.5
    Ck = 8.
    Cd = 5.
    Cf = 1./3.
    K = 10./3.
    Cv = 1.

    newDroplets = []

    for droplet in droplets:

        Wec = droplet.We(freestream)*Cf/(Ck*Cb)
        td = droplet.rho*droplet.diameter()**2/(2.*Cd*droplet.mu)
        omega = (8.*Ck*droplet.sigma/(droplet.rho*droplet.diameter()**3) - (1./td**2))**0.5

        yn = droplet.y

        droplet.y = Wec + exp(-dt/td)* \
                    ((yn - Wec)*cos(omega*dt) + 1./omega*(droplet.dydt + (yn - Wec)/td)*sin(omega*dt))

        droplet.dydt = (Wec - droplet.y)/td + \
                       omega*exp(-dt/td)* \
                       (1./omega*(droplet.dydt + (yn - Wec)/td)*cos(omega*dt) - (yn - Wec)*sin(omega*dt))

        if droplet.y > 1.:

            print "A droplet has broken!"


