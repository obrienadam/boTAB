# -*- coding: utf-8 -*-

"""

boTAB
=====

This module contains the functions necessary for TAB breakup

Author: Adam O'Brien

"""

from math import exp, sin, cos, fabs
import numpy as np
import copy as cp
from fluid import *

# Model constants

Cb = 0.5
Ck = 8.
Cd = 5.
Cf = 1./3.
Cd = 5.
K = 10./3.
Cv = 1.

def getSMR(droplet):

    # This function computes the Sauter Mean Radius (SMR) of the child droplets after a break-up

    rOverRmean = 1. + (K/5.)*Ck*Cb**2 + droplet.rho*droplet.radius**3/droplet.sigma*Cb**2*droplet.dydt**2*(6.*K - 5.)/30.

    return droplet.radius/rOverRmean

def breakupTab(freestream, droplets, dt):

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

        if droplet.y >= 1.:

            rSmr = getSMR(droplet)

            # Begin sampling droplets

            volOfNewDrops = 0.

            while True:
                
                newRadius = np.random.normal(rSmr, 0.2*rSmr)
                
                if newRadius < 0.:
                    
                    continue
                
                randomNo = np.random.uniform(-1, 1)                
                
                newVelocity = droplet.velocity + droplet.velocity.normalVector().unitVector()*Cb*droplet.radius*droplet.dydt*(randomNo/fabs(randomNo))                
                
                newDroplet = Droplet(newRadius, cp.copy(droplet.position), cp.copy(newVelocity))
                
                newDroplets.append(newDroplet)
                
                volOfNewDrops += newDroplets[-1].volume()
                
                if volOfNewDrops >= droplet.volume():
                    
                    break
                
    # Remove the old parent droplets
                
    droplets[:] = [droplet for droplet in droplets if droplet.y < 1.]
    
    # Add the newly created child droplets
    
    for droplet in newDroplets:
        
        droplets.append(cp.deepcopy(droplet))
        
    # Return the number of droplets created
        
    return len(newDroplets)
