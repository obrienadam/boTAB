#!/bin/usr/env python2
# -*- coding: utf-8 -*-

"""

boTAB
=====

This solver uses the popular TAB model to simulate the atomization of sprays.

Author: Adam O'Brien

"""

from math import exp, cos, sin
from fluid import Freestream, Droplet, Vector

def computeY(Wec, dt, td, t, omega, yn, dydtn):
    return Wec + exp(-dt/td)*((yn - Wec)*cos(omega*t) + 1/omega*(dydtn + \
            (yn - Wec)/td)*sin(omega*t))
            
def computdYdt(Wec, dt, td, t, omega, yn, dydtn):
    pass

def main():
    
    freestream = Freestream(1.205,             # density
                            18.27e-6,          # viscosity
                            Vector(200., 0.))  # velocity
        
    droplet = Droplet(0.0005,           # radius
                      998.,             # density
                      8.94e-4,          # viscosity
                      0.07262,          # surface tension coefficient
                      Vector(0., 0.),   # position
                      Vector(0., 0.))   # velocity
    
# Execute the main function    
    
if __name__ == "__main__":
    main()