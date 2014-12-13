# -*- coding: utf-8 -*-

"""

boTAB
=====

This module is for various plotting outputs

Author: Adam O'Brien

"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation

def plotDroplets(droplets):

    xcoords = [droplet.position.x for droplet in droplets]
    ycoords = [droplet.position.y for droplet in droplets]
    radii = [50000.*droplet.radius for droplet in droplets]

    plt.axis('equal')
    plt.grid(True)
    plt.title('Droplet Distribution', fontsize=20)
    plt.xlabel('x (m)', fontsize=16)
    plt.ylabel('y (m)', fontsize=16)

    plt.scatter(xcoords, ycoords, s=radii, alpha=0.5)

    plt.show()