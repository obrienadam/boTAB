# -*- coding: utf-8 -*-

"""

boTAB
=====

This module is for various plotting outputs

Author: Adam O'Brien

"""

import matplotlib.pyplot as plt

def plotDroplets(droplets):

    xcoords = []
    ycoords = []

    for droplet in droplets:

        xcoords.append(droplet.position.x)
        ycoords.append(droplet.position.y)

    plt.plot(xcoords, ycoords, 'o')
    plt.show()