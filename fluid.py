#/usr/bin/env python2
# -*- coding: utf-8 -*-

"""

boTAB
=====

This module contains classes and function for the modelling of droplets in
freestream flows

Author: Adam O'Brien

"""

from math import sqrt, pi
import copy as cp
import random

# Vector class for position/velocity

class Vector(object):

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):

        return "%s, %s"%(self.x, self.y)

    def __add__(self, other):

        return Vector(self.x + other.x, self.y + other.y)

    def __sub__(self, other):

        return Vector(self.x - other.x, self.y - other.y)

    def __mul__(self, other):

        return Vector(self.x*other, self.y*other)

    def mag(self):
        return sqrt(self.x**2 + self.y**2)

    def normalVector(self):
        return Vector(self.y, -self.x)

    def scale(self, other):
        return Vector(self.x*other, self.y*other)

    def rVector(self, other):

        return other - self

    def unitVector(self):

        return self.scale(1./self.mag())

def dot(u, v):

    return u.x*v.x + u.y*v.y

# Freestream class for representing the freestream flow

class Freestream(object):

    def __init__(self, velocity, gravity):

        # constructed properties

        self.velocity = velocity
        self.gravity = gravity
        
        # default air proplerties (can be changed)

        self.rho = 1.205
        self.mu = 18.27e-6
        self.Tambient = 298.
        self.Cp = 1.0005
        self.k = 0.0257
        self.Pr = self.Cp*self.mu/self.k

# Droplet class for representing droplets

class Droplet(object):

    def __init__(self, radius, position, velocity):
        
        # constructed properties

        self.radius = radius
        self.position = position
        self.velocity = velocity

        # default water properties (can be changed)

        self.rho = 998.
        self.mu = 8.94e-4        
        self.sigma = 0.07262
        self.Tboil = 373.
        self.L = 2257.
        self.Cp = 4.183
        self.k = 0.58
        
        # computed quantities

        self.diameter = 2.*self.radius        
        self.volume = (4./3.)*pi*self.radius**3
        self.area = pi*self.radius**2
        self.mass = self.rho*self.volume

        # TAB properties

        self.y = 0.
        self.dydt = 0.

    def __repr__(self):

        return "Radius: %s\nPosition: %s\nVelocity: %s"%(str(self.radius), \
        self.position, self.velocity)

    def weberNo(self, freestream):

        vRel = self.velocity.rVector(freestream.velocity)

        return freestream.rho*dot(vRel, vRel)*self.radius/self.sigma

    def weberNoCrit(self, freestream, Cf, Ck, Cb):

        return Cf/(Ck*Cb)*self.weberNo(freestream)

    def td(self, Cd):

        return 2.*self.rho*self.radius ** 2/(Cd*self.mu)

    def omega(self, Cd, Ck):

        td = self.td(Cd)

        return sqrt(Ck*self.sigma/(self.rho*self.radius ** 3) - 1./td**2)

    def Re(self, freestream):

        return freestream.rho*(freestream.velocity - self.velocity).mag()*2.*self.radius/freestream.mu

    def dragCoefficient(self, freestream):

        # This drag coefficient for a sphere is based on the correlation of
        # F.A. Morrison in "An Introduction to Fluid Mechanics"

        Re = self.Re(freestream)

        return 24./Re + 2.6*(Re/5.)/(1. + (Re/5.)**1.52) + 0.411*(Re/263000.)**-7.94/(1. + (Re/263000.)**-8.) + (Re**0.8/461000.)

    def dragForce(self, freestream):

        Cd = self.dragCoefficient(freestream)

        vRel = freestream.velocity - self.velocity

        return vRel.unitVector()*0.5*freestream.rho*dot(vRel, vRel)*Cd*self.area

    def acceleration(self, freestream):

        return self.dragForce(freestream).scale(1./self.mass) + freestream.gravity

    def advectEuler(self, freestream, dt):

        a = self.acceleration(freestream)

        self.position += self.velocity*dt + a*0.5*dt**2
        self.velocity += a*dt

    def advectPredictorCorrector(self, freestream, dt):

        originalPosition = self.position

        a = self.acceleration(freestream)
        f1 = self.velocity + a*0.5*dt

        self.position += f1*dt

        a = self.acceleration(freestream)
        f2 = self.velocity + a*0.5*dt

        self.position = originalPosition + (f1 + f2)*0.5*dt
        self.velocity += a*dt

    def checkBreakup(self):

        if self.y > 1.:
            return True
        else:
            return False

    def printAll(self):

        print "Radius:", self.radius
        print "Rho:", self.rho
        print "mu:", self.mu
        print "sigma:", self.sigma
        print "Boiling Temp:", self.boilingTemp
        print "Latent Heat:", self.latentHeat
        print "Specific Heat:", self.specificHeat
        print "k:", self.k
        print "Position:", self.position
        print "Velocity:", self.velocity

class DropletInlet(object):

    def __init__(self, newDropletFrequency, inletWidth, velocityDeviation):

        self.newDropletFrequency = newDropletFrequency
        self.inletWidth = inletWidth
        self.velocityDeviation = velocityDeviation
        self.timeSinceLastDroplet = 0.
        self.newDropletPeriod = 1./self.newDropletFrequency
        self.dropsAdded = 0

    def addDrops(self, initialDroplet, droplets, dt):

        self.timeSinceLastDroplet += dt

        nDropsToAdd = int(self.timeSinceLastDroplet/self.newDropletPeriod)
        self.timeSinceLastDroplet -= float(nDropsToAdd)*self.newDropletPeriod

        if nDropsToAdd > 0:
            self.dropsAdded += nDropsToAdd
            print "Drops added:", self.dropsAdded

        for i in range(0, nDropsToAdd):

            droplets.append(cp.deepcopy(initialDroplet))
            droplets[-1].position.x += random.uniform(-0.5*self.inletWidth, 0.5*self.inletWidth)
            droplets[-1].velocity.x += random.normalvariate(0., self.velocityDeviation.x)
            droplets[-1].velocity.y += random.normalvariate(0., self.velocityDeviation.y)



