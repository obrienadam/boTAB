#/usr/bin/env python2
# -*- coding: utf-8 -*-

"""

boTAB
=====

This module contains classes and function for the modelling of droplets in
freestream flows

Author: Adam O'Brien

"""

# Vector class for position/velocity

class Vector(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
    def __add__(self, other):
        return Vector(self.x + other.x, self.y + self.y)
            
    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y)
        
    def rVector(self, other):
        return other - self
            
def dot(u, v):
    return u.x*v.x + u.y*v.y
    
# Freestream class for representing the freestream flow

class Freestream(object):
    def __init__(self, rho, mu, velocity):
        self.rho = rho
        self.mu = mu
        self.velocity = velocity
    
# Droplet class for representing droplets    
    
class Droplet(object):
    def __init__(self, radius, rho, mu, sigma, position, velocity):
        self.radius = radius
        self.rho = rho
        self.mu = mu
        self.sigma = sigma
        self.position = position
        self.velocity = velocity
        
    def weberNo(self, freestream):
        vRel = self.velocity.rVector(freestream.velocity)
        
        return freestream.rho*dot(vRel, vRel)*self.radius/self.sigma
        
    def weberNoCrit(self, freestream, Cf, Ck, Cb):
        return Cf/(Ck*Cb)*self.weberNo(freestream)
        