# -*- coding: utf-8 -*-

"""
This module is for simple file input and configuration

Author: Adam O'Brien

"""

from fluid import Vector

def str2num(string):
    
    if string.partition("/")[1] == "/":
        
        string = string.partition("/")
        
        return float(string[0])/float(string[2])
        
    elif string.partition(",")[1] == ",":
        
        string = string.partition(",")
        
        return Vector(float(string[0]), float(string[2]))
        
    else:
        
        try:
            
            return int(string)
            
        except ValueError:
            
            return float(string)
        

def process(inputDict, line):

    line = line.replace("\n", "")
    line = line.replace(" ", "")
    line = line.partition("#")
    line = line[0]
    
    if line == "":
        return
    
    line = line.partition("=")
    
    if not line[1] == "=":
        print "Warnining, potentially bad input:", line
    
    inputDict[line[0]] = str2num(line[2])
    
    return

def readInputFile(filename = "config.in"):
    
    userInput = {}
    
    with open(filename) as inFile:
        
        for line in inFile:
            
            process(userInput, line)
            
    print "The following input parameters have been loaded from \"{}\"".format(filename), "\n"    
    
    for parameter in userInput:
        
        print parameter, "=", userInput[parameter]
            
    return userInput
            