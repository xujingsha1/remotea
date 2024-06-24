import time, json, itertools
from fenics import *
import numpy as np
# import mpi4py
from ufl import nabla_div
from mshr import *
from terminalPrinting import terminalPrint as tP
# import terminalPrinting.terminalPrint as tP

def buildMesh(meshType, meshInfo):
    
    if meshType[0].lower() == 'c':
        print('Need to implement')
    elif meshType[0].lower() == 'r':
        if len(meshInfo[0]) == 1:
            tP.printNamedLine('Building Interval Mesh')
            return buildIntervalMesh(meshInfo)
        elif len(meshInfo[1]) == 2:
            tP.printNamedLine('Building Rectangle Mesh')
            return buildRectangleMesh(meshInfo)
        elif len(meshInfo[2]) == 3:
            tP.printNamedLine('Building Box Mesh')
            return buildBoxMesh(meshInfo)
            
def buildBoxMesh(meshInfo):
    
    return BoxMesh(Point(meshInfo[0]), Point(meshInfo[1]), *meshInfo[2])

def buildRectangleMesh(meshInfo):
    
    return RectangleMesh(Point(meshInfo[0]), Point(meshInfo[1]), *meshInfo[2])

def buildIntervalMesh(meshInfo):
    
    return IntervalMesh(meshInfo[2], meshInfo[0][0], meshInfo[1][0])