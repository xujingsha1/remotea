import time, json, itertools
from fenics import *
import numpy as np
# import mpi4py
from ufl import nabla_div
from mshr import *

class generalRoutines():

    def mapLoop(x,y,dims):
        for d in range(len(x)):
            if d not in dims:
                y[d] = x[d]
        return y

    def mapArbitaryLoop(x,y):
        for d in range(len(x)):
            y[d] = -1000
        return y

class setupPeriodicBoundarys():

    def periodic_1D(dim, O=0, L=1):

        class PeriodicBC_1D(SubDomain):
            def inside(self, x, on_boundary):
                return bool((near(x[dim],O)) and (not near(x[dim],L)) and on_boundary)

            def map(self,x,y):
                if near(x[dim],L):
                    y[dim] = x[dim] - L
                    y = generalRoutines.mapLoop(x,y,[dim])
                else:
                    y = generalRoutines.mapArbitaryLoop(x,y)

        return PeriodicBC_1D()

    def periodic_2D(dims, O=[0,0], L=[1,1]):

        class PeriodicBC_2D(SubDomain):
            def inside(self, x, on_boundary):
                return bool((near(x[dims[0]],O[0]) or near(x[dims[1]],O[1])) and
                (not ((near(x[dims[0]], O[0]) and near(x[dims[1]], L[1])) or
                (near(x[dims[0]], L[0]) and near(x[dims[1]], O[1])))) and on_boundary)

            def map(self,x,y):
                if near(x[dims[0]], L[0]) and near(x[dims[1]], L[1]):
                    for dim in dims:
                        y[dim] = x[dim] - L[dim]
                        y = generalRoutines.mapLoop(x,y,dims)
                elif near(x[dims[0]], L[0]):
                    y[dims[0]] = x[dims[0]] - L[0]
                    y = generalRoutines.mapLoop(x,y,[dims[0]])
                elif near(x[dims[1]], L[1]):
                    y[dims[1]] = x[dims[1]] - L[1]
                    y = generalRoutines.mapLoop(x,y,[dims[1]])
                else:
                    y = generalRoutines.mapArbitaryLoop(x,y)
        return PeriodicBC_2D()

class setupDirichletBoundarys():

    def createDirichlet(mesh, dim, value, loc=0):

        def dirBoundary(x, on_boundary):
            return on_boundary and near(x[dim], loc)

        return DirichletBC(mesh, value, dirBoundary)
