# Elasticity Solver Test

from fenics import *
import numpy as np
from ufl import nabla_div

import warnings
warnings.filterwarnings("ignore")

#QPOP Specific Routines
import boundaryUtilities as qpopBU
import solvers as qpopSolve
import tensorUtilities as qpopTU
import terminalPrinting as qpopTP
import meshUtilities as qpopMesh

qpopTP.terminalPrint.printIntroMessage()

L, W = 1, 0.2
mu = 1
rho = 1
delta = W / L
gamma = 0.4 * delta**2
beta = 1.25
meshInfo = ([0,0,0], [L,W,W], [10, 3, 3])

mesh = qpopMesh.buildMesh('R', meshInfo)

V = VectorFunctionSpace(mesh, 'P', 1)

bc = qpopBU.setupDirichletBoundarys.createDirichlet(V, 0, Constant((0,0,0)), loc=0)

f = Constant((0, 0, -rho * gamma))
T = Constant((0, 0, 0))

elastEqObj = qpopSolve.setup_Elasticity(V, F=f, T=T, bc=bc, C=None, lame=(beta, mu))

qpopSolve.solve_Elasticity(elastEqObj)

qpopTP.terminalPrint.printOutroMessage()

