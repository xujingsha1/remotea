# Allen-Cahn Equation Test

from fenics import *
import numpy as np
import matplotlib.pyplot as plt

#QPOP Specific Routines
import boundaryUtilities as qpopBU
import solvers as qpopSolve
import tensorUtilities as qpopTU
import terminalPrinting as qpopTP
import meshUtilities as qpopMesh

qpopTP.terminalPrint.printIntroMessage()

dt = 5.0e-6
meshInfo = ([0,0], [1,1], [96,96])
eps = 1.0e-2
theta = 0.5

mesh = qpopMesh.buildMesh('R', meshInfo)

P1 = FiniteElement('CG', mesh.ufl_cell(), 1)
ME = FunctionSpace(mesh, P1*P1)

u, u0 = Function(ME), Function(ME)
c, mu = split(u)

u.vector()[:] = 0.63 + 0.05 * (0.5 - np.random.random(len(u.vector()[:])))
u0.vector()[:] = 0.63 + 0.05 * (0.5 - np.random.random(len(u.vector()[:])))

c = variable(c)
f = 100 * c**2 * (1-c)**2
dfdc = diff(f, c)

chEqObj = qpopSolve.setup_CahnHilliardSolver(ME, u, u0, dt, eps, dfdc, theta=0.5)

ttime = 0
while ttime <= 50 * dt:
    qpopSolve.solve_CahnHilliard(chEqObj)
    ttime += dt

qpopTP.terminalPrint.printOutroMessage()