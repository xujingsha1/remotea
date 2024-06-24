#Linear Poisson Equation Test

from fenics import *

#QPOP Specific Routines
import boundaryUtilities as qpopBU
import solvers as qpopSolve
import tensorUtilities as qpopTU
import terminalPrinting as qpopTP
import meshUtilities as qpopMesh

qpopTP.terminalPrint.printIntroMessage()

meshInfo = ([0,0],[1,1],[32,32])

mesh = qpopMesh.buildMesh('R', meshInfo)
# mesh = UnitSquareMesh(32, 32)

V = FunctionSpace(mesh, "CG", 1)

bc = [qpopBU.setupDirichletBoundarys.createDirichlet(V, 0, Constant(0.0), loc=1.0), qpopBU.setupDirichletBoundarys.createDirichlet(V, 0, Constant(0.0), loc=0.0)]

f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=1)
g = Expression("sin(5*x[0])", degree=1)
u = Function(V)

lPoissonObj = qpopSolve.setup_linearPoissonSolve(V, u, bc=bc, F=f, G=g)

qpopSolve.solve_LinearPoisson(lPoissonObj)

qpopTP.terminalPrint.printOutroMessage()