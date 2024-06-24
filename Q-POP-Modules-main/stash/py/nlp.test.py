#Nonlinear Poisson Solver Test

from fenics import *

#QPOP Specific Routines
import boundaryUtilities as qpopBU
import solvers as qpopSolve
import tensorUtilities as qpopTU
import terminalPrinting as qpopTP
import meshUtilities as qpopMesh

qpopTP.terminalPrint.printIntroMessage()

meshInfo = ([0,0], [1,1], [32,32])

mesh = qpopMesh.buildMesh('R', meshInfo)
# mesh = UnitSquareMesh(32, 32)

V = FunctionSpace(mesh, "CG", 1)

bc = qpopBU.setupDirichletBoundarys.createDirichlet(V, 0, Constant(1.0), loc=1.0)

f = Expression("x[0] * sin(x[1])", degree=1)
u = Function(V)

nlPoissonObj = qpopSolve.setup_nonlinearPoisson(V, u, NLF=(1 + u**2), bc=bc,
                           F=f, G=Constant(0.0))

qpopSolve.solveNonlinearPoisson(nlPoissonObj)

qpopTP.terminalPrint.printOutroMessage()