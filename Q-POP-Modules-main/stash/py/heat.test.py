#Transient Heat Equation Test

from fenics import *
import numpy as np

#QPOP Specific Routines
import boundaryUtilities as qpopBU
import solvers as qpopSolve
import tensorUtilities as qpopTU
import terminalPrinting as qpopTP
import meshUtilities as qpopMesh

qpopTP.terminalPrint.printIntroMessage()

T = 2.0
ttime = 0.0
num_steps = 10
dt = T / num_steps
alpha = 3
beta = 1.2

meshInfo = ([0,0],[1,1],[8,8])

mesh = qpopMesh.buildMesh('R', meshInfo)
V = FunctionSpace(mesh, "CG", 1)

u_D = Expression('1 + x[0] * x[0] + alpha*x[1]*x[1] + beta*t',
                 degree=2, alpha=alpha, beta=beta, t=0)

bc = [qpopBU.setupDirichletBoundarys.createDirichlet(V, i, u_D, loc=j) for i in range(2) for j in [0.0, 1.0]]

f = Constant(beta - 2 - 2 * alpha)

u = Function(V)
u0 = interpolate(u_D, V)

heatEqObj = qpopSolve.setup_heatSolve(V, u0, dt, K=Constant(1.0), F=f, bc=bc)

while ttime < T:
    ttime += dt
    u_D.t = ttime
    qpopSolve.solve_Heat(heatEqObj)
    
    u_e = interpolate(u_D, V)
    error = np.abs(u_e.vector()[:] - heatEqObj.u.vector()[:]).max()
    print('T = {}, Error = {}'.format(ttime, error))

qpopTP.terminalPrint.printOutroMessage()