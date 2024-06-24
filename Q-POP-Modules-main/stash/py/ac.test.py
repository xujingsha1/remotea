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

dt = 5.0e-3
meshInfo = ([-1], [1], 1000)

pbc = qpopBU.setupPeriodicBoundarys.periodic_1D(0, O=meshInfo[0][0], L=meshInfo[1][0]-meshInfo[0][0])

mesh = qpopMesh.buildMesh('R', meshInfo)
V = FunctionSpace(mesh, "CG", 1, constrained_domain=pbc)

u_init = Expression("pow(x[0],2) * sin(2*pi*x[0])", degree=2)

gamma1 = 0.0001
gamma2 = 4

u, u0 = Function(V), Function(V)
u.interpolate(u_init)
u0.interpolate(u_init)

drivingForce = gamma2/gamma1 * (u**3 - u)

acEqObj = qpopSolve.setupScalarAllenCahnSolver(V, u0, u, drivingForce, gamma1, dt, bc=None)

ttime = 0.0
n = 0

dof_coords = V.tabulate_dof_coordinates()
u_dofs = V.dofmap().dofs()
dofs = np.squeeze(dof_coords[u_dofs])
asc_order = np.argsort(dofs)

while ttime < 1.0:
    qpopSolve.solveScalarAllenCahn(acEqObj)
    
    if n % round((1.0/dt)/6,0) == 0:
        plt.plot(dofs[asc_order], acEqObj.u.vector()[asc_order], label='t={}'.format(ttime))
    
    ttime += dt
    n += 1
    
plt.legend()
plt.savefig('AC_Test.png')
plt.close()

qpopTP.terminalPrint.printOutroMessage()