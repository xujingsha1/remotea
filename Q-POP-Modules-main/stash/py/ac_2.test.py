# Allen-Cahn Equation Test
# In This test, we demonstrate that you can define your driving force in a seperate file
# Then you can import into your code. This therefore eliminates the need for FEniCS to perform
# any type of symbolic differentiation for driving forces.

from fenics import *
import numpy as np
import matplotlib.pyplot as plt

#QPOP Specific Routines
import boundaryUtilities as qpopBU
import solvers as qpopSolve
import tensorUtilities as qpopTU
import terminalPrinting as qpopTP
import meshUtilities as qpopMesh

#Import the energy function/driving force
from ac_ener import drivingForce, init_condition

# Start Up QPOP and Print Message
qpopTP.terminalPrint.printIntroMessage()

# Timestep and iteration information
dt = 5.0e-3
ttime = 0.0
n = 0
T = 1.0

# Information regarding the mesh
# ==============================================================================
meshInfo = ([-1], [1], 1000)

# Create Periodic Boundary Conditions
pbc = qpopBU.setupPeriodicBoundarys.periodic_1D(0, O=meshInfo[0][0], L=meshInfo[1][0]-meshInfo[0][0])

#Build the mesh using qpop interface utilities
mesh = qpopMesh.buildMesh('R', meshInfo)
#===============================================================================

#Setup the Function Space, the solution functions, and initial Conditions
# ==============================================================================

# Create Function Space on our mesh with our periodic condition
V = FunctionSpace(mesh, "CG", 1, constrained_domain=pbc)

# Create solution functions (u is the timestep we are solving for, u0 is the previous time solution)
u, u0 = Function(V), Function(V)

# Create the initial condition and apply to our solution function
u_init = Expression(init_condition(u), degree=2)
u.interpolate(u_init)
u0.interpolate(u_init)

# ==============================================================================

# Setup the allen-cahn Equation
# ==============================================================================
gamma1 = 0.0001
gamma2 = 4

acEqObj = qpopSolve.setupScalarAllenCahnSolver(V, u0, u, drivingForce, gamma1, dt, bc=None, args=(gamma1, gamma2))
# ==============================================================================

# Information for easy in-loop visualization
# ==============================================================================
dof_coords = V.tabulate_dof_coordinates()
u_dofs = V.dofmap().dofs()
dofs = np.squeeze(dof_coords[u_dofs])
asc_order = np.argsort(dofs)
# ==============================================================================

# Main Timestepping loop
# ==============================================================================
while ttime < T:

    # Solve Allen-Cahn
    # ==========================================================================
    qpopSolve.solveScalarAllenCahn(acEqObj)

    # Plot order parameters
    # ==========================================================================
    if n % round((T/dt)/6,0) == 0:
        plt.plot(dofs[asc_order], acEqObj.u.vector()[asc_order], label='t={}'.format(round(ttime,2)))

    # Increment Time
    # ==========================================================================
    ttime += dt
    n += 1

# ==============================================================================

# Finish Plotting and Save Figure
# ==============================================================================
plt.legend()
plt.savefig('AC_2_Test.png')
plt.close()
# ==============================================================================

# Exit QPOP and Print Message
# ==============================================================================
qpopTP.terminalPrint.printOutroMessage()
