import time, json, itertools
from fenics import *
import numpy as np
# import mpi4py
from ufl import nabla_div
from mshr import *
from terminalPrinting import terminalPrint as tP
from tensorUtilities import stressStrainCalculators
import inspect

# comm = MPI.comm_world

class solverObject():
    def __init__(self, name):
        self.name = name
        
def setup_Elasticity(V, F=Constant(0.0), T=Constant(0.0), bc=None, C=None, lame=(1.25, 1)):
    
    v = TestFunction(V)
    u = TrialFunction(V)
    
#     F = inner(stressStrainCalculators.sigma(U, C=C, lame=lame),
#               stressStrainCalculators.epsilon(U))*dx - dot(F,v)*dx - dot(T,v)*ds

    a = inner(stressStrainCalculators.sigma(u, C=C, lame=lame), 
              stressStrainCalculators.epsilon(v)) * dx
    L = dot(F,v)*dx + dot(T,v)*ds
    
    solverObj = solverObject('Elasticity Equation')
    solverObj.a = a
    solverObj.L = L
    solverObj.u = Function(V)
    solverObj.bc = bc
    
    return solverObj

def solve_Elasticity(solveObject):
    toc = time.time()
    if solveObject.bc == None:
#         solve(solveObject.F == 0, solveObject.u)
        solve(solveObject.a == solveObject.L, solveObject.u)
    else:
#         solve(solveObject.F == 0, solveObject.u, solveObject.bc)
        solve(solveObject.a == solveObject.L, solveObject.u, solveObject.bc)
    tic = time.time()
    
    solveMsg = 'Solved {}'.format(solveObject.name)
    tP.printSolveUpdate(solveMsg, tic - toc)

def setup_heatSolve(V, U0, dt, K=Constant(1.0), F=Constant(0.0), bc=None):
    
    v = TestFunction(V)
    U = Function(V)
    
    F = U*v*dx + K*dt*dot(grad(U),grad(v))*dx - (U0 + dt*F)*v*dx
    
    solverObj = solverObject('Heat Equation')
    solverObj.F = F
    solverObj.u = U
    solverObj.u0 = U0
    solverObj.bc = bc
    
    return solverObj

def solve_Heat(solveObject):
    toc = time.time()
    if solveObject.bc == None:
        solve(solveObject.F == 0, solveObject.u)
    else:
        solve(solveObject.F == 0, solveObject.u, solveObject.bc)
    tic = time.time()
    solveObject.u0.assign(solveObject.u)
    
    solveMsg = 'Solved {}'.format(solveObject.name)
    tP.printSolveUpdate(solveMsg, tic - toc)

def scalarWaveSolve(V,U0,U,dt,C=Constant(1.0),F=Constant(0.0),bc=None,prms={'linear_solver': 'gmres','preconditioner': 'hypre_amg', \
                             "krylov_solver":{"nonzero_initial_guess":True, \
                                               "absolute_tolerance": 1e-10,\
                                               "relative_tolerance":1e-10, \
                                               "monitor_convergence":True
                                             }}):

    u = TrialFunction(V)
    v = TestFunction(V)

    a = u*v*dx + dt*dt*C*C*inner(grad(u),grad(v))*dx - dt*dt*F*v*dx
    L = 2*U*v*dx - U0*v*dx

    if bc is not None:
        A, b = assemble_system(a, L, bc)
    else:
        A, b = assemble_system(a, L)

    u = Function(V)

    solve(A, u.vector(), b)
    U0.assign(U)
    U.assign(u)

    return U, U0

def setup_linearPoissonSolve(V, U, bc=None, F=Constant(0.0), G=Constant(0.0)):
    v = TrialFunction(V)
    F = dot(grad(U), grad(v))*dx - F*v*dx - G*v*ds
    
    solverObj = solverObject('Linear Poisson')
    solverObj.F = F
    solverObj.u = U
    solverObj.bc = bc
    
    return solverObj

def solve_LinearPoisson(solveObject):
    toc = time.time()
    if solveObject.bc == None:
        solve(solveObject.F == 0, solveObject.u)
    else:
        solve(solveObject.F == 0, solveObject.u, solveObject.bc)
    tic = time.time()
    
    solveMsg = 'Solved {}'.format(solveObject.name)
    tP.printSolveUpdate(solveMsg, tic - toc)

def setup_nonlinearPoisson(V, U, NLF=Constant(1.0), bc=None,
                           F=Constant(0.0), G=Constant(0.0),
                           prms={'linear_solver': 'gmres','preconditioner': 'hypre_amg', \
                        "krylov_solver":{"nonzero_initial_guess":True, \
                                          "absolute_tolerance": 1e-10,\
                                          "relative_tolerance":1e-10, \
                                          "monitor_convergence":True}}):

    #Define variation problem
    v   = TestFunction(V)
    F   = inner(NLF * grad(U), grad(v))*dx - F*v*dx

    solverObj           = solverObject('Nonlinear Poisson')
    solverObj.F         = F
    solverObj.params    = prms
    solverObj.u         = U
    solverObj.bc        = bc

    return solverObj

def solveNonlinearPoisson(solveObject):
    toc = time.time()
    if solveObject.bc == None:
        solve(solveObject.F == 0, solveObject.u)
    else:
        solve(solveObject.F == 0, solveObject.u, solveObject.bc)
    tic = time.time()
    solveMsg = 'Solved {}'.format(solveObject.name)
    tP.printSolveUpdate(solveMsg, tic - toc)

def setupScalarAllenCahnSolver(V, U0, U, DF, kineticCoef, dt, bc=None, args=None):

    class AllenCahnEquation(NonlinearProblem):
        def __init__(self, a, L):
            NonlinearProblem.__init__(self)
            self.L = L
            self.a = a
        def F(self, b, x):
            assemble(self.L, tensor=b)
        def J(self, A, x):
            assemble(self.a, tensor=A)

    du, v   = TrialFunction(V), TestFunction(V)
    
    if inspect.isfunction(DF):
        if args != None:
            DF = DF(U,*args)
        else:
            DF = DF(U)

    F       = U*v*dx - U0*v*dx + kineticCoef*(dt * dot(grad(U), grad(v))*dx + dt * DF * v * dx)
    J       = derivative(F, U, du)

    problem = AllenCahnEquation(J, F)
    solver  = NewtonSolver()

    solverObj           = solverObject('AC')
    solverObj.problem   = problem
    solverObj.solver    = solver
    solverObj.u         = U
    solverObj.u0        = U0
    solverObj.bc        = bc

    return solverObj

def solveScalarAllenCahn(solveObject):
    toc = time.time()
    if solveObject.bc == None:
        solveObject.solver.solve(solveObject.problem, solveObject.u.vector())
    else:
        solveObject.solver.solve(solveObject.problem, solveObject.u.vector())
    tic = time.time()
    solveMsg = 'Solved {}'.format(solveObject.name)
    tP.printSolveUpdate(solveMsg, tic - toc)
    solveObject.u0.vector()[:]      = solveObject.u.vector()
    
def solve_CahnHilliard(solveObject):
    toc = time.time()
    if solveObject.bc == None:
        solveObject.solver.solve(solveObject.problem, solveObject.u.vector())
    else:
        solveObject.solver.solve(solveObject.problem, solveObject.u.vector())
    tic = time.time()
    solveMsg = 'Solved {}'.format(solveObject.name)
    tP.printSolveUpdate(solveMsg, tic - toc)
    solveObject.u0.vector()[:] = solveObject.u.vector()
    
def setup_CahnHilliardSolver(ME, U, U0, dt, eps, DF, theta=0.5, bc=None):
    class CahnHilliardEquation(NonlinearProblem):
        def __init__(self, a, L):
            NonlinearProblem.__init__(self)
            self.L = L
            self.a = a
        def F(self, b, x):
            assemble(self.L, tensor=b)
        def J(self, A, x):
            assemble(self.a, tensor=A)
            
    du = TrialFunction(ME)
    v, q = TestFunction(ME)
    
    dc, dmu = split(du)
    c, mu = split(U)
    c0, mu0 = split(U0)
    
    mu_mid = (1.0 - theta)*mu0 + theta * mu
    
    F0 = c*q*dx - c0*q*dx + dt*dot(grad(mu_mid), grad(q))*dx
    F1 = mu*v*dx - DF*v*dx - eps * dot(grad(c), grad(v))*dx
    F = F0 + F1
    
    J       = derivative(F, U, du)

    problem = CahnHilliardEquation(J, F)
    solver  = NewtonSolver()

    solverObj           = solverObject('Cahn-Hilliard')
    solverObj.problem   = problem
    solverObj.solver    = solver
    solverObj.u         = U
    solverObj.u0        = U0
    solverObj.bc        = bc

    return solverObj
