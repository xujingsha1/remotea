from fenics import *

class Problem(NonlinearProblem):
    def __init__(self, J, J_pc, F, bcs):
        self.bilinear_form = J
        self.preconditioner_form = J_pc
        self.linear_form = F
        self.bcs = bcs
        NonlinearProblem.__init__(self)

    def F(self, b, x):
        pass

    def J(self, A, x):
        pass

    def form(self, A, P, b, x):
        assemble(self.linear_form, tensor=b)
        assemble(self.bilinear_form, tensor=A)
        assemble(self.preconditioner_form, tensor=P)
        if self.bcs != None :
            for bc in self.bcs:
                bc.apply(b, x)
                bc.apply(A)
                bc.apply(P)

class CustomSolver(NewtonSolver):
    def __init__(self,mesh):
        NewtonSolver.__init__(self, mesh.mpi_comm(),
                              PETScKrylovSolver(), PETScFactory.instance())

    def solver_setup(self, A, P, problem, iteration):
        self.linear_solver().set_operators(A, P)
        # self.linear_solver().set_operator(A)
        
        #1. mumps direct solver 
        # PETScOptions.set("ksp_type", "preonly") 
        # PETScOptions.set("pc_type", 'lu')
        # PETScOptions.set("pc_factor_mat_solver_type", 'mumps')

        #2. gmres, mumps preconditioner
        PETScOptions.set("ksp_type", "gmres") 
        # PETScOptions.set('pc_type', 'bjacobi')
        PETScOptions.set('pc_type', 'lu')
        PETScOptions.set("pc_factor_mat_solver_type", 'mumps')
        # to change to other preconditioner. e.g ilu, sor

        #3. gmres, hypre_amg preconditioner
        # PETScOptions.set("ksp_type", "gmres") 
        # PETScOptions.set("pc_type",'hypre')
        # PETScOptions.set('pc_hypre_type', 'boomeramg')  #boomeramgeuclid

        #4. gmres, other preconditioners
        # PETScOptions.set("ksp_type", "gmres") 
        # PETScOptions.set("pc_type",'ilu') #bjacobi, sor

        # PETScOptions.set('ksp_rtol', 1e-10)
        # PETScOptions.set('ksp_atol', 1e-10)
        PETScOptions.set('ksp_max_it', 800)
        #PETScOptions.set("ksp_monitor")
        PETScOptions.set('ksp_monitor_true_residual')
        PETScOptions.set('ksp_error_if_not_converged',0)
        # PETScOptions.set('ksp_view') # view solver and matrix info
        
        self.linear_solver().set_from_options()