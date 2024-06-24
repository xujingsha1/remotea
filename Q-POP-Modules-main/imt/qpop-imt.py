"""
Last modified: 2024-04-20
        Added the interface for block GS preoconditioner 
"""

from __future__ import print_function
from fenics import *
import math
import numpy as np
# from numpy.random import rand
import time
import xml.etree.ElementTree as ET
from customSolver import Problem, CustomSolver # user-defined class for custom solver 

start_time = time.time()
comm = MPI.comm_world

# Using mpi4py commands
rank = comm.Get_rank()  #  MPI.rank(comm)
psize = comm.Get_size()   # MPI.size(comm)

if rank == 0:
    list_linear_algebra_backends()
    list_linear_solver_methods()

print(f'User message ===> Processors report: rank/size = {rank:5d}/{psize:5d}', flush=True)

# Set some global parameters
parameters["mesh_partitioner"] = "ParMETIS"  # Default SCOTCH might have some side effects and bugs
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True

# Define parameter and physical parameter classes
class Parameter:
    def __init__(self, value):
        self.value = value

    def read(self, fileroot, fileroute, unit=1.0, comm=None, rank=None):
        isroot = True
        bfind = False
        if comm is not None:
            isroot = (rank == 0)
        if isroot:
            routelist = fileroute.replace(' ', '').split('/')
            nroute = len(routelist)
            route0 = fileroot
            bfind = True
            routename = ''
            for i in range(nroute - 1):
                route0 = route0.find(routelist[i])
                routename = routename + routelist[i] + '/'
                if route0 is None:
                    print('User message ===> No', routename[:-1], 'specified, using default value(s)', flush=True)
                    bfind = False
                    break
                
            if bfind:
                if '.' in routelist[nroute - 1]:
                    r = routelist[nroute - 1].split('.')
                    route0 = route0.find(r[0])
                    if route0 is None:
                        print('User message ===> No', routename + r[0], 'specified, using default value(s)', flush=True)
                        bfind = False
                    elif r[1] in route0.attrib:
                        intext = route0.get(r[1])
                    else:
                        print('User message ===> No', routename + routelist[nroute - 1], 'specified, using default value(s)', flush=True)
                        bfind = False
                else:
                    route0 = route0.find(routelist[nroute - 1])
                    if route0 is None:
                        print('User message ===> No', routename + routelist[nroute - 1], 'specified, using default value(s)', flush=True)
                        bfind = False
                    else:
                        intext = route0.text

            if bfind:
                if type(self.value) is str:
                    self.value = intext
                elif type(self.value) is int:
                    self.value = int(intext)
                elif type(self.value) is float:
                    self.value = float(intext) / unit
                elif type(self.value) is bool:
                    if intext == 'True':
                        self.value = True
                    elif intext == 'False':
                        self.value = False
                    else:
                        raise ValueError('Not inputing a valid boolean value')
                else:
                    raise TypeError('Not inputing a valid data type')
                print('User message ===> Input', fileroute, 'is', self.value, flush=True)

        if comm is not None:
            bfind = comm.bcast(bfind, root=0)
            if bfind:
                self.value = comm.bcast(self.value, root=0)
                if isroot:
                    print('User message ===> Broadcasted', fileroute, flush=True)
                # Test
                # if not isroot:
                #     print('!!!!!!!!!!!!!!!!! Rank', str(rank)+'\'s', fileroute, 'is', self.value, flush=True)
                
class PhysParam(Parameter, Constant):
    def __init__(self, value):
        Parameter.__init__(self, float(value))
        Constant.__init__(self, float(value))

    def read(self, fileroot, fileroute, unit=1.0, comm=None, rank=None):
        Parameter.read(self, fileroot, fileroute, unit, comm, rank)   # Calling class method, so instance self needs to be explicitly passed to the function
        self.assign(Constant(self.value))   # Calling instance method, so instance self is automatically passed to the function
        # Test
        # print('!!!!!!!!!!!!!!!!! Rank', str(rank)+'\'s', fileroute, 'is', self.value, self.values(), flush=True)

    def update(self, value):
        self.value = float(value)
        self.assign(Constant(self.value))


#-----------------------------------------------
# Define units (in SI units) to renormalize quantities
#-----------------------------------------------
LUNIT = 1e-9
TUNIT = 1e-9
TEMPUNIT = 338.0
EUNIT = 1.3806504e-23 * TEMPUNIT
VUNIT = 1e-3
CUNIT = EUNIT/VUNIT
MUNIT = EUNIT*(TUNIT/LUNIT)**2
RUNIT = VUNIT/(CUNIT/TUNIT)


intree = ET.parse('input.xml')
inroot = intree.getroot()
#-----------------------------------------------
# Define parameters for geometry and time
#-----------------------------------------------
Lx = PhysParam(50e-9/LUNIT)
Lx.read(inroot, 'external/Lx', unit=LUNIT/1e-9, comm=comm, rank=rank)
Ly = PhysParam(20e-9/LUNIT)
Ly.read(inroot, 'external/Ly', unit=LUNIT/1e-9, comm=comm, rank=rank)
nx = Parameter(50)   # 96
nx.read(inroot, 'external/Lx.mesh', comm=comm, rank=rank)
ny = Parameter(20)   # 32
ny.read(inroot, 'external/Ly.mesh', comm=comm, rank=rank)
Lz = PhysParam(36e-9/LUNIT)
Lz.read(inroot, 'external/Lz', unit=LUNIT/1e-9, comm=comm, rank=rank)

tf = Parameter(5000e-9/TUNIT)         # Final time in nanoseconds  5000.0
tf.read(inroot, 'time/endtime', unit=TUNIT/1e-9, comm=comm, rank=rank)
dt = PhysParam(1e-6)       # Initial time step
savemethod = Parameter('fix')  # or 'auto'
savemethod.read(inroot, 'time/savemethod', comm=comm, rank=rank)
saveperiod = Parameter(5.0)
saveperiod.read(inroot, 'time/saveperiod', comm=comm, rank=rank)
t_out = np.arange(0, tf.value+1e-12*tf.value, saveperiod.value) # np.linspace(0, tf, 201)  # Output time sequence
len_t_out = len(t_out)

alpha = Parameter(1e-6) #epsilon = alpha * h
alpha.read(inroot, 'solverparameters/Nitschefactor', comm=comm, rank=rank)

#adaptive time stepping parameters
tol_tstep = Parameter(0.01)  # Goal relative error for time-stepping
tol_tstep.read(inroot, 'solverparameters/timesteptolerance', comm=comm, rank=rank)
r_t_max = Parameter(4.0)   # 5.0
r_t_min = Parameter(0.2)
s_t = Parameter(0.9)
max_div = Parameter(100)  # Allowed maximum number of successive divergence

#newton solver parameters
rel_par = Parameter(0.9)
rel_par.read(inroot, 'solverparameters/Newtonrelaxation', comm=comm, rank=rank)
max_iter = Parameter(30)
max_iter.read(inroot, 'solverparameters/Newtonmaxiteration', comm=comm, rank=rank)
# kr_solver_iter = 2000
newton_abs_tol = Parameter(1e-6)
newton_abs_tol.read(inroot, 'solverparameters/Newtonabsolutetolerance', comm=comm, rank=rank)
newton_rel_tol = Parameter(1e-3)   #  1E-2
newton_rel_tol.read(inroot, 'solverparameters/Newtonrelativetolerance', comm=comm, rank=rank)

# linear solver parameters
use_GS_block_preconditioner = Parameter(0)
use_GS_block_preconditioner.read(inroot, 'solverparameters/useGSblockpreconditioner', comm=comm, rank=rank)  

#not tuned parameters
# num_steps = 1000     # Number of time steps
# dt_out = tf / num_steps  # Time interval for output

# tol = 1E-14   # Tolerance for geometry identifying

directsolver = Parameter('superlu_dist')
directsolver.read(inroot, 'solverparameters/directsolver', comm=comm, rank=rank)
fenicslog = Parameter('INFO')
fenicslog.read(inroot, 'solverparameters/loglevel', comm=comm, rank=rank)

#-------------------------------------------
# Define physical parameters
#-------------------------------------------
KB = PhysParam(1.3806504e-23/EUNIT*TEMPUNIT)
ECHARGE = PhysParam(1.602176487e-19/CUNIT)
UCVOL = PhysParam(59e-30/LUNIT**3)
UCVOL.read(inroot, 'internal/unitcellvol', unit=LUNIT**3/1e-27, comm=comm, rank=rank)
TC = PhysParam(338.0/TEMPUNIT)
TC.read(inroot, 'internal/Tc', unit=TEMPUNIT, comm=comm, rank=rank)
T1 =  PhysParam(275.0/TEMPUNIT)   # 0.81361
T1.read(inroot, 'internal/T1', unit=TEMPUNIT, comm=comm, rank=rank)
AN1 = PhysParam(2.05714*KB.value*TC.value/UCVOL.value)  # 34.867
AN1.read(inroot, 'internal/a1', unit=UCVOL.value/(KB.value*TC.value), comm=comm, rank=rank)
AN2 = PhysParam((-0.623108 + 0.121228)/2*KB.value*TC.value/UCVOL.value)
AN2.read(inroot, 'internal/a2', unit=UCVOL.value/(KB.value*TC.value), comm=comm, rank=rank)
AN3 = PhysParam((0.330568 + 4.18947)/4*KB.value*TC.value/UCVOL.value)
AN3.read(inroot, 'internal/a3', unit=UCVOL.value/(KB.value*TC.value), comm=comm, rank=rank)
T2 = PhysParam(270.0/TEMPUNIT)
T2.read(inroot, 'internal/T2', unit=TEMPUNIT, comm=comm, rank=rank)
AU1 = PhysParam(3.94286*KB.value*TC.value/UCVOL.value)
AU1.read(inroot, 'internal/b1', unit=UCVOL.value/(KB.value*TC.value), comm=comm, rank=rank)
AU2 = PhysParam((1.36767 - 3.67915)/2*KB.value*TC.value/UCVOL.value)
AU2.read(inroot, 'internal/b2', unit=UCVOL.value/(KB.value*TC.value), comm=comm, rank=rank)
AU3 = PhysParam((0.4 + 2)/4*KB.value*TC.value/UCVOL.value)
AU3.read(inroot, 'internal/b3', unit=UCVOL.value/(KB.value*TC.value), comm=comm, rank=rank)
GNU1 = PhysParam(0.3*KB.value*TC.value/UCVOL.value)
GNU1.read(inroot, 'internal/c1', unit=UCVOL.value/(KB.value*TC.value), comm=comm, rank=rank)
GNU2 = PhysParam((0.2 - 1.5 + 0.3/2)/2*KB.value*TC.value/UCVOL.value)
GNU2.read(inroot, 'internal/c2', unit=UCVOL.value/(KB.value*TC.value), comm=comm, rank=rank)
GNU3 = PhysParam((0.05 + 2.0)/2*KB.value*TC.value/UCVOL.value)
GNU3.read(inroot, 'internal/c3', unit=UCVOL.value/(KB.value*TC.value), comm=comm, rank=rank)
CHI = PhysParam(0.286*1.602176487e-19/EUNIT)
CHI.read(inroot, 'internal/gapcoeff', unit=EUNIT/1.602176487e-19, comm=comm, rank=rank)
tau_latt = PhysParam(1e-12/TUNIT)
tau_latt.read(inroot, 'internal/SPTtimeconst', unit=TUNIT/1e-9, comm=comm, rank=rank)
KN = PhysParam(1.0/(tau_latt.value*AN1.value*(TC.value-T1.value)/TC.value))
KAPPAN = PhysParam(1.0*(1.602176487e-19/1E-9)/EUNIT*LUNIT)
KAPPAN.read(inroot, 'internal/SOPgradcoeff', unit=(EUNIT/1.602176487e-19)/(LUNIT/1E-9), comm=comm, rank=rank)
tau_elec = PhysParam(10e-15/TUNIT)
tau_elec.read(inroot, 'internal/IMTtimeconst', unit=TUNIT/1e-9, comm=comm, rank=rank)
nexcite = PhysParam(0.16/UCVOL.value)    # Density of excited electrons to close the gap
nexcite.read(inroot, 'internal/IMTtimeconst.excitdensity', unit=UCVOL.value, comm=comm, rank=rank)
KU = PhysParam(1.0/(tau_elec.value*CHI.value*2*nexcite.value))
KAPPAU = PhysParam(1.0*(1.602176487e-19/1E-9)/EUNIT*LUNIT)
KAPPAU.read(inroot, 'internal/EOPgradcoeff', unit=(EUNIT/1.602176487e-19)/(LUNIT/1E-9), comm=comm, rank=rank)
NC = PhysParam(2*( 65*(9.10938215e-31/MUNIT) * KB.value*TC.value / (2*math.pi * (1.054571628e-34/EUNIT/TUNIT)**2) )**1.5)
NC.read(inroot, 'internal/conductbandeffDOS', unit=UCVOL.value, comm=comm, rank=rank)
NV = PhysParam(NC.value)
NC.read(inroot, 'internal/valencebandeffDOS', unit=UCVOL.value, comm=comm, rank=rank)

MEC = PhysParam(5e-5/LUNIT**2*VUNIT*TUNIT)
MEC.read(inroot, 'internal/elecmobilityy', unit=LUNIT**2/(VUNIT*TUNIT*1e-4), comm=comm, rank=rank)
MEA = PhysParam(MEC.value*0.5)
MEA.read(inroot, 'internal/elecmobilityx', unit=LUNIT**2/(VUNIT*TUNIT*1e-4), comm=comm, rank=rank)
MHC = PhysParam(MEC.value/1.2)
MHC.read(inroot, 'internal/holemobilityy', unit=LUNIT**2/(VUNIT*TUNIT*1e-4), comm=comm, rank=rank)
MHA = PhysParam(MHC.value*0.5)
MHA.read(inroot, 'internal/holemobilityx', unit=LUNIT**2/(VUNIT*TUNIT*1e-4), comm=comm, rank=rank)
KEH0 = PhysParam(1.0/( 2 * math.sqrt(NC.value*NV.value)*math.exp(-CHI.value*0.827987**2/(KB.value*322.0/TEMPUNIT)) * (14.235e-6/TUNIT) ))  # 1.0/( 2 * math.sqrt(NC*NV)*math.exp(-CHI/(2*KB*TC)) * (14e-6/TUNIT) )
KEH0.read(inroot, 'internal/ehrecrate', unit=(1e-9/TUNIT)/UCVOL.value, comm=comm, rank=rank)
rel_perm = PhysParam(60.0)
rel_perm.read(inroot, 'internal/relativepermitt', comm=comm, rank=rank)
PERMITTIVITY = PhysParam(rel_perm.value*(8.854187817e-12*VUNIT/CUNIT*LUNIT))
CPV = PhysParam(690*4340/EUNIT*LUNIT**3*TEMPUNIT)
CPV.read(inroot, 'internal/volheatcapacity', unit=EUNIT/(TEMPUNIT*LUNIT**3), comm=comm, rank=rank)
THETA = PhysParam(6/EUNIT*TUNIT*LUNIT*TEMPUNIT)
THETA.read(inroot, 'internal/thermalconduct', unit=EUNIT/(TUNIT*LUNIT*TEMPUNIT), comm=comm, rank=rank)
HTRAN = PhysParam(3e6/EUNIT*TUNIT*LUNIT**2*TEMPUNIT)
HTRAN.read(inroot, 'external/heatdiss', unit=EUNIT/(TUNIT*LUNIT**2*TEMPUNIT), comm=comm, rank=rank)
# ETA_IN = 0
# MU_IN = 0
# ETA_IN = 1.0
# MU_IN = -1.0
# Lian changed
CHP_IN = PhysParam(0.0)   # Intrinsic chemical potential

#----------------------------------------------------------------
# Define boundary and initial values
#----------------------------------------------------------------
Ts = PhysParam(300.0 / TEMPUNIT)
Ts.read(inroot, 'external/temperature', unit=TEMPUNIT, comm=comm, rank=rank)
#Lian change
#Ts = 300 / 338.0             # K / [Chosen temperature unit]
delV = PhysParam(0.07 / VUNIT)                  # V
delV.read(inroot, 'external/voltage', unit=VUNIT, comm=comm, rank=rank)
Resistor = PhysParam(8.0E3 / RUNIT)  # Ohm / [Chosen resistance unit]
Resistor.read(inroot, 'external/resistor', unit=RUNIT, comm=comm, rank=rank)
Capacitor = PhysParam(1e-9 / (CUNIT/VUNIT))
Capacitor.read(inroot, 'external/capacitor', unit=CUNIT/VUNIT/1e-9, comm=comm, rank=rank)
tramp = Parameter(10e-9 / TUNIT)
tramp.read(inroot, 'time/rampt', unit=TUNIT/1e-9, comm=comm, rank=rank)
etas = PhysParam(1.0)
mus = PhysParam(-1.0)
deff = PhysParam(10e-9 / LUNIT)

T_i = PhysParam(Ts.value)   # TC
T_i.read(inroot, 'initialization/temperature', unit=TEMPUNIT, comm=comm, rank=rank)
eta_i = PhysParam(0.791296*math.sqrt(2))   # 1.0
eta_i.read(inroot, 'initialization/SOP', comm=comm, rank=rank)
mu_i = PhysParam(-0.914352*math.sqrt(2))   # -1.0
mu_i.read(inroot, 'initialization/EOP', comm=comm, rank=rank)
phi_i = PhysParam(-1000.0)  # negative here just to raise a flag. See below.
phi_i.read(inroot, 'initialization/voltage', unit=VUNIT, comm=comm, rank=rank)

Tcvarmethod = Parameter('random')
Tcvarmethod.read(inroot, 'initialization/Tcvariance.method', comm=comm, rank=rank)
sigma = Parameter(0.1 * TC.value)
sigma.read(inroot, 'initialization/Tcvariance/sigma', unit=TEMPUNIT, comm=comm, rank=rank)
Tcvar0 = Parameter(2.0 * sigma.value)
Tcvar0.read(inroot, 'initialization/Tcvariance/mean', unit=TEMPUNIT, comm=comm, rank=rank)
corr_len = Parameter(0.2 * Ly.value)
corr_len.read(inroot, 'initialization/Tcvariance/correlationlength', unit=LUNIT/1e-9, comm=comm, rank=rank)
rseed = Parameter(11793)
rseed.read(inroot, 'initialization/Tcvariance/randomseed', comm=comm, rank=rank)
TCIMP0 = Parameter(-20.0 / TEMPUNIT)
TCIMP0.read(inroot, 'initialization/Tcvariance/Tcshift', unit=TEMPUNIT, comm=comm, rank=rank)
RIMPDIS = Parameter(3e-9 / LUNIT)
RIMPDIS.read(inroot, 'initialization/Tcvariance/radius', unit=LUNIT/1e-9, comm=comm, rank=rank)
bump = Parameter(0.02)
bump.read(inroot, 'initialization/Tcvariance/bump', comm=comm, rank=rank)

#filenames
fig_file = 'solution/'
allin1file = 'log.txt'

if savemethod.value == 'fix':
    t_out = np.arange(0, tf.value+1e-12*tf.value, saveperiod.value)
    len_t_out = len(t_out)
elif savemethod.value == 'auto':
    saveperiod.value = int(saveperiod.value)

gamma_ei = PhysParam(-(CHI.value * mu_i.value**2/2 - CHP_IN.value)/(KB.value * T_i.value))
gamma_hi = PhysParam(-(CHI.value * mu_i.value**2/2 + CHP_IN.value)/(KB.value * T_i.value))
RVO2_i = Ly.value/( ECHARGE.value * (NC.value * math.exp(gamma_ei.value) * MEC.value + NV.value * math.exp(gamma_hi.value) * MHC.value) * Lx.value * Lz.value )
Vfrac_i = RVO2_i/(RVO2_i + Resistor.value)

if phi_i.value < -500:
    phi_i.update(delV.value * Vfrac_i)   # This is default phi_i
    delV_i = delV.value
else:
    delV_i = phi_i.value / Vfrac_i

delVr = PhysParam(delV_i)   # This is the ramping voltage
integral_phi_i = PhysParam(phi_i.value / RVO2_i / Lz.value)
# print('!!!!!!!!!!!!! Vfrac_i = ', Vfrac_i, '!!!!!!!!!!!!!', flush=True)

# Ramp up delV
def ramped_delV(t):
    if t <= tramp.value:
        _delVr = delV_i + (delV.value - delV_i) * t / tramp.value
    else:
        _delVr = delV.value
    return _delVr

comm.Barrier()

#----------------------------------------------------
# Generate mesh
#----------------------------------------------------
mesh = RectangleMesh(Point(0, 0), Point(Lx.value, Ly.value), nx.value, ny.value, 'crossed')
meshdim = mesh.topology().dim()

#--------------------------------------------------------------
# Define function space and functions for multiple variables
#--------------------------------------------------------------
P1 = FiniteElement('Lagrange', triangle, 1)
R0 = FiniteElement('R', triangle, 0)
V1 = FunctionSpace(mesh, P1)
VR = FunctionSpace(mesh, R0)

# n_f = 9  # Number of functions that need to be solved
# First 6 are for unknown variables and last 3 are for Lagrange multipliers
element = MixedElement([P1,P1,P1,P1,P1,P1,R0])
V = FunctionSpace(mesh, element)

# Define test functions
v_1, v_2, v_3, v_4, v_5, v_6, v_10 = TestFunctions(V)


# Define trial functions
du = TrialFunction(V)   # Newton iteration step

# Define functions
u = Function(V)     # The most recently computed solution
# Here we don't bother including in the previous solution only the 
# dynamical functions (that need update), because this will make updating previous solution 
# convenient and won't require other temporary functions to store the solution 
u_n = Function(V)  # The previous solution

# Split system functions to access components
eta, mu, gamma_e, gamma_h, phi, T, integral_phi = split(u)
eta_n, mu_n, gamma_en, gamma_hn, phi_n, T_n, integral_phin = split(u_n)

#------------------------------------------------------
# Assign initial values to the functions u_n and initial
# guess to u
#------------------------------------------------------
eta_0 = project(eta_i, V1)
mu_0 = project(mu_i, V1)
gamma_e0 = project(gamma_ei, V1)
gamma_h0 = project(gamma_hi, V1)
phi_0 = project(Expression('phi_i - phi_i*x[1]/Ly', degree=1, phi_i=phi_i.value, Ly=Ly.value), V1)
T_0 = project(T_i, V1)
# lam_10 = project(Constant(0), V1)
# lam_20 = project(Constant(0), V1)
# lam_30 = project(Constant(0), V1)
integral_phi_0 = project(integral_phi_i, VR)

# lam_1n = project(Expression('x[1]>=tol && Ly-x[1]>=tol ? 0.0 : 0.1', \
#                             degree=1, tol=tol, Ly=Ly), V1)
# lam_2n = project(Expression('x[1]>=tol && Ly-x[1]>=tol ? 0.0 : 0.1', \
#                             degree=1, tol=tol, Ly=Ly), V1)
# lam_3n = project(Expression('x[1]>=tol ? 0.0 : 0.1', \
#                             degree=1, tol=tol), V1)

fa = FunctionAssigner(V, [V1,V1,V1,V1,V1,V1,VR])
fa.assign(u_n, [eta_0, mu_0, gamma_e0, gamma_h0, phi_0, T_0, integral_phi_0])
u.assign(u_n)

# Set the Tc variation.
if Tcvarmethod.value == 'random':
    if rank == 0:
        print('User message ===> Generating correlated random Tc variance field...', flush=True)
    # coord = mesh.coordinates()  # Obtain local mesh coordinates which may contain "ghost dof"
    coordflat = V1.tabulate_dof_coordinates().flatten()   # Obtain local coordinates ordered by their local dof  (no ghost dof)
    ncoordflat = len(coordflat)
    # print(coordflat)
    # ncoordtotal = int(MPI.sum(comm, ncoord))
    sendcounts = np.array(comm.gather(ncoordflat, root=0))

    if rank == 0:
        # print("sendcounts: {}, total: {}".format(sendcounts, sum(sendcounts)))
        coordflat_g = np.empty(sum(sendcounts), dtype=np.float64)
    else:
        coordflat_g = None
    # coord_g = np.array([[(i%(nx+1))*(Lx/nx), (i//(nx+1))*(Ly/ny)] for i in range((nx+1)*(ny+1))])
    comm.Gatherv(sendbuf=coordflat, recvbuf=(coordflat_g, sendcounts), root=0)   # Gather the global coordinates
    # Calculate the correlated field in processor 0
    if rank == 0:
        coord_g = coordflat_g.reshape((-1, meshdim))
        ncoord_total = len(coord_g)
        # print(ncoord_total, 'coordinates:')
        # print(coord_g)
        # Set the global covariance matrix.
        cov_g = [[sigma.value**2*math.exp(-np.linalg.norm(coord_g[i]-coord_g[j])/corr_len.value) for i in range(ncoord_total)] for j in range(ncoord_total)]  
        rng = np.random.default_rng(rseed.value)   # Fix seed
        Tcvar_g = rng.multivariate_normal(np.ones(ncoord_total)*Tcvar0.value, cov_g)  # Obtain the global Tc variation field in numpy array
        maxTcvar = np.amax(Tcvar_g)
        Tcvar_g = Tcvar_g - maxTcvar
        # Tcvar_g = np.array([coord_g[i,0]+coord_g[i,1] for i in range(ncoord_total)])  # Test
    else:
        Tcvar_g = None

    # Scatter correlated field back to each processor
    Tcvar_ = np.empty(ncoordflat//meshdim, dtype=np.float64)
    if rank == 0:
        sendcounts2 = sendcounts // meshdim
    else:
        sendcounts2 = np.zeros(psize, dtype=int) # Initialize sendcounts2 on worker processors
    comm.Scatterv(sendbuf=(Tcvar_g, sendcounts2), recvbuf=Tcvar_, root=0)
    # print(Tcvar_)

    # dm = V1.dofmap()
    # local_range = dm.ownership_range()
    # local_dim = local_range[1] - local_range[0]

    # coord_test = dm.tabulate_all_coordinates(mesh).reshape((-1, meshdim))  # New fenics version no longer supports tabulate_all_coordinates()

    Tcvar = Function(V1)

    # Get mapping from dof to vertices. From manual: Only works for FunctionSpace with dofs exclusively on vertices. 
    # For MixedFunctionSpaces vertex index is offset with the number of dofs per vertex. 
    # In parallel the returned map only maps local (to processor) dofs.
    # d2v = dof_to_vertex_map(V1)  
    # v2d = vertex_to_dof_map(V1)
    # print('!!!!!!!!!!!', ncoord, local_range, local_dim, len(Tcvar.vector().get_local()), len(v2d), '!!!!!!!!!!!')
    # Tcvar.vector().set_local(Tcvar_[d2v])  # Set the Function values. Must use .set_local() for parallel run

    # gathertest = np.array([rank, rank+1])
    # gathertest_new = comm.gather(gathertest, root=0)

    Tcvar.vector().set_local(Tcvar_)   # Set the function local dof values
    Tcvar.vector().apply("insert")

    if rank == 0:
        print('User message ===> Done.', flush=True)

elif Tcvarmethod.value == 'nucleus1':
    # Set a small region of impurity where the transition temperature is lower than pristine case
    Tcvar = Expression('TCIMP0*(-tanh(2*(sqrt(pow(x[0] - Lx/2, 2) + pow(x[1], 2)) - RIMPDIS)/DWWIDTH) + 1)/2', \
                       degree=1, Lx=Lx.value, TCIMP0=TCIMP0.value, RIMPDIS=RIMPDIS.value, DWWIDTH=10e-9/LUNIT)
else:
    # Set a small region of impurity where the transition temperature is lower than pristine case
    Tcvar = Expression('TCIMP0*(-tanh(2*(sqrt(pow(x[0] - Lx/2, 2) + pow(x[1], 2) + bump*(x[0] - Lx/2)*x[1]) - RIMPDIS)/DWWIDTH) + 1)/2', \
                       degree=1, Lx=Lx.value, TCIMP0=TCIMP0.value, RIMPDIS=RIMPDIS.value, DWWIDTH=10e-9/LUNIT, bump=bump.value)


vtkfile_Tcvar = File(fig_file+'Tcvar.pvd') # Plot and save the Tc variance
_Tcvar = project(Tcvar*Constant(TEMPUNIT), V1)
_Tcvar.rename("Tcvar","Tcvar")
vtkfile_Tcvar << _Tcvar

comm.Barrier()

#-----------------------------------------------------------------
# Mark boundaries and define normal boundary conditions
#-----------------------------------------------------------------
boundary_markers = MeshFunction('size_t', mesh, meshdim - 1)
boundary_markers.set_all(9)

class BoundaryX0(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0)
bx0 = BoundaryX0()
bx0.mark(boundary_markers, 0)   # 0 marks y = 0 boundary

class BoundaryX1(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], Ly.value)
bx1 = BoundaryX1()
bx1.mark(boundary_markers, 1)   # 1 marks y = Ly boundary

ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)  # Measure ds

bc_phi_1 = DirichletBC(V.sub(4), Constant(0), bx1)  # Dirichlet boundary condition for phi

# def domain_lam_1_2(x):
#     return x[1] > tol and x[1] < Ly - tol

# lam_1_constraint = DirichletBC(V.sub(6), Constant(0), domain_lam_1_2) # Set unused lam_1 to zero
# lam_2_constraint = DirichletBC(V.sub(7), Constant(0), domain_lam_1_2) # Set unused lam_2 to zero

# def domain_lam_3(x):
#     return x[1] > tol

# lam_3_constraint = DirichletBC(V.sub(6), Constant(0), domain_lam_3) # Set unused lam_3 to zero

#bcs = [bc_phi_1, lam_1_constraint, lam_2_constraint, lam_3_constraint] # Collect boundary conditions
bcs = [bc_phi_1] # Nitsche used for gamma_e,gamma_h

#-------------------------------------------------------------------
# Define python functions used in the weak form  
#-------------------------------------------------------------------
def dfb_deta(T, eta, mu):
    return AN1*(T - T1 - Tcvar)/TC*eta + AN2*eta**3 + AN3*eta**5 + GNU1*mu \
    - GNU2*eta*mu**2 + 1.5*GNU3*eta**2*mu

def dfb_dmu(T, eta, mu):
    return AU1*(T - T2 - Tcvar)/TC*mu + AU2*mu**3 + AU3*mu**5 + GNU1*eta \
    - GNU2*eta**2*mu + 0.5*GNU3*eta**3

def Fermi(gamma):
    return 1 / (exp(-gamma) + 3*sqrt(DOLFIN_PI)/4*(4 + gamma**2)**(-0.75))

def dFermi(gamma):   # Derivative of Fermi(gamma)
    return Fermi(gamma)**2 * (exp(-gamma) + 9*sqrt(DOLFIN_PI)/8 \
    *(4 + gamma**2)**(-1.75)*gamma)

# Branched version of Fermi and dFermi
def Fermi_b(gamma):
    return conditional(lt(gamma, -100), 0, conditional(gt(gamma, 100), \
                4/(3*sqrt(DOLFIN_PI))*(4 + gamma**2)**0.75, Fermi(gamma)))

def dFermi_b(gamma):
    return conditional(lt(gamma, -100), 0, conditional(gt(gamma, 100), \
                2/sqrt(DOLFIN_PI)*gamma*(4 + gamma**2)**(-0.25), dFermi(gamma)))

# gam5half = Constant(math.gamma(5.0/2))
# # Approximate Fermi integral of order 3/2
# def Fermi3half(gamma):
#     a = Constant(14.9)
#     b = Constant(2.64)
#     c = Constant(9.0/4)
#     return 1 / ( 2.5 * 2**2.5 / ( b + gamma + ( abs( gamma - b )**c + a )**( 1 / c ) )**2.5 \
#             + exp( -gamma ) / gam5half ) / gam5half

# n = FacetNormal(mesh)

#----------------------------------------------------------------
# Define intermediate variables for the variational problem
#----------------------------------------------------------------
nd_e = NC*Fermi(gamma_e)    # Electron density
nd_h = NV*Fermi(gamma_h)    # Hole density
# nd_e0 = NC*Fermi(gamma_e0)  # Corresponding initial values
# nd_h0 = NV*Fermi(gamma_h0)

# nd_en = NC*Fermi_b(gamma_en)    # Electron density
# nd_hn = NV*Fermi_b(gamma_hn)    # Hole density

j_ex = -nd_e*(MEA/ECHARGE)*(KB*T*gamma_e.dx(0) + gamma_e*KB*T.dx(0) \
                          + CHI*mu*mu.dx(0) - ECHARGE*phi.dx(0))
j_ey = -nd_e*(MEC/ECHARGE)*(KB*T*gamma_e.dx(1) + gamma_e*KB*T.dx(1) \
                          + CHI*mu*mu.dx(1) - ECHARGE*phi.dx(1))
j_e = as_vector([j_ex, j_ey])   # Electron flux
# Corresponding initial values
# j_ex0 = -nd_e0*(MEA/ECHARGE)*(KB*T_0*gamma_e0.dx(0) + gamma_e0*KB*T_0.dx(0) \
#                           + CHI*mu_0*mu_0.dx(0) - ECHARGE*phi_0.dx(0))
# j_ey0 = -nd_e0*(MEC/ECHARGE)*(KB*T_0*gamma_e0.dx(1) + gamma_e0*KB*T_0.dx(1) \
#                           + CHI*mu_0*mu_0.dx(1) - ECHARGE*phi_0.dx(1))

# j_exn = -nd_en*(MEA/ECHARGE)*(KB*T_n*gamma_en.dx(0) + gamma_en*KB*T_n.dx(0) \
#                           + CHI*mu_n*mu_n.dx(0) - ECHARGE*phi.dx(0))
# j_eyn = -nd_en*(MEC/ECHARGE)*(KB*T_n*gamma_en.dx(1) + gamma_en*KB*T_n.dx(1) \
#                           + CHI*mu_n*mu_n.dx(1) - ECHARGE*phi.dx(1))
# j_en = as_vector([j_exn, j_eyn])   # Electron flux



j_hx = -nd_h*(MHA/ECHARGE)*(KB*T*gamma_h.dx(0) + gamma_h*KB*T.dx(0) \
                          + CHI*mu*mu.dx(0) + ECHARGE*phi.dx(0))
j_hy = -nd_h*(MHC/ECHARGE)*(KB*T*gamma_h.dx(1) + gamma_h*KB*T.dx(1) \
                          + CHI*mu*mu.dx(1) + ECHARGE*phi.dx(1))
j_h = as_vector([j_hx, j_hy])   # Hole flux
# Corresponding initial values
# j_hx0 = -nd_h0*(MHA/ECHARGE)*(KB*T_0*gamma_h0.dx(0) + gamma_h0*KB*T_0.dx(0) \
#                           + CHI*mu_0*mu_0.dx(0) + ECHARGE*phi_0.dx(0))
# j_hy0 = -nd_h0*(MHC/ECHARGE)*(KB*T_0*gamma_h0.dx(1) + gamma_h0*KB*T_0.dx(1) \
#                           + CHI*mu_0*mu_0.dx(1) + ECHARGE*phi_0.dx(1))

# j_hxn = -nd_hn*(MHA/ECHARGE)*(KB*T_n*gamma_hn.dx(0) + gamma_hn*KB*T_n.dx(0) \
#                           + CHI*mu_n*mu_n.dx(0) + ECHARGE*phi.dx(0))
# j_hyn = -nd_hn*(MHC/ECHARGE)*(KB*T_n*gamma_hn.dx(1) + gamma_hn*KB*T_n.dx(1) \
#                           + CHI*mu_n*mu_n.dx(1) + ECHARGE*phi.dx(1))
# j_hn = as_vector([j_hxn, j_hyn])   # Hole flux

nd_in = NC*Fermi(-(CHI*mu**2/2 - CHP_IN)/(KB*T))   # Intrinsic carrier density
nd_eeq = NC*Fermi(-(CHI*mu**2/2 - ECHARGE*phi - CHP_IN)/(KB*T))  # Equilibrium electron density
nd_heq = NV*Fermi(-(CHI*mu**2/2 + ECHARGE*phi + CHP_IN)/(KB*T))  # Equilibrium hole density
# nd_in0 = NC*Fermi(-(CHI*mu_0**2/2 - CHP_IN)/(KB*T_0))   #  Corresponding initial values
# nd_eeq0 = NC*Fermi(-(CHI*mu_0**2/2 - ECHARGE*phi_0 - CHP_IN)/(KB*T_0))  
# nd_heq0 = NV*Fermi(-(CHI*mu_0**2/2 + ECHARGE*phi_0 + CHP_IN)/(KB*T_0)) 

# nd_inn = NC*Fermi_b((-CHI*mu_n**2/2 + CHP_IN)/(KB*T_n))   # Intrinsic carrier density
# nd_eeqn = NC*Fermi_b((-CHI*mu_n**2/2 + ECHARGE*phi + CHP_IN)/(KB*T_n))  # Equilibrium electron density
# nd_heqn = NV*Fermi_b((-CHI*mu_n**2/2 - ECHARGE*phi - CHP_IN)/(KB*T_n))  # Equilibrium hole density

Jy = ECHARGE*(j_hy - j_ey)    # y component of the total current
# Jy_n = ECHARGE*(j_hyn - j_eyn)    # y component of the total current


# Jy_testn = -NV*dFermi_b(gamma_h)*MHC*(KB*T*gamma_h.dx(1) + gamma_h*KB*T.dx(1) \
#            + CHI*mu*mu.dx(1) + ECHARGE*phi.dx(1)) \
#            -nd_h*MHC*(KB*gamma_h.dx(1) + KB*T + KB*T.dx(1) \
#                      + gamma_h*KB + CHI*mu.dx(1) \
#                      + CHI*mu + ECHARGE) \
#            + NC*dFermi_b(gamma_e)*MEC*(KB*T*gamma_e.dx(1) + gamma_e*KB*T.dx(1) \
#                                         + CHI*mu*mu.dx(1) - ECHARGE*phi.dx(1)) \
#            + nd_e*MEC*(KB*gamma_e.dx(1) + KB*T + KB*T.dx(1) \
#                       + gamma_e*KB + CHI*mu.dx(1) \
#                       + CHI*mu - ECHARGE)  # Test function of Jy

# Jy_testn = -v_5.dx(1) * NV*dFermi_b(gamma_hn)*MHC*(  ECHARGE) \
#            +v_5.dx(1) * NC*dFermi_b(gamma_en)*MEC*(- ECHARGE)
# 
# Jy_testn2 = -v_5.dx(1) * assemble(NV*dFermi_b(gamma_hn)*MHC*(  ECHARGE)  * ds(0))\
#            +v_5.dx(1) * assemble(NC*dFermi_b(gamma_en)*MEC*(- ECHARGE) * ds(0))

dUdt = dfb_deta(Constant(0.0), eta, mu)*(eta - eta_n)/dt + dfb_dmu(Constant(0.0), eta, mu)*(mu - mu_n)/dt

# Define free energy densities
# fb = AN1*(T - T1 - Tcvar)/(2*TC)*eta**2 + AN2/4*eta**4 + AN3/6*eta**6 \
#      + AU1*(T - T2 - Tcvar)/(2*TC)*mu**2 + AU2/4*mu**4 + AU3/6*mu**6 \
#      + GNU1*eta*mu - GNU2/2*eta**2*mu**2 + GNU3/2*eta**3*mu
# gamma_ein = (-CHI*mu**2/2 + CHP_IN)/(KB*T)   # Intrinsic electron equilibrium gamma
# gamma_hin = (-CHI*mu**2/2 - CHP_IN)/(KB*T)  # Intrinsic hole equilibrium gamma
# feh = CHI*mu**2/2*(NC*Fermi(gamma_e) + NV*Fermi(gamma_h)) \
#       + ECHARGE*phi*(NV*Fermi(gamma_h) - NC*Fermi(gamma_e)) \
#       + KB*T*(NC*Fermi(gamma_e)*gamma_e + NV*Fermi(gamma_h)*gamma_h) \
#       - KB*T*(NC*Fermi3half(gamma_e) + NV*Fermi3half(gamma_h)) \
#       + KB*T*(NC*Fermi3half(gamma_ein) + NV*Fermi3half(gamma_hin))
# ftot = fb + KAPPAN/2*dot(grad(eta), grad(eta)) + KAPPAU/2*dot(grad(mu), grad(mu)) + feh

#-------------------------------------------------------------------------
# Define variational problem
#-------------------------------------------------------------------------
qd = 6

# Lian times a factor 2
Feta = ((eta - eta_n)/dt + 2.0 * KN*dfb_deta(T, eta, mu))*v_1*dx(metadata={'quadrature_degree': qd}) \
       + (KN*KAPPAN)*dot(grad(eta), grad(v_1))*dx(metadata={'quadrature_degree': qd}) \
       - (KN*KAPPAN/deff)*(etas - eta)*v_1*ds(metadata={'quadrature_degree': qd})   # Boundary condition: interact with surroundings having an effective order parameter of etas
# detadt_0 = project(-KN*2*dfb_deta(T_0, eta_0, mu_0), V1)   # Calculate initial time derivatives of fields according to the evolution equations

# Lian times a factor 2
Fmu = ((mu - mu_n)/dt + 2.0 * KU*(dfb_dmu(T, eta, mu) \
       + CHI*mu*(nd_e + nd_h - 2*nd_in)))*v_2*dx(metadata={'quadrature_degree': qd}) + (KU*KAPPAU)*dot(grad(mu), grad(v_2))*dx(metadata={'quadrature_degree': qd}) \
      - (KU*KAPPAU/deff)*(mus - mu)*v_2*ds(metadata={'quadrature_degree': qd})  # Boundary condition: interact with surroundings having an effective order parameter of mus
# dmudt_0 = project(-KU*2*(dfb_dmu(T_0, eta_0, mu_0) + CHI*mu_0*(nd_e0 + nd_h0 - 2*nd_in0)), V1)


Fe = (NC*dFermi(gamma_e)*(gamma_e - gamma_en)/dt \
      - KEH0*mu**2*(nd_eeq*nd_heq - nd_e*nd_h))*v_3*dx(metadata={'quadrature_degree': qd}) - dot(j_e, grad(v_3))*dx(metadata={'quadrature_degree': qd})
# dgammaedt_0 = project(KEH0*mu_0**2*(nd_eeq0*nd_heq0 - nd_e0*nd_h0)/(NC*dFermi(gamma_e0)), V1)

Fh = (NV*dFermi(gamma_h)*(gamma_h - gamma_hn)/dt \
      - KEH0*mu**2*(nd_eeq*nd_heq - nd_e*nd_h))*v_4*dx(metadata={'quadrature_degree': qd}) - dot(j_h, grad(v_4))*dx(metadata={'quadrature_degree': qd})
# dgammahdt_0 = project(KEH0*mu_0**2*(nd_eeq0*nd_heq0 - nd_e0*nd_h0)/(NV*dFermi(gamma_h0)), V1)

Fphi = dot(grad(phi), grad(v_5))*dx(metadata={'quadrature_degree': qd}) - (ECHARGE/PERMITTIVITY)*(nd_h - nd_e)*v_5*dx(metadata={'quadrature_degree': qd})
# Fphi = PERMITTIVITY*dot(grad(phi), grad(v_5))*dx - ECHARGE*(nd_h - nd_e)*v_5*dx - PERMITTIVITY*delV/Constant(Ly)*v_5*ds(0)
# dphidt_0 = project(Constant(0.0), V1)

FT = (CPV*(T - T_n)/dt \
      - ECHARGE*((j_hx - j_ex)**2/(nd_e*MEA + nd_h*MHA) \
                 + (j_hy - j_ey)**2/(nd_e*MEC + nd_h*MHC)) \
      + dUdt \
      + (HTRAN/Lz)*(T - Ts))*v_6*dx(metadata={'quadrature_degree': qd}) \
     + THETA*dot(grad(T), grad(v_6))*dx(metadata={'quadrature_degree': qd})
# dUdt_0 = dfb_deta(Constant(0.0), eta_0, mu_0)*detadt_0 + dfb_dmu(Constant(0.0), eta_0, mu_0)*dmudt_0
# dTdt_0 = project((ECHARGE*((j_hx0 - j_ex0)**2/(nd_e0*MEA + nd_h0*MHA) + (j_hy0 - j_ey0)**2/(nd_e0*MEC + nd_h0*MHC)) - dUdt_0 - (HTRAN/Lz_do)*(T_0 - Ts))/CPV, V1)

# Define nonstandard boundary conditions using Lagrange multiplier
# Fbc_e = ((gamma_e*KB*T + CHI*mu**2/2 - CHP_IN)*v_7 \
#          + (v_3*KB*T + gamma_e*KB*v_6 + CHI*mu*v_2)*lam_1)*ds(0) \
#         + ((gamma_e*KB*T + CHI*mu**2/2 - CHP_IN)*v_7 \
#            + (v_3*KB*T + gamma_e*KB*v_6 + CHI*mu*v_2)*lam_1)*ds(1)


# Fbc_h = ((gamma_h*KB*T + CHI*mu**2/2 + CHP_IN)*v_8 \
#          + (v_4*KB*T + gamma_h*KB*v_6 + CHI*mu*v_2)*lam_2)*ds(0) \
#         + ((gamma_h*KB*T + CHI*mu**2/2 + CHP_IN)*v_8 \
#            + (v_4*KB*T + gamma_h*KB*v_6 + CHI*mu*v_2)*lam_2)*ds(1)


# Fbc_en = ((gamma_e*KB*T_n + CHI*mu_n**2/2 - CHP_IN)*v_7 \
#          + (v_3*KB*T_n)*lam_1)*ds(0) \
#         + ((gamma_e*KB*T_n + CHI*mu_n**2/2 - CHP_IN)*v_7 \
#            + (v_3*KB*T_n)*lam_1)*ds(1)


# Fbc_hn = ((gamma_h*KB*T_n + CHI*mu_n**2/2 + CHP_IN)*v_8 \
#          + (v_4*KB*T_n)*lam_2)*ds(0) \
#         + ((gamma_h*KB*T_n + CHI*mu_n**2/2 + CHP_IN)*v_8 \
#            + (v_4*KB*T_n)*lam_2)*ds(1)

# Fbc_phi_0 = (delV - phi - Resistor*Lz_do*assemble(Jy_n*ds(0)))*v_9*ds(0) \
#             - lam_3 *  v_5 * ds(0) \
#             - Resistor*Lz_do*Jy_testn*assemble(lam_3 * ds(0))*ds(0)

# Fbc_phi_0 = (delV - phi - integral_phi )*v_9*ds(0) \
#             - lam_3  *  v_5 * ds(0) \
#             + (integral_phi /60.0 - Resistor * Lz_do * Jy ) * v_10 *ds(0) \
            # - lam_3  *  v_10 * ds(0)\
            # - lam_3  *  v_10 * ds(0)\

# -------------------------------------------------------
# Nitsche's trick: Define nonstandard boundary conditions
# -------------------------------------------------------

# Yin corrected expression of Fbc_h_Nitche
epsilon = Constant(alpha.value*(Lx.value/nx.value))
Fbc_e_Nitsche = -1.0/epsilon*(gamma_e*KB*T +CHI*mu**2/2 - CHP_IN )*v_3*ds(0, metadata={'quadrature_degree': qd}) \
                -1.0/epsilon*(gamma_e*KB*T +CHI*mu**2/2 - CHP_IN )*v_3*ds(1, metadata={'quadrature_degree': qd})
Fbc_h_Nitsche = -1.0/epsilon*(gamma_h*KB*T +CHI*mu**2/2 + CHP_IN )*v_4*ds(0, metadata={'quadrature_degree': qd}) \
                -1.0/epsilon*(gamma_h*KB*T +CHI*mu**2/2 + CHP_IN )*v_4*ds(1, metadata={'quadrature_degree': qd})
Fbc_phi_Nitsche = -1.0/epsilon*(phi + (Resistor*Lz)*integral_phi - delVr \
                                + (Resistor*Capacitor)*(phi - phi_n)/dt)*v_5*ds(0, metadata={'quadrature_degree': qd}) \
                  + (integral_phi - Lx * Jy)*v_10*ds(0, metadata={'quadrature_degree': qd})
# dintegralphidt_0 = project(Constant(0.0), VR)

F = Feta + Fmu + Fe + Fh + Fphi + FT + Fbc_e_Nitsche + Fbc_h_Nitsche + Fbc_phi_Nitsche

# Gateaux derivative in direction of du (variational form for solving for Jacobian of F)
Jac = derivative(F, u, du)

#----------------------------------------------------------
# Define Block Gauss Seidel preconditioner weak form definition
#----------------------------------------------------------


#eta 
Feta_pc = ((eta - eta_n)/dt + 2.0 * KN*dfb_deta(T_n, eta, mu))*v_1*dx(metadata={'quadrature_degree': qd}) \
       + (KN*KAPPAN)*dot(grad(eta), grad(v_1))*dx(metadata={'quadrature_degree': qd}) \
       - (KN*KAPPAN/deff)*(etas - eta)*v_1*ds(metadata={'quadrature_degree': qd})   # Boundary condition: interact with surroundings having an effective order parameter of etas

#mu 
nd_e_pc = NC*Fermi_b(gamma_en)    # Electron density
nd_h_pc = NV*Fermi_b(gamma_hn) 
nd_in_pc = NC*Fermi_b((-CHI*mu**2/2 + CHP_IN)/(KB*T_n)) # contain only mu,T
Fmu_pc = ((mu - mu_n)/dt + 2.0 * KU*(dfb_dmu(T_n, eta, mu) \
       + CHI*mu*(nd_e_pc + nd_h_pc - 2*nd_in_pc)))*v_2*dx(metadata={'quadrature_degree': qd}) + (KU*KAPPAU)*dot(grad(mu), grad(v_2))*dx(metadata={'quadrature_degree': qd}) \
      - (KU*KAPPAU/deff)*(mus - mu)*v_2*ds(metadata={'quadrature_degree': qd})  # Boundary condition: interact with surroundings having an effective order parameter of mus

#gamma_e
nd_eeq_pc = NC*Fermi_b((-CHI*mu**2/2 + ECHARGE*phi + CHP_IN)/(KB*T_n))
nd_heq_pc = NV*Fermi_b((-CHI*mu**2/2 - ECHARGE*phi - CHP_IN)/(KB*T_n)) 
j_ex_pcpnp = -nd_e*MEA/ECHARGE*(KB*T_n*gamma_e.dx(0) + gamma_e*KB*T_n.dx(0) \
                          + CHI*mu*mu.dx(0) - ECHARGE*phi.dx(0))
j_ey_pcpnp = -nd_e*MEC/ECHARGE*(KB*T_n*gamma_e.dx(1) + gamma_e*KB*T_n.dx(1) \
                          + CHI*mu*mu.dx(1) - ECHARGE*phi.dx(1))
j_e_pcpnp = as_vector([j_ex_pcpnp, j_ey_pcpnp])  
Fe_pc = (NC*dFermi(gamma_e)*(gamma_e - gamma_en)/dt \
      - KEH0*mu**2*(nd_eeq_pc*nd_heq_pc - nd_e*nd_h))*v_3*dx(metadata={'quadrature_degree': qd}) - dot(j_e_pcpnp, grad(v_3))*dx(metadata={'quadrature_degree': qd}) \
        -1.0/epsilon*(gamma_e*KB*T_n +CHI*mu**2/2 - CHP_IN )*v_3*ds(0,metadata={'quadrature_degree': qd}) \
        -1.0/epsilon*(gamma_e*KB*T_n +CHI*mu**2/2 - CHP_IN )*v_3*ds(1,metadata={'quadrature_degree': qd})

#gamma_h
j_hx_pcpnp = -nd_h*MHA/ECHARGE*(KB*T_n*gamma_h.dx(0) + gamma_h*KB*T_n.dx(0) \
                          + CHI*mu*mu.dx(0) + ECHARGE*phi.dx(0))
j_hy_pcpnp = -nd_h*MHC/ECHARGE*(KB*T_n*gamma_h.dx(1) + gamma_h*KB*T_n.dx(1) \
                          + CHI*mu*mu.dx(1) + ECHARGE*phi.dx(1))
j_h_pcpnp = as_vector([j_hx_pcpnp, j_hy_pcpnp]) 
Fh_pc = (NV*dFermi(gamma_h)*(gamma_h - gamma_hn)/dt \
        - KEH0*mu**2*(nd_eeq_pc*nd_heq_pc - nd_e*nd_h))*v_4*dx(metadata={'quadrature_degree': qd}) - dot(j_h_pcpnp, grad(v_4))*dx(metadata={'quadrature_degree': qd}) \
        - 1.0/epsilon*(gamma_h*KB*T_n +CHI*mu**2/2 + CHP_IN )*v_4*ds(0, metadata={'quadrature_degree': qd}) \
        - 1.0/epsilon*(gamma_h*KB*T_n +CHI*mu**2/2 + CHP_IN )*v_4*ds(1, metadata={'quadrature_degree': qd})

#phi
Jy_pcpnp =ECHARGE*(j_hy_pcpnp - j_ey_pcpnp) 
Fbc_phi_Nitsche_pc = -1.0/epsilon*(phi + (Resistor*Lz)*integral_phi - delVr \
                                + (Resistor*Capacitor)*(phi - phi_n)/dt)*v_5*ds(0, metadata={'quadrature_degree': qd}) \
                  + (integral_phi - Lx * Jy_pcpnp)*v_10*ds(0, metadata={'quadrature_degree': qd})

Fphi_pc = dot(grad(phi), grad(v_5))*dx(metadata={'quadrature_degree': qd}) - (ECHARGE/PERMITTIVITY)*(nd_h - nd_e)*v_5*dx(metadata={'quadrature_degree': qd}) \
        + Fbc_phi_Nitsche_pc 
 
 # T 
FT_pc = FT

F_pc = Feta_pc + Fmu_pc + Fe_pc + Fh_pc + Fphi_pc + FT_pc 

Jac_pc = derivative(F_pc,u,du) 

#----------------------------------------------------------
# Solve the problem and save the solution
#----------------------------------------------------------
# Encapsulate the process of updating the previous solution "u_n"
# def update_u_n(u):
#     assign(u_n.sub(0), u.sub(0))
#     assign(u_n.sub(1), u.sub(1))
#     assign(u_n.sub(2), u.sub(2))
#     assign(u_n.sub(3), u.sub(3))
#     assign(u_n.sub(4), u.sub(4))
#     assign(u_n.sub(5), u.sub(5))

# Encapsulate the solving process for the solution at "t" (u_n not updated)
def solve4dt(t):
    delVr.update(ramped_delV(t))

    # problem = NonlinearVariationalProblem(F, u, bcs, Jac)
    # solver = NonlinearVariationalSolver(problem)

    if use_GS_block_preconditioner: 
        num, converged = solver.solve(problem, u.vector())    
    else:
        num, converged = solver.solve()   # Return whether the solver converges

    # Update previous solution
    # update_u_n(u)
    return converged

# Create VTK files for visualization output
vtkfile_eta = File(fig_file+'eta.pvd')
vtkfile_mu = File(fig_file+'mu.pvd')
vtkfile_nd_e = File(fig_file+'n.pvd')
vtkfile_nd_h = File(fig_file+'p.pvd')
vtkfile_phi = File(fig_file+'phi.pvd')
vtkfile_T = File(fig_file+'T.pvd')
# vtkfile_gamma_e = File(fig_file+'gamma_e.pvd')
# vtkfile_gamma_h = File(fig_file+'gamma_h.pvd')
# Encapsulate the exportation of solutions to files
def save_sol(u, t):
    _u = u.split()   # Not a deep copy, which is fast and memory saving
 
    # Save solution to file (VTK)
    _u[0].rename("eta", "eta")
    _u[1].rename("mu", "mu")
    vtkfile_eta << (_u[0], t)
    vtkfile_mu << (_u[1], t)
    # _u[2].rename("gamma_e", "gamma_e")
    # vtkfile_gamma_e << (_u[2], t)
    _nd_e = project(NC*Fermi(_u[2])*UCVOL, V1)
    _nd_e.rename("n","n")
    vtkfile_nd_e << (_nd_e, t)
    # _u[3].rename("gamma_h", "gamma_h")
    # vtkfile_gamma_h << (_u[3], t)
    _nd_h = project(NV*Fermi(_u[3])*UCVOL, V1)
    _nd_h.rename("p","p")
    vtkfile_nd_h << (_nd_h, t)
    _phi = project(_u[4]*Constant(VUNIT), V1)
    _phi.rename("phi","phi")
    vtkfile_phi << (_phi, t)
    _T = project(_u[5]*Constant(TEMPUNIT), V1)
    _T.rename("T","T")
    vtkfile_T << (_T, t)

# # hdf5 version of saving solution
# hdf5_eta = HDF5File(comm, fig_file+'eta.hdf5', 'w')
# hdf5_mu = HDF5File(comm, fig_file+'mu.hdf5', 'w')
# hdf5_nd_e = HDF5File(comm, fig_file+'n.hdf5', 'w')
# hdf5_nd_h = HDF5File(comm, fig_file+'p.hdf5', 'w')
# hdf5_phi = HDF5File(comm, fig_file+'phi.hdf5', 'w')
# hdf5_T = HDF5File(comm, fig_file+'T.hdf5', 'w')
# # Encapsulate the exportation of solutions to files
# def save_sol_hdf5(u, t):
#     _u = u.split()   # Not a deep copy, which is fast and memory saving
#  
#     # Save solution to file (VTK)
#     _u[0].rename("eta", "eta")
#     hdf5_eta.write(_u[0], 'eta', t)
# 
#     _u[1].rename("mu", "mu")
#     hdf5_mu.write(_u[1], 'mu', t)
# 
#     _nd_e = project(NC*Fermi(_u[2])*UCVOL, V1)
#     _nd_e.rename("n","n")
#     hdf5_nd_e.write(_nd_e, 'n', t)
# 
#     _nd_h = project(NV*Fermi(_u[3])*UCVOL, V1)
#     _nd_h.rename("p","p")
#     hdf5_nd_h.write(_nd_h, 'p', t)
# 
#     _phi = project(_u[4]*Constant(VUNIT), V1)
#     _phi.rename("phi","phi")
#     hdf5_phi.write(_phi, 'phi', t)
# 
#     _T = project(_u[5]*Constant(TEMPUNIT), V1)
#     _T.rename("T","T")
#     hdf5_T.write(_T, 'T', t)

# Encapsulate the calculation of the relative error for adaptive time-stepping
# def rel_err(u1, u2):
#     # rerr = 0.0
#     # for i in range(n_f):
#     #     _u1_v = _u1[i].vector().get_local()
#     #     _u2_v = _u2[i].vector().get_local()
# 
#     #     rerr_i = np.nanmax( np.abs(_u1_v - _u2_v) / (np.abs(_u2_v) + 1E-15) )
#     #     rerr_i_global = MPI.max(comm, rerr_i)
# 
#     #     rerr += rerr_i_global**2
#     # rerr = math.sqrt(rerr)  
#     u1_v = u1.vector().get_local()
#     u2_v = u2.vector().get_local()
# 
#     rerr = np.nanmax( np.abs(u1_v - u2_v) / (np.abs(u2_v) + 1E-15) )
#     #rerr_global = MPI.max(comm, rerr)
# 
#     return rerr
#     #return rerr_global


def rel_err_L2(u1, u2): 
    # Calculate relative error for each field
    #rerror1 = errornorm(u1.sub(0), u2.sub(0), 'L2') / (errornorm(Constant(0.0), u1.sub(0), 'L2') + 1e-15)
    #rerror2 = errornorm(u1.sub(1), u2.sub(1), 'L2') / (errornorm(Constant(0.0), u1.sub(1), 'L2') + 1e-15)
    #rerror3 = errornorm(u1.sub(2), u2.sub(2), 'L2') / (errornorm(Constant(0.0), u1.sub(2), 'L2') + 1e-15)
    #rerror4 = errornorm(u1.sub(3), u2.sub(3), 'L2') / (errornorm(Constant(0.0), u1.sub(3), 'L2') + 1e-15)
    #rerror5 = errornorm(u1.sub(4), u2.sub(4), 'L2') / (errornorm(Constant(0.0), u1.sub(4), 'L2') + 1e-15)
    #rerror6 = errornorm(u1.sub(5), u2.sub(5), 'L2') / (errornorm(Constant(0.0), u1.sub(5), 'L2') + 1e-15)

    delu = Function(V)
    delu.vector().set_local(u2.vector().get_local() - u1.vector().get_local())  # Difference of u1 and u2
    _u1 = u1.split()   # Not a deep copy, which is fast and memory saving
    _delu = delu.split()
    n = len(_u1)

    r = 0.0
    for i in range(n):
        r += assemble(_delu[i] ** 2 * dx) / (assemble(_u1[i] ** 2 * dx) + 1e-15)

    return math.sqrt(r / n)
    
# Define reservoir functions for adaptive time-stepping
# dudt_n = Function(V)  # Stores o(dt) approximation of time direvative of u at last time step
# fa.assign(dudt_n, [detadt_0, dmudt_0, dgammaedt_0, dgammahdt_0, dphidt_0, dTdt_0, dintegralphidt_0])
# u_guess = Function(V) 
# u_guess.assign(u)

# Create progress bar
# progress = Progress('Time-stepping')
if fenicslog.value == 'DBG':
    set_log_level(LogLevel.DBG)
elif fenicslog.value == 'TRACE':
    set_log_level(LogLevel.TRACE)
elif fenicslog.value == 'PROGRESS':
    set_log_level(LogLevel.PROGRESS)
elif fenicslog.value == 'INFO':
    set_log_level(LogLevel.INFO)
elif fenicslog.value == 'WARNING':
    set_log_level(LogLevel.WARNING)
elif fenicslog.value == 'ERROR':
    set_log_level(LogLevel.ERROR)
elif fenicslog.value == 'CRITICAL':
    set_log_level(LogLevel.CRITICAL)
elif fenicslog.value == 'FALSE':
    set_log_active(False)
else:
    set_log_level(LogLevel.INFO)


#------------------------------------------------------------------
# Time-stepping
#------------------------------------------------------------------
# t_out = dt_out
t = 0.0
n_step = 0
n_out = 0  # Mark the number of the output time
Tfail = 0
Nfail = 0
otherfail = 0
successive_div = 0
# save_sol(u) Xu:deleted

# dtt = dt 
# dt = 1E-10

logfile = open(allin1file, 'w')
if use_GS_block_preconditioner.value: 
    print("Using Gauss-Seidel block preconditioner")
    problem = Problem(Jac,Jac_pc,F, bcs) #user defined problem 
    solver = CustomSolver(mesh)
    prm = solver.parameters
    prm['relaxation_parameter'] = rel_par.value 
    prm['maximum_iterations'] = max_iter.value 
    prm['absolute_tolerance'] = newton_abs_tol.value
    prm['relative_tolerance'] = newton_rel_tol.value
    prm['error_on_nonconvergence'] = False

else: 
    problem = NonlinearVariationalProblem(F, u, bcs, Jac)
    solver = NonlinearVariationalSolver(problem)

    prm = solver.parameters
    prm['newton_solver']['relaxation_parameter'] = rel_par.value
    prm['newton_solver']['maximum_iterations'] = max_iter.value
    prm['newton_solver']['absolute_tolerance'] = newton_abs_tol.value
    prm['newton_solver']['relative_tolerance'] = newton_rel_tol.value
    prm['newton_solver']['error_on_nonconvergence'] = False   # Not stop program if Newton solver does not converge

    prm['newton_solver']['linear_solver'] = directsolver.value   # 'superlu_dist', 'mumps', 'gmres', 'bicgstab'
    # prm['newton_solver']['krylov_solver']['monitor_convergence'] = True
    # prm['newton_solver']['preconditioner'] = 'hypre_amg' #'hypre_euclid'   # 'hypre_euclid'
    # prm["newton_solver"]["krylov_solver"]["maximum_iterations"] = kr_solver_iter
    # prm["newton_solver"]["krylov_solver"]['error_on_nonconvergence'] = False  # Not stop program if linear solver does not converge

# Write log file header. Tfail is the accumulative number of time refinement, Nfail is the accumulative number of Newton solver nonconvergence, and Terr is the final L2 error of time stepping.
if rank == 0:
    logfile.write(f'          #Step            Time       Time step           Tfail           Nfail      Other fail    Av. EOP norm       Av. T (K)           V (V)         R (Ohm)\n')
    logfile.flush()

while t < tf.value + 1e-9*tf.value:
    if rank == 0: 
        print('------------------------------------------------------------------------', flush=True)
        print(f'User message ===> Start time step refinement at t = {t:15.6e} with dt = {dt.value:15.6e} (tol = {tol_tstep.value:15.6f})', flush=True)

    rerror = 1000.0 # set big to initiate time step refinement
    dt_factor = 0.999  # set to smaller than 1 to initiate time step refinement

    while rerror > tol_tstep.value:
        if successive_div > max_div.value:
            break

        # Solve for the current time step
        try:
            converged = solve4dt(t+dt.value) # u_n not updated
        except:
            otherfail += 1   # Count the number of other fails (probably linear solver fail)
            successive_div += 1
            u.assign(u_n)  # Restore initial guess for newton iteration 
            dt.update(dt.value * r_t_min.value)
            if rank == 0:
                print(f'User message ===> Solver did not converge for other reason! Refined time step dt = {dt.value/r_t_min.value:15.6e} --> {dt.value:15.6e} and test again', flush=True)
            continue

        if not converged: 
            Nfail += 1   # Count the number of Newton solver nonconvergence
            successive_div += 1
            u.assign(u_n)  # Restore initial guess for newton iteration 
            dt.update(dt.value * r_t_min.value)
            if rank == 0:
                print(f'User message ===> Newton solver did not converge! Refined time step dt = {dt.value/r_t_min.value:15.6e} --> {dt.value:15.6e} and test again', flush=True)
            continue

        successive_div = 0  # set back to zero because solve is successful

        # Compute the second-order accurate estimate of the solution; the current solution u is first-order accurate (backward Euler)
        if n_step > 0:
            u2acc = project(u_n + (dudt_n + (u - u_n)/dt)*dt/2.0, V)
        else:
            break

        # Calculate the backward Euler time integration error and time step changing factor
        rerror = rel_err_L2(u2acc, u)
        dt_factor = min(max(s_t.value*math.sqrt(tol_tstep.value/max(rerror, 1e-10)), r_t_min.value), r_t_max.value)
        
        # Adjust things if the time stepping accuracy is not met
        if rerror > tol_tstep.value:
            Tfail += 1  # Count the number of time refinement
            u.assign(u_n) # Restore initial guess for newton iteration 
            dt.update(dt.value * dt_factor)
            if rank == 0: 
                print(f'User message ===> Time stepping error {rerror:15.6e} is too big, refined time step dt = {dt.value/dt_factor:15.6e} --> {dt.value:15.6e} and test again', flush=True)

    if successive_div > max_div.value:
        break
    
    # Time step refinement is done, now save solution, update solution, and output
    n_step += 1
    if savemethod.value == 'auto':
        if n_step % saveperiod.value == 1:
            save_sol(u, t+dt.value)
    # Output at specified times by linear interpolation
    elif savemethod.value == 'fix':
        while t_out[n_out] <= t + dt.value:
            u_out = project(u_n + (u-u_n)*Constant(t_out[n_out]-t)/dt, V)  # Compute linear interpolation between current solution and previous solution
            save_sol(u_out, t_out[n_out])
            if rank == 0:
                print(f'User message ===> Out: {n_out:15d}, t = {t_out[n_out]:15.6e}', flush=True)
                logfile.write(f' - {n_out:12d} {t_out[n_out]:15.6e} - out\n')
                logfile.flush()
            n_out += 1
            if n_out > len_t_out - 1:
                break
            
    mu_norm_av = math.sqrt(assemble(mu**2*dx)/(Lx.value*Ly.value))
    V_VO2 = assemble(phi*ds(0)) / Lx.value * VUNIT  # Calculate voltage drop across VO2
    Tav = assemble(T*dx) / (Lx.value*Ly.value) * TEMPUNIT  # Calculate average temperature across VO2
    R_VO2 = V_VO2 / (Lz.value*assemble(Jy*ds(0)) * CUNIT/TUNIT)  # Calculate VO2 resistance
    if rank == 0:
        print(f'User message ===> Completed refinement: dt = {dt.value:15.6e}, t = {t:15.6e} --> {t+dt.value:15.6e}', flush=True)
        logfile.write(f'{n_step:15d} {t+dt.value:15.6e} {dt.value:15.6e} {Tfail:15d} {Nfail:15d} {otherfail:15d} {mu_norm_av:15.6f} {Tav:15.6e} {V_VO2:15.6e} {R_VO2:15.6e}\n') 
        logfile.flush()

    dudt_n = project((u - u_n)/dt, V)
    u_n.assign(u)
    t += dt.value 

    if dt_factor > 1: 
        dt.update(dt.value * dt_factor)
        if rank == 0:
            print(f'User message ===> dt is increased: dt = {dt.value/dt_factor:15.6e} --> {dt.value:15.6e}', flush=True)
        


comm.Barrier()
elapsed_time = time.time() - start_time
if successive_div > max_div.value:
    if rank == 0:
        print(f'User message ===> !!! Solving process diverged, stop the process !!! Computation time: {elapsed_time:14.3f} s.', flush=True)
        logfile.write(f'!!! Solving process diverged, stop the process !!! Computation time: {elapsed_time:14.3f} s.')
else:
    if rank == 0:
        print(f'User message ===> Finished computation, computation time: {elapsed_time:14.3f} s.', flush=True)
        logfile.write(f'Finished computation, computation time: {elapsed_time:14.3f} s.')

logfile.close()


# Hold plot
# interactive()
