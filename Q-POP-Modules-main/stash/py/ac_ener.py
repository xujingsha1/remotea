import sympy as sp
from fenics import *

def drivingForce(u, gamma1, gamma2):
    return gamma2 / gamma1 * (u**3 - u)

def init_condition(u):
    x = sp.symbols(['x[{}]'.format(i) for i in range(u.geometric_dimension())])
    funcInit = x[0]**2 * sp.Sin(2 * sp.pi * x[0])
    return sp.ccode(funcInit)
