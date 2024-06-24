import time, json, itertools
from fenics import *
import numpy as np
# import mpi4py
from ufl import nabla_div
from mshr import *
from terminalPrinting import terminalPrint as tP
# import terminalPrinting.terminalPrint as tP

class voigtTransfers():

    def voigt3stress(vec):
        return as_tensor(
            [[vec[0], vec[5], vec[4]], [vec[5], vec[1], vec[3]], [vec[4], vec[3], vec[2]]]
        )

    def strain3voigt(ten):
        # FEniCS does not know anything about Voigt notation,
        # so one need to access the components directly as eps[0, 0] etc.
        return as_vector(
            [ten[0, 0], ten[1, 1], ten[2, 2], 2 * ten[1, 2], 2 * ten[0, 2], 2 * ten[0, 1]]
        )

    def voigt3strain(vec):
        return as_tensor(
        [[vec[0], vec[5]/2, vec[4]/2],
        [vec[5]/2, vec[1], vec[3]/2],
        [vec[4]/2, vec[3]/2, vec[2]]
        ]
        )

    def stress3voigt(ten):
        return as_vector(
        [ten[0,0], ten[1,1], ten[2,2], ten[1,2], ten[0,2], ten[0,1]]
        )

    def strain2voigt(ten):
        return as_vector(
            [ten[0,0], ten[1,1], 2*ten[0,1]]
        )

    def voigt2stress(vec):
        return as_tensor(
            [[vec[0], vec[2]], [vec[2], vec[1]]]
        )

class stressStrainCalculators():

    def epsilon(u):
        return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    
    def sigma(u, C=None, lame=(1.25, 1)):
        if C != None:
            stressStrainCalculators.sigma_calc(C, u)
        else: #defined via lame constants
            return lame[0] * nabla_div(u)*Identity(u.geometric_dimension()) + 2 * lame[1] * stressStrainCalculators.epsilon(u)

    def sigma_calc(C, u):
        if u.geometric_dimension() == 2:
            stressStrainCalculators.sigma2D(C,u)
        elif u.geometric_dimension() == 3:
            stressStrainCalculators.sigma3D(C,u)
        else:
            msg = 'The Geometric Dimension of Displacement Vector is not recognized. Please fix'
            tP.printErrorLine(msg)

    def sigma3D(C, u):
        return voigtTransfers.voigt3stress(dot(C, voigtTransfers.strain3voigt(stressStrainCalculators.epsilon(u))))

    def sigma2D(C, u):
        return voigtTransfers.voigt2stress(
            dot(C, voigtTransfers.strain2voigt(stressStrainCalculators.epsilon(u)))
        )

class tensorGenerator():

    #4th rank tensors

    # Cubic has 3 Independent Coefficients
    # Hexagonal has 5 Independent Coefficients
    # Trigonal (32,3M,-3M) has 6 Independent Coefficients
    # Trigonal (3,-3) has has 7 Independent Coefficients
    # Tetragonal (4MM, -42M, 422, 4/MMM) has 6 Independent Coefficients (C11,C12,C13,C33,C44,C66)
    # Tetragonal (4,-4,4/M) has 7 Independent Coefficients (C11,C12,C13,C16,C44,C66)
    # Orthorhombic (222,MM2,MMM) has 9 Independent Coefficients (C11,C12,C13,C22,C23,C33,C44,C55,C66)
    # Monoclinic
    # Triclinic

    #This creates tensor representations in Matrix form. A separate module will be issued for converting from Matrix
    #form into traditional tensor form

    def parse4thTensorCoefficients(tensorDict):

        #Get the number of independent coefficients in the dictionary that describes the properties of the material
        indCoef = len(list(tensorDict.keys()))

        #leader is the label of the coefficient, kept general here so we can parse all kinds of properties
        leader = list(tensorDict.keys())[0][0]

        if indCoef == 3:
            tP.printNamedLine('Setting up 4th Rank Cubic Tensor')
            tensor4th = rank4_cubic(tensorDict['{}{}'.format(leader,11)], tensorDict['{}{}'.format(leader,12)], tensorDict['{}{}'.format(leader,44)])
        elif indCoef == 5:
            errMsg = 'Hexagonal 4th Rank Tensors have not been implemented yet'
            tP.printErrorLine(errMsg)
        elif indCoef == 6:
            errMsg = '4th Rank Tensors for Trigonal and Tetragonal Crystals with 6 Independent Coffiecients is not yet implemented'
            tP.printErrorLine(errMsg)
        elif indCoef == 7:
            errMsg = '4th Rank Tensors for Trigonal/Tetragonal Crystals with 7 Independent Coefficients is not yet implemented'
            tP.printErrorLine(errMsg)
        elif indCoef == 9:
            errMsg = '4th Rank Tensors for Orthorhombic Crystals have not yet been implemented'
            tP.printErrorLine(errMsg)

        return tensor4th


    def rank4_cubic(a11, a12, a44):

        c_11 = c_22 = c_33 = a11
        c_44 = c_55 = c_66 = a44
        c_12 = c_13 = c_23 = a12
        c_14 = c_15 = c_16 = c_24 = c_25 = c_26 = c_34 = c_35 = c_36 =  0
        c_45 = c_46 = c_56 = 0
        C = as_tensor(
            [
                [c_11, c_12, c_13, c_14, c_15, c_16],
                [c_12, c_22, c_23, c_24, c_25, c_26],
                [c_13, c_23, c_33, c_34, c_35, c_36],
                [c_14, c_24, c_34, c_44, c_45, c_46],
                [c_15, c_25, c_35, c_45, c_55, c_56],
                [c_16, c_26, c_36, c_46, c_56, c_66],
            ]
        )

        return C

    def rank2_fromVec(a1, a2, a3, a4, a5, a6):

        return as_tensor([[a1,a6,a5],[a6,a2,a4],[a5,a4,a3]])
