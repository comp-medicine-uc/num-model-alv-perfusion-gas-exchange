'''Boundary classes for alveolar perfusion and gas exchange simulations.
'''

__author__ = 'pzuritas'
__email__ = 'pzurita@uc.cl'

from dolfin import *
import numpy as np


class GammaIn(SubDomain):
    '''Subdomain class for boundary conditions on the inlet of a rectangular
    prism.
    '''
    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of flow direction. (float)
        dir_max: maximum value of flow direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.dir_min = dir_min
        self.dir_max = dir_max
        self.tol = tol
    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and near(x[0], self.dir_min, self.tol)

class GammaOut(SubDomain):
    '''Subdomain class for boundary conditions on the outlet of a rectangular
    prism.
    '''
    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of flow direction. (float)
        dir_max: maximum value of flow direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.dir_min = dir_min
        self.dir_max = dir_max
        self.tol = tol
    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and near(x[0], self.dir_max, self.tol)

class GammaAir(SubDomain):
    '''Subdomain class for boundary conditions on the blood-air barrier of a
    rectangular prism.
    '''
    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of flow direction. (float)
        dir_max: maximum value of flow direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.dir_min = dir_min
        self.dir_max = dir_max
        self.tol = tol

    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and not (
            near(x[0], self.dir_min, self.tol/2) \
                or near(x[0], self.dir_max, self.tol/2)
        )

class GammaSheetPi(SubDomain):
    '''Subdomain class for periodic boundary conditions in sheet mesh.
    '''
    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of periodic direction (z). (float)
        dir_max: maximum value of periodic direction (z). (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.dir_min = dir_min
        self.dir_max = dir_max
        self.tol = tol

    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and near(x[2], self.dir_max, self.tol)

    def map(self, x, y):
        '''Maps opposite faces for periodic boundary conditions.
        
        x: position in subdomain.
        y: position in opposite subdomain.
        '''
        y[0] = x[0]
        y[1] = x[1]
        y[2] = x[2] + 6.0  # (self.dir_max - self.dir_min)

class GammaAirSheetPi(GammaAir):
    '''Alternative subdomain class for blood-air barrier when using
    periodic boundary conditions on sheet mesh.
    '''
    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of height direction. (float)
        dir_max: maximum value of height direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__(dir_min, dir_max, tol)

    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and (
            near(x[1], self.dir_max, self.tol) \
                or near(x[1], self.dir_min, self.tol)
        )

class GammaInSphereV2(SubDomain):
    '''Subdomain class for inlet boundary conditions in hollow sphere model.'''
    def __init__(self, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of flow direction. (float)
        dir_max: maximum value of flow direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.tol = tol
    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and near(
            -(106/20)*np.sqrt(x[1]**2 + x[2]**2), x[0],
            self.tol*abs(20/106*x[0])
        ) and x[0] < 0 and not (
            np.sqrt(x[1]**2 + x[2]**2) > 20.7
        )

class GammaOutSphereV2(SubDomain):
    '''Subdomain class for outlet boundary conditions in hollow sphere model.'''
    def __init__(self, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of flow direction. (float)
        dir_max: maximum value of flow direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.tol = tol
    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and near(
            (106/20)*np.sqrt(x[1]**2 + x[2]**2), x[0],
            self.tol*abs(20/106*x[0])
        ) and x[0] > 0 and not (
            np.sqrt(x[1]**2 + x[2]**2) > 20.7
        )

class GammaAirSphereV2(SubDomain):
    '''Subdomain class for blood-air barrier boundary conditions in hollow
    sphere model.'''
    def __init__(self, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of flow direction. (float)
        dir_max: maximum value of flow direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.tol = tol

    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and not (near(
            -(106/20)*np.sqrt(x[1]**2 + x[2]**2), x[0],
            self.tol*abs(20/106*x[0])
        ) and x[0] < 0) and not (near(
            (106/20)*np.sqrt(x[1]**2 + x[2]**2), x[0], self.tol*abs(20/106*x[0])
        ) and x[0] > 0) or (
            np.sqrt(x[1]**2 + x[2]**2) > 20.15
        )

class GammaTKDPiV2(SubDomain):
    '''Subdomain class for periodic boundary conditions in TKD mesh.'''
    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of periodic direction (z or y). (float)
        dir_max: maximum value of periodic direction (z or y). (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.dir_min = dir_min
        self.dir_max = dir_max
        self.tol = tol

    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and (
            near(x[1], self.dir_max, self.tol) \
                or near(x[2], self.dir_max, self.tol)
        )

    def map(self, x, y):
        '''Maps opposite faces for periodic boundary conditions.
        
        x: position in subdomain.
        y: position in opposite subdomain.
        '''
        y[0] = x[0]
        if near(x[1], self.dir_max, self.tol):
            y[1] = x[1] - (self.dir_max - self.dir_min)
            y[2] = x[2]
        else:
            y[1] = x[1]
            y[2] = x[2] - (self.dir_max - self.dir_min)

class GammaAirTKDV2(SubDomain):
    '''Subdomain class for blood-air barrier when using
    periodic boundary conditions on TKD mesh.
    '''
    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of periodic direction (z or y). (float)
        dir_max: maximum value of periodic direction (z or y). (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.dir_min = dir_min
        self.dir_max = dir_max
        self.tol = tol

    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and not (
            x[0] < 0 and near(
                x[1]**2 + x[2]**2, (100/5/sqrt(2))**2, 1E1
            )
        ) and not (
            x[0] > 0 and near(
                x[1]**2 + x[2]**2, (100/5/sqrt(2))**2, 1E1
            )
        ) and not (
            near(x[1], self.dir_min, self.tol) \
                or near(x[1], self.dir_max, self.tol/10)
        ) and not (
            near(x[2], self.dir_min, self.tol) \
                or near(x[2], self.dir_max, self.tol/10)
        )

class GammaTKDInV2(SubDomain):
    '''Subdomain class for inlet when using periodic boundary
    conditions on TKD mesh.
    '''
    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of flow direction. (float)
        dir_max: maximum value of flow direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.dir_min = dir_min
        self.dir_max = dir_max
        self.tol = tol
    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and x[0] < 0 and near(
           x[1]**2 + x[2]**2, (100/5/sqrt(2))**2, self.tol
        )

class GammaTKDOutV2(SubDomain):
    '''Alternative subdomain class for outlet when using periodic boundary
    conditions on TKD mesh.
    '''
    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of flow direction. (float)
        dir_max: maximum value of flow direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.dir_min = dir_min
        self.dir_max = dir_max
        self.tol = tol
    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        return on_boundary and x[0] > 0 and near(
           x[1]**2 + x[2]**2, (100/5/sqrt(2))**2, self.tol
        )