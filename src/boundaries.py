'''Boundary classes for alveolar perfusion and gas exchange simulations.
'''

__author__ = 'pzuritas'
__email__ = 'pzurita@uc.cl'

from dolfin import *


class GammaIn(SubDomain):
    '''Subdomain class for boundary conditions.'''
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
        return near(x[0], self.dir_min, self.tol)

class GammaOut(SubDomain):
    '''Subdomain class for boundary conditions.'''
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
        return near(x[0], self.dir_max, self.tol)

class GammaAir(SubDomain):
    '''Subdomain class for boundary conditions.'''
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

class GammaSlabPi(SubDomain):
    '''Subdomain class for periodic boundary conditions in slab mesh.'''
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

class GammaAirSlabPi(GammaAir):
    '''Alternative subdomain class for \Gamma_{\text{air}} when using
    periodic boundary conditions on slab mesh.
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

class GammaTKDPi(SubDomain):
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

class GammaAirTKD(SubDomain):
    '''Subdomain class for periodic boundary conditions in slab mesh.'''
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
            near(x[0], self.dir_max, self.tol) \
                or near(x[0], self.dir_max, self.tol)
        ) and not (
            near(x[1], self.dir_min, self.tol) \
                or near(x[1], self.dir_max, self.tol)
        ) and not (
            near(x[2], self.dir_min, self.tol) \
                or near(x[2], self.dir_max, self.tol)
        )