'''
Model file for alveolar perfusion and gas exchange simulations.

@author pzuritas
'''
from fenics import *
from dolfin import *
from mshr import Cylinder, generate_mesh
from datetime import date
import time as tm
import numpy as np
import os


class PerfusionGasExchangeModel():
    '''FEniCS simulater class for microscale alveolar perfusion and gas exchange
    model.
    '''
    def __init__(self, folder_path, params):
        '''Instance the model.

        folder_path: path to folder for storing solution files. (string)
        params: dictionary for parameter values. (dict)
        '''
        self.folder_path = folder_path
        self.params = params

        # Save details
        if not os.path.exists(self.folder_path):
            os.mkdir(self.folder_path)
        with open(self.folder_path+'/info.txt', 'w') as file:
            file.write(f'Simulation done on {date.today()}\n\n')
            file.write(f'Parameters used:\n\n')
            for param in self.params:
                file.write(f'Parameter {param}: {self.params[param]}\n')

    def import_h5mesh(self, mesh_path):
        '''Imports mesh from .h5 file for use in simulations.

        mesh_path: path to .h5 file. (string)
        '''
        hdf5 = HDF5File(MPI.comm_world, mesh_path, 'r')
        self.mesh = Mesh()
        hdf5.read(self.mesh, 'mesh', False)

        dir_arr = np.array(
                [coords[0] for coords in self.mesh.coordinates()]
            )  # Node coordinates for principal (x) direction
        self.dir_max = np.max(dir_arr)  # Maximum principal direction coordinate
        self.dir_min = np.min(dir_arr)  # Minimum principal direction coordinate

    def generate_slab_mesh(self, dims, elems, save=True):
        '''Generates a rectangular prism mesh for simulations on a slab.

        dims: dimensions of the mesh. (tuple)
        elems: number of elements in each direction. (tuple)
        save: save the generated mesh on a .pvd file. (bool)
        '''
        geometry = BoxMesh(
            Point(0, 0, 0), Point(*dims), *elems
        )
        self.mesh = geometry
        mesh_file = File(self.folder_path+'/mesh.pvd')
        mesh_file << self.mesh
        
        # Assign values to flow direction range

        self.dir_min = 0
        self.dir_max = dims[0]

    def generate_cylinder_mesh(self, end, r, save=True):
        '''Generates a cylindrical mesh for simulations on a tube.

        end: endpoint. (tuple)
        r: radius. (tuple)
        save: save the generated mesh on a .pvd file. (bool)
        '''
        geometry = Cylinder(
            Point(0, 0, 0), Point(*end), r, r, 200
        )
        self.mesh = generate_mesh(geometry, 50)  # What is 25?
        mesh_file = File(self.folder_path+'/mesh.pvd')
        mesh_file << self.mesh
        
        # Assign values to flow direction range

        self.dir_min = 0
        self.dir_max = end[0]

    def sim_p(self):
        '''Solves the perfusion (P) problem of the model.'''
        # Instance the funciton spaces por scalar and vector functions.

        self.W_h = FunctionSpace(self.mesh, 'Lagrange', 1)
        self.V_h = VectorFunctionSpace(self.mesh, 'Lagrange', 1)

        # Instance the relevant boundaries

        self.gamma_in = GammaIn(
            self.dir_min, self.dir_max, DOLFIN_EPS
        )
        self.gamma_out = GammaOut(
            self.dir_min, self.dir_max, DOLFIN_EPS
        )
        self.gamma_air = GammaAir(
            self.dir_min, self.dir_max, DOLFIN_EPS
        )

        # Declare the boundaries in the mesh and tag them

        self.boundaries = MeshFunction('size_t', self.mesh, dim=2)
        self.boundaries.set_all(3)
        self.gamma_in.mark(self.boundaries, 1)
        self.gamma_out.mark(self.boundaries, 2)
        self.gamma_air.mark(self.boundaries, 3)

        # Declare Dirichlet boundary conditions for (P)

        self.p_dbc = [
            DirichletBC(self.W_h, self.params['p_max'], self.gamma_in),
            DirichletBC(self.W_h, self.params['p_min'], self.gamma_out),
        ]

        # Assemble problem

        p = TrialFunction(self.W_h)
        v = TestFunction(self.W_h)
        f = Constant(0)
        a = inner(grad(p), grad(v))*dx
        F = f*v*dx

        # Solve problem

        self.p = Function(self.W_h)
        solve(
            a == F, self.p, self.p_dbc,
            solver_parameters={
                'linear_solver': 'gmres',
                'preconditioner': 'ilu'
            }
        )

        self.u = project(
            -1/self.params['mu']*self.params['kappa']*grad(self.p),
            self.V_h
        )

        # Save solutions

        u_file = File(self.folder_path+'/p/u.pvd')
        u_file << self.u
        p_file = File(self.folder_path+'/p/p.pvd')
        p_file << self.p

    def f(self, X, p_X, c_HbX, c_HbY):
        '''Generation rate as defined in the model.
        
        X: gas species. (string)
        p_X: partial pressure of X.
        c_HbX: concentration of Hb(X).
        c_HbY: concentration of Hb(Y)
        '''
        c_t = self.params['c_t']
        if X == 'O2':
            beta_O2 = self.params['beta_O2']
            k_O2 = self.params['k_O2']
            k_prime_O2 = self.params['k_prime_O2']
            first_term = k_prime_O2*(c_t - c_HbX - c_HbY)*p_X
            alpha = 1
            second_term = -k_O2/beta_O2*c_HbX
            beta = 1
            return alpha*first_term + beta*second_term
        elif X == 'CO2':
            beta_CO2 = self.params['beta_CO2']
            k_CO2 = self.params['k_CO2']
            k_prime_CO2 = self.params['k_prime_CO2']
            first_term = k_prime_CO2*(c_t - c_HbX - c_HbY)*p_X
            alpha = 1
            second_term = -k_CO2/beta_CO2*c_HbX
            beta = 1
            return alpha*first_term + beta*second_term
        else:
            raise ValueError('Gas species in f must be O2 or CO2.')

    def g(self, X, p_X, c_HbX, c_HbY):
        '''Scaled generation rate as defined in the model.
        
        X: gas species. (string)
        p_X: partial pressure of X.
        c_HbX: concentration of Hb(X).
        c_HbY: concentration of Hb(Y)
        '''
        if X == 'O2':
            beta_O2 = self.params['beta_O2']
            return beta_O2*self.f(X, p_X, c_HbX, c_HbY)
        elif X == 'CO2':
            beta_CO2 = self.params['beta_CO2']
            return beta_CO2*self.f(X, p_X, c_HbX, c_HbY)
        else:
            raise ValueError('Gas species in f must be O2 or CO2.')

    def sim_bst(self, final_time, num_steps, hb=True):
        '''Solves the blood-side transport (BST) problem of the model.
        
        final_time: time to simulate in seconds. (float)
        num_steps: number of time steps. (int)
        hb: toggle effects of hemoglobin. (bool)
        '''
        # Instance parameters

        self.T = final_time
        self.N = num_steps
        dt = self.T/self.N

        beta_O2 = self.params['beta_O2']
        k_O2 = self.params['k_O2']
        k_prime_O2 = self.params['k_prime_O2']
        p_air_O2 = self.params['p_air_O2']
        d_ba_O2 = self.params['d_ba_O2']
        d_pla_O2 = self.params['d_pla_O2']

        beta_CO2 = self.params['beta_CO2']
        p_air_CO2 = self.params['p_air_CO2']
        d_ba_CO2 = self.params['d_ba_CO2']
        d_pla_CO2 = self.params['d_pla_CO2']

        h_ba = self.params['h_ba']

        # Instance function space for the multi-field problem

        element = VectorElement('P', tetrahedron, 1, dim=4)
        self.M_h = FunctionSpace(self.mesh, element)

        # Declare functions and test functions

        x = Function(self.M_h)
        x_nm1 = Function(self.M_h)
        self.p_O2, self.p_CO2, self.c_HbO2, self.c_HbCO2 = split(x)

        v, w, eta, xi = TestFunctions(self.M_h)

        n = FacetNormal(self.mesh)

        # Set initial conditions

        x_nm1 = project(
            Expression(
                ('p_O2_0', 'p_CO2_0', 'c_HbO2_0', 'c_HbCO2_0'),
                degree=1,
                p_O2_0=self.params['p_O2_in'],
                p_CO2_0=self.params['p_CO2_in'],
                c_HbO2_0=0,
                c_HbCO2_0=0
            ), self.M_h
        )
        p_O2_nm1, p_CO2_nm1, c_HbO2_nm1, c_HbCO2_nm1 = split(x_nm1)

        # Define residuals

        G_p_O2 = self.p_O2*v*dx
        G_p_O2 += dt*d_pla_O2*inner(grad(self.p_O2), grad(v))*dx
        G_p_O2 += dt*d_ba_O2/h_ba*self.p_O2*v*ds(3)
        G_p_O2 += dt*inner(self.u*self.p_O2, grad(v))*dx
        G_p_O2 += -dt*inner(self.u*self.p_O2, n)*v*ds(2)
        if hb:
            G_p_O2 += dt*self.f('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)*v*dx
        G_p_O2 += -dt*d_ba_O2/h_ba*p_air_O2*v*ds(3)
        G_p_O2 += -p_O2_nm1*v*dx

        G_p_CO2 = self.p_CO2*w*dx
        G_p_CO2 += dt*d_pla_CO2*inner(grad(self.p_CO2), grad(w))*dx
        G_p_CO2 += dt*d_ba_CO2/h_ba*self.p_CO2*w*ds(3)
        G_p_CO2 += dt*inner(self.u*self.p_CO2, grad(w))*dx
        G_p_CO2 += -dt*inner(self.u*self.p_CO2, n)*w*ds(2)
        if hb:
            G_p_CO2 += dt*self.f(
                'CO2', self.p_CO2, self.c_HbCO2, self.c_HbO2
            )*w*dx
        G_p_CO2 += -dt*d_ba_CO2/h_ba*p_air_CO2*w*ds(3)
        G_p_CO2 += -p_CO2_nm1*w*dx

        G_c_O2 = self.c_HbO2*eta*dx
        G_c_O2 += -c_HbO2_nm1*eta*dx
        if hb:
            G_c_O2 += dt*inner(self.u*self.c_HbO2, grad(eta))*dx
            G_c_O2 += -dt*inner(self.u*self.c_HbO2, n)*eta*ds(2)
            G_c_O2 += -dt*self.g(
                'O2', self.p_O2, self.c_HbO2, self.c_HbCO2
            )*eta*dx

        G_c_CO2 = self.c_HbCO2*xi*dx
        G_c_CO2 += -c_HbCO2_nm1*xi*dx
        if hb:
            G_c_CO2 += dt*inner(self.u*self.c_HbCO2, grad(xi))*dx
            G_c_CO2 += -dt*inner(self.u*self.c_HbCO2, n)*xi*ds(2)
            G_c_CO2 += -dt*self.g('CO2', self.p_CO2, self.c_HbCO2, self.c_HbO2)*xi*dx

        G = G_p_O2 + G_p_CO2 + G_c_O2 + G_c_CO2

        #  Write useful comments for the simulation info file

        with open(self.folder_path+'/info.txt', 'w') as file:
            file.write('Full implicit')  # This should have a better interface

        # Create files for output

        p_O2_file = File(self.folder_path+'/bst/pO2.pvd')
        p_CO2_file = File(self.folder_path+'/bst/pCO2.pvd')
        c_HbO2_file = File(self.folder_path+'/bst/cHbO2.pvd')
        c_HbCO2_file = File(self.folder_path+'/bst/cHbCO2.pvd')

        # Time-stepping

        t = 0
        for n in range(self.N):
            if n == 0:
                # Save initial conditions
                
                _p_O2, _p_CO2, _c_HbO2, _c_HbCO2 = x.split()
                _p_O2.assign(project(p_O2_nm1, self.W_h))
                _p_CO2.assign(project(p_CO2_nm1, self.W_h))
                _c_HbO2.assign(project(c_HbO2_nm1, self.W_h))
                _c_HbCO2.assign(project(c_HbCO2_nm1, self.W_h))
                p_O2_file << (_p_O2, t)
                p_CO2_file << (_p_CO2, t)
                c_HbO2_file << (_c_HbO2, t)
                c_HbCO2_file << (_c_HbCO2, t)

            # Update current time

            t += dt

            # Declare Dirichlet boundary conditions for (BST)

            self.bst_dbc = [
                DirichletBC(
                    self.M_h.sub(0), self.params['p_O2_in'], self.gamma_in
                ),
                DirichletBC(
                    self.M_h.sub(1), self.params['p_CO2_in'], self.gamma_in
                ),
                DirichletBC(
                    self.M_h.sub(2), Constant(0), self.gamma_in
                ),
                DirichletBC(
                    self.M_h.sub(3), Constant(0), self.gamma_in
                )
            ]

            # Solve variational problem for time step

            solve(G == 0, x, self.bst_dbc)

            # Save solution to files

            _p_O2, _p_CO2, _c_HbO2, _c_HbCO2 = x.split()
            p_O2_file << (_p_O2, t)
            p_CO2_file << (_p_CO2, t)
            c_HbO2_file << (_c_HbO2, t)
            c_HbCO2_file << (_c_HbCO2, t)

            # Update previous solution

            x_nm1.assign(x)

            # Update progress

            print(
                f'Finished time step {n+1}/{self.N} ({round(t/self.T*100)}%)\n'
            )

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
            near(x[0], self.dir_min, self.tol) \
                or near(x[0], self.dir_max, self.tol)
        )

class GammaPi(SubDomain):
    '''Subdomain class for periodic boundary conditions.'''
    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of periodic direction. (float)
        dir_max: maximum value of periodic direction. (float)
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
        return on_boundary and near(x[1], self.dir_max, self.tol)

    def map(self, x, y):
        '''Maps opposite faces for periodic boundary conditions.
        
        x: position in subdomain.
        y: position in opposite subdomain.
        '''
        y[0] = x[0]
        y[1] = x[1] - (self.dir_max - self.dir_min)
        y[2] = x[2]