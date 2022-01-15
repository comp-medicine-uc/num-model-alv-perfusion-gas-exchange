'''Model file for alveolar perfusion and gas exchange simulations.
'''

__author__ = 'pzuritas'
__email__ = 'pzurita@uc.cl'

import os
import time as tm
from datetime import date
import numpy as np
from fenics import *
from dolfin import *
#from mshr import Cylinder, generate_mesh
from src.boundaries import *


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

    def import_mesh(self, mesh_path, meshtype, type="h5", periodic=False):
        '''Imports mesh from .h5 file for use in simulations.

        mesh_path: path to file. (string)
        periodic: use periodic boundary conditions. (bool)
        '''
        if type == "h5":
            hdf5 = HDF5File(MPI.comm_world, mesh_path, 'r')
            self.mesh = Mesh()
            hdf5.read(self.mesh, 'mesh', False)
        elif type == "xml":
            self.mesh = Mesh(mesh_path)
        elif type == "xdmf":
            self.mesh = Mesh()
            with XDMFFile(mesh_path) as infile:
                infile.read(self.mesh)
        else:
            raise ValueError("type of mesh should be h5 or xml")

        if meshtype == 'sphere':
            marker = MeshFunction('bool', self.mesh, dim=3)
            marker.set_all(False)
            for cell in cells(self.mesh):
                coords = cell.get_vertex_coordinates()
                for i, coord in enumerate(coords):
                    if i % 3 == 0 and (coord < -40 or coord > 45):
                        marker[cell] = True
            self.mesh = refine(self.mesh, marker)

        if meshtype == 'tkd':
            marker = MeshFunction('bool', self.mesh, dim=3)
            marker.set_all(False)
            for cell in cells(self.mesh):
                coords = cell.get_vertex_coordinates()
                for i, coord in enumerate(coords):
                    if i % 3 == 0 and (coord < -43 or coord > 45):
                        marker[cell] = True
            self.mesh = refine(self.mesh, marker)

        dir_arr = np.array(
                [coords[0] for coords in self.mesh.coordinates()]
            )  # Node coordinates for principal (x) direction
        self.dir_max = np.max(dir_arr)  # Maximum principal direction coordinate
        self.dir_min = np.min(dir_arr)  # Minimum principal direction coordinate

        # Save mesh dims for other uses

        self.dims = (self.dir_max, self.dir_max, self.dir_max)

        # Flag for periodicity

        self.periodic = periodic

    def generate_slab_mesh(
        self, dims, elems, save=True, periodic=False, refined=True
    ):
        '''Generates a rectangular prism mesh for simulations on a slab.

        dims: dimensions of the mesh. (tuple)
        elems: number of elements in each direction. (tuple)
        save: save the generated mesh on a .pvd file. (bool)
        periodic: use periodic boundary conditions on the z direction. (bool)
        '''
        geometry = BoxMesh(
            Point(0, 0, 0), Point(*dims), *elems
        )
        self.mesh = geometry

        if refined:
            marker = MeshFunction('bool', self.mesh, dim=3)
            marker.set_all(False)
            for cell in cells(self.mesh):
                coords = cell.get_vertex_coordinates()
                for i, coord in enumerate(coords):
                    if i % 3 == 0 and coord < -50:
                        marker[cell] = True
            self.mesh = refine(self.mesh, marker)

        mesh_file = File(self.folder_path+'/mesh.pvd')
        mesh_file << self.mesh
        
        # Assign values to flow direction range

        self.dir_min = 0
        self.dir_max = dims[0]
        self.dims = dims

        self.periodic = periodic

    def instance_boundaries(self, mesh=None):
        '''Instances the relevant boundaries for boundary conditions.
        
        mesh: type of mesh created. None, "slab" or "tkd". (None or str)
        '''

        self.meshtype = mesh

        # Instance the relevant boundaries

        if not self.periodic:
            if not mesh == "sphere":
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
                self.boundaries.set_all(0)
                self.gamma_in.mark(self.boundaries, 1)
                self.gamma_out.mark(self.boundaries, 2)
                self.gamma_air.mark(self.boundaries, 3)

            elif mesh == "sphere":

                ## large sphere
                #self.gamma_in = GammaInSphereV2(220)
                #self.gamma_out = GammaOutSphereV2(220)
                #self.gamma_air = GammaAirSphereV2(211.99)

                ## small sphere
                self.gamma_in = GammaInSphereV2(60)
                self.gamma_out = GammaOutSphereV2(60)
                self.gamma_air = GammaAirSphereV2(58)

                # Declare the boundaries in the mesh and tag them

                self.boundaries = MeshFunction('size_t', self.mesh, dim=2)
                self.boundaries.set_all(0)
                self.gamma_in.mark(self.boundaries, 1)
                self.gamma_out.mark(self.boundaries, 2)
                self.gamma_air.mark(self.boundaries, 3)

        else:
            if mesh == "slab":
                self.gamma_in = GammaIn(
                    self.dir_min, self.dir_max, DOLFIN_EPS
                )
                self.gamma_out = GammaOut(
                    self.dir_min, self.dir_max, DOLFIN_EPS
                )
                self.gamma_pi = GammaSlabPi(0, self.dims[2], DOLFIN_EPS)
                self.gamma_air = GammaAirSlabPi(0, self.dims[1], DOLFIN_EPS)
            elif mesh == "tkd":
                #self.gamma_in = GammaIn(
                #    self.dir_min, self.dir_max, 1
                #)
                #self.gamma_out = GammaOut(
                #    self.dir_min, self.dir_max, 1
                #)
                self.gamma_in = GammaTKDIn(
                    self.dir_min, self.dir_max, 1E0
                )
                self.gamma_out = GammaTKDOut(
                    self.dir_min, self.dir_max, 1E0
                )
                self.gamma_pi = GammaTKDPi(
                    self.dir_min, self.dir_max, 1E0
                )
                self.gamma_air = GammaAirTKD(
                    self.dir_min, self.dir_max, 3.9599E-8
                )
                gamma_outish = GammaTKDOut(
                    self.dir_min, self.dir_max, 1E1
                )
                refine_out = MeshFunction('bool', self.mesh, dim=3)
                refine_out.set_all(False)
                gamma_outish.mark(refine_out, True)
                self.mesh = refine(self.mesh, refine_out)
            else:
                raise ValueError(
                    "Mesh type must be slab or tkd for periodicity."
                )
            
            # Declare the boundaries in the mesh and tag them

            self.boundaries = MeshFunction('size_t', self.mesh, dim=2)
            self.boundaries.set_all(0)
            self.gamma_in.mark(self.boundaries, 1)
            self.gamma_out.mark(self.boundaries, 2)
            self.gamma_air.mark(self.boundaries, 3)
            self.gamma_pi.mark(self.boundaries, 4)

            boundaries = File(self.folder_path+'/bnd/bnd.pvd')
            boundaries << self.boundaries

    def instance_function_spaces(self):
        '''Instances the relevant function spaces.'''

        if not self.periodic:
            self.W_h = FunctionSpace(self.mesh, 'Lagrange', 2)
            self.V_h = VectorFunctionSpace(self.mesh, 'Lagrange', 1)
        else:
            self.W_h = FunctionSpace(
                self.mesh, 'Lagrange', 2,
                constrained_domain=self.gamma_pi
            )
            self.V_h = VectorFunctionSpace(
                self.mesh, 'Lagrange', 1,
                constrained_domain=self.gamma_pi
            )

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

        self.periodic = False

    def sim_p(self, save=True, meshtype=None):
        '''Solves the perfusion (P) problem of the model.
        
        save: saves to vtk. (bool)
        meshtype: type of mesh. None, "slab" or "tkd". (None or str)
        '''
        self.instance_boundaries(mesh=meshtype)
        self.instance_function_spaces()

        # Declare Dirichlet boundary conditions for (P)

        self.p_dbc = [
            #DirichletBC(self.W_h, self.params['p_max'], self.gamma_in)#,
            DirichletBC(self.W_h, self.params['p_min'], self.gamma_out)
        ]

        # Assemble problem

        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)

        p = TrialFunction(self.W_h)
        v = TestFunction(self.W_h)
        a = inner(grad(p), grad(v))*dx
        F = Constant(0)*v*dx
        #if self.meshtype == 'sphere':
        F += 200*self.params["mu"]/self.params["kappa"]*v*ds(1)
        '''
        elif self.meshtype == 'tkd':
            c = 100/3
            x_hat_1 = project(Expression(
                'x[1]', degree=2
            ), self.W_h)
            file = File(self.folder_path+'/test-bc/x_hat_1.pvd')
            file << x_hat_1
            x_hat_2 = project(Expression(
                'x[2]', degree=2
            ), self.W_h)
            file = File(self.folder_path+'/test-bc/x_hat_2.pvd')
            file << x_hat_2
            x_hat_1_norm = project(abs(x_hat_1) + abs(x_hat_2), self.W_h)
            file = File(self.folder_path+'/test-bc/x_hat_1_norm.pvd')
            file << x_hat_1_norm
            n_expr = project(
                Expression(
                    ('1', '0', '0'), degree=2
                ), self.V_h
            )
            x_hat = project(Expression(
                ('0', 'x[1]', 'x[2]'), degree=2
            ), self.V_h)
            x_hat_normed = project(x_hat/inner(x_hat, x_hat)**0.5, self.V_h)
            file = File(self.folder_path+'/test-bc/x_hat_normed.pvd')
            file << x_hat_normed
            n_min_x_hat = project(n_expr - x_hat_normed, self.V_h)
            file = File(self.folder_path+'/test-bc/x_hat.pvd')
            file << x_hat
            file = File(self.folder_path+'/test-bc/n_min_x_hat.pvd')
            file << n_min_x_hat
            norm_n_x_hat = project(
                inner(n_min_x_hat, n_min_x_hat)**(1/2),
                self.W_h
            )
            weight = project(
                1/(
                    x_hat_1_norm/c/norm_n_x_hat + 1 - x_hat_1_norm/c
                )*(
                    x_hat_1_norm/c/norm_n_x_hat
                ), self.W_h
            )
            F += weight*200*self.params["mu"]/self.params["kappa"]*v*ds(1)
        '''

        # Solve problem

        self.p = Function(self.W_h)
        solve(
            a == F, self.p, self.p_dbc#,
            #solver_parameters={
            #    'linear_solver': 'gmres',
            #    'preconditioner': 'ilu'
            #}
        )

        self.u = project(
            -1/self.params['mu']*self.params['kappa']*grad(self.p),
            self.V_h
        )

        self.u.rename("u", "blood velocity [um/s]")
        self.p.rename("p", "blood pressure [mmHg]")

        if save:

            # Save solutions

            u_file = File(self.folder_path+'/p/u.pvd')
            u_file << self.u
            p_file = File(self.folder_path+'/p/p.pvd')
            p_file << self.p

    def set_u(self, value=(0, 0, 0), mesh='slab', save=True):
        '''Prescribes a velocity field u to the mesh instead of solving (P).

        value: uniform velocity field. (tuple)
        save: saves to vtk. (bool)
        '''
        self.instance_boundaries(mesh=mesh)
        self.instance_function_spaces()

        self.u = project(Expression(tuple(map(str, value)), degree=1), self.V_h)
        self.u.rename("u", "Blood velocity")

        if save:

            u_file = File(self.folder_path+'/p/u.pvd')
            u_file << self.u

    def f(self, X, p_X, c_HbX, c_HbY):
        '''Generation rate as defined in the model.
        
        X: gas species. (string)
        p_X: partial pressure of X. (FEniCS Function)
        c_HbX: concentration of Hb(X). (FEniCS Function)
        c_HbY: concentration of Hb(Y). (FEniCS Function)
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
        p_X: partial pressure of X. (FEniCS Function)
        c_HbX: concentration of Hb(X). (FEniCS Function)
        c_HbY: concentration of Hb(Y). (FEniCS Function)
        '''
        if X == 'O2':
            beta_O2 = self.params['beta_O2']
            return beta_O2*self.f(X, p_X, c_HbX, c_HbY)
        elif X == 'CO2':
            beta_CO2 = self.params['beta_CO2']
            return beta_CO2*self.f(X, p_X, c_HbX, c_HbY)
        else:
            raise ValueError('Gas species in f must be O2 or CO2.')

    def sim_bst(self, final_time, num_steps, hb=True, save=True):
        '''Solves the blood-side transport (BST) problem of the model.
        
        final_time: time to simulate in seconds. (float)
        num_steps: number of time steps. (int)
        hb: toggle effects of hemoglobin. (bool)
        save: saves to vtk. (bool)
        '''
        # Instance parameters

        self.T = final_time
        self.N = num_steps
        dt = self.T/self.N

        p_air_O2 = self.params['p_air_O2']
        d_ba_O2 = self.params['d_ba_O2']
        d_pla_O2 = self.params['d_pla_O2']

        p_air_CO2 = self.params['p_air_CO2']
        d_ba_CO2 = self.params['d_ba_CO2']
        d_pla_CO2 = self.params['d_pla_CO2']

        h_ba = self.params['h_ba']

        # Instance function space for the multi-field problem

        element = VectorElement('P', tetrahedron, 1, dim=4)
        self.M_h = FunctionSpace(self.mesh, element)
        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)

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

        if save:

            # Create files for output

            p_O2_file = File(self.folder_path+'/bst/pO2.pvd')
            p_CO2_file = File(self.folder_path+'/bst/pCO2.pvd')
            c_HbO2_file = File(self.folder_path+'/bst/cHbO2.pvd')
            c_HbCO2_file = File(self.folder_path+'/bst/cHbCO2.pvd')

        # Time-stepping

        t = 0
        for n in range(self.N):
            if n == 0:
                if save:

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

            if save:

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

    def sim_sbst(self, hb=True, save=True, guess=None):
        '''Solves the steady state blood-side transport (SBST) problem of the
        model.
        
        hb: toggle effects of hemoglobin. (bool)
        save: saves to vtk. (bool)
        '''
        # Instance parameters

        p_air_O2 = self.params['p_air_O2']
        d_ba_O2 = self.params['d_ba_O2']
        d_pla_O2 = self.params['d_pla_O2']

        p_air_CO2 = self.params['p_air_CO2']
        d_ba_CO2 = self.params['d_ba_CO2']
        d_pla_CO2 = self.params['d_pla_CO2']

        h_ba = self.params['h_ba']

        # Instance function space for the multi-field problem

        element = VectorElement('P', tetrahedron, 1, dim=4)
        self.M_h = FunctionSpace(self.mesh, element)
        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)

        # Declare functions and test functions

        if not guess:
            x = Function(self.M_h)
        else:
            x = project(guess, self.M_h)
        self.p_O2, self.p_CO2, self.c_HbO2, self.c_HbCO2 = split(x)

        v, w, eta, xi = TestFunctions(self.M_h)

        n = FacetNormal(self.mesh)

        # Define residuals

        G_p_O2 = d_pla_O2*inner(grad(self.p_O2), grad(v))*dx
        G_p_O2 += -d_ba_O2/h_ba*Constant(p_air_O2)*v*ds(3)
        G_p_O2 += d_ba_O2/h_ba*self.p_O2*v*ds(3)
        G_p_O2 += inner(self.p_O2*self.u, n)*v*ds(2)
        G_p_O2 += -inner(self.p_O2*self.u, grad(v))*dx
        if hb:
            G_p_O2 += self.f('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)*v*dx

        G_p_CO2 = d_pla_CO2*inner(grad(self.p_CO2), grad(w))*dx
        G_p_CO2 += -d_ba_CO2/h_ba*Constant(p_air_CO2)*w*ds(3)
        G_p_CO2 += d_ba_CO2/h_ba*self.p_CO2*w*ds(3)
        G_p_CO2 += inner(self.p_CO2*self.u, n)*w*ds(2)
        G_p_CO2 += -inner(self.p_CO2*self.u, grad(w))*dx
        if hb:
            G_p_CO2 += self.f('CO2', self.p_CO2, self.c_HbCO2, self.c_HbO2)*w*dx

        if hb:
            G_c_O2 = self.g('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)*eta*dx
            G_c_O2 += inner(self.c_HbO2*self.u, grad(eta))*dx
            G_c_O2 += -inner(self.c_HbO2*self.u, n)*eta*ds(2)

            G_c_CO2 = self.g(
                'CO2', self.p_CO2, self.c_HbCO2, self.c_HbO2
            )*xi*dx
            G_c_CO2 += inner(self.c_HbCO2*self.u, grad(xi))*dx
            G_c_CO2 += -inner(self.c_HbCO2*self.u, n)*xi*ds(2)
        else:
            G_c_O2 = self.c_HbO2*eta*dx
            G_c_CO2 = self.c_HbCO2*xi*dx

        '''
        if self.meshtype == 'tkd':
            # Define outlet boundary condition for gases
            c = 100/3
            x_hat_1 = project(Expression(
                'x[1]', degree=2
            ), self.W_h)
            file = File(self.folder_path+'/test-bc/x_hat_1.pvd')
            file << x_hat_1
            x_hat_2 = project(Expression(
                'x[2]', degree=2
            ), self.W_h)
            file = File(self.folder_path+'/test-bc/x_hat_2.pvd')
            file << x_hat_2
            x_hat_1_norm = project(abs(x_hat_1) + abs(x_hat_2), self.W_h)
            file = File(self.folder_path+'/test-bc/x_hat_1_norm.pvd')
            file << x_hat_1_norm
            n_expr = project(
                Expression(
                    ('1', '0', '0'), degree=2
                ), self.V_h
            )
            x_hat = project(Expression(
                ('0', 'x[1]', 'x[2]'), degree=2
            ), self.V_h)
            x_hat_normed = project(x_hat/inner(x_hat, x_hat)**0.5, self.V_h)
            file = File(self.folder_path+'/test-bc/x_hat_normed.pvd')
            file << x_hat_normed
            n_min_x_hat = project(n_expr - x_hat_normed, self.V_h)
            file = File(self.folder_path+'/test-bc/x_hat.pvd')
            file << x_hat
            file = File(self.folder_path+'/test-bc/n_min_x_hat.pvd')
            file << n_min_x_hat
            norm_n_x_hat = project(
                inner(n_min_x_hat, n_min_x_hat)**(1/2),
                self.W_h
            )
            weight = project(
                1/(
                    x_hat_1_norm/c/norm_n_x_hat + 1 - x_hat_1_norm/c
                )*(
                    x_hat_1_norm/c/norm_n_x_hat
                ), self.W_h
            )
            G_p_O2 += weight*d_pla_O2*inner(grad(self.p_O2), x_hat)*v*ds(2)
            G_p_CO2 += weight*d_pla_CO2*inner(grad(self.p_CO2), x_hat)*w*ds(2)
        '''

        G = G_p_O2 + G_p_CO2 + G_c_O2 + G_c_CO2

        if save:

            # Create files for output

            p_O2_file = File(self.folder_path+'/sbst/pO2.pvd')
            p_CO2_file = File(self.folder_path+'/sbst/pCO2.pvd')
            c_HbO2_file = File(self.folder_path+'/sbst/cHbO2.pvd')
            c_HbCO2_file = File(self.folder_path+'/sbst/cHbCO2.pvd')

        # Declare Dirichlet boundary conditions for (SBST)

        self.sbst_dbc = [
            DirichletBC(
                self.M_h.sub(0), Constant(self.params['p_O2_in']), self.gamma_in
            ),
            DirichletBC(
                self.M_h.sub(1), Constant(self.params['p_CO2_in']),
                self.gamma_in
            ),
            DirichletBC(
                self.M_h.sub(2), Constant(self.params['c_HbO2_in']),
                self.gamma_in
            ),
            DirichletBC(
                self.M_h.sub(3), Constant(self.params['c_HbCO2_in']),
                self.gamma_in
            )
        ]

        # Solve variational problem

        solve(
            G == 0, x, self.sbst_dbc,
            solver_parameters={"newton_solver": {
                "relative_tolerance": 1E-8,
                "absolute_tolerance": 1E-8#,
                #"linear_solver": "gmres",
                #"preconditioner": "ilu"
            }}
        )

        if save:

            # Save solution to files

            _p_O2, _p_CO2, _c_HbO2, _c_HbCO2 = x.split()
            _p_O2.rename("p_O2", "oxygen partial pressure [mmHg]")
            _p_CO2.rename("p_CO2", "carbon dioxide partial pressure [mmHg]")
            _c_HbO2.rename("c_HbO2", "oxyhemoglobin concentration [mole/um^3]")
            _c_HbCO2.rename(
                "c_HbCO2", "carboaminohemoglobin concentration [mole/um^3]"
            )
            p_O2_file << _p_O2
            p_CO2_file << _p_CO2
            c_HbO2_file << _c_HbO2
            c_HbCO2_file << _c_HbCO2

        return x

    def avg_field_along_dir(self, u_0, n):
        '''Calculates an average field along a certain direction for a slab
        mesh. Specifically, defines the function v(x) = int_{S(x)} u_0 ds where
        S(x) = {(x1, y2, y3) | y2, y3 in [0, 6]}.

        u_0: field. (FEniCS Function)
        n: number of slices (more than 1) (discretization of x). (int)
        '''
        # https://mscroggs.github.io/qa/3257/
        # building-2d-function-g-x-y-from-a-3d-function-f-x-y-z-constant/

        x_arr = np.linspace(0.0, self.dir_max - self.dir_min, n)
        v = np.zeros((n,))

        # Create boundary mesh
        bmesh = BoundaryMesh(self.mesh, "exterior")

        # Create SubMesh for side at z=0
        # This will be a UnitSquareMesh with topology dimension 2 in 3 space dimensions
        cc = MeshFunction('size_t', bmesh, dim=2)
        xyplane = AutoSubDomain(lambda x: x[0] < 1e-8)
        xyplane.mark(cc, 6)
        submesh = SubMesh(bmesh, cc, 6)

        for i in range(n):

            # Move slice/submesh to z=0.5
            x = submesh.coordinates()
            x[:, 2] += x_arr[i]

            # Create a FunctionSpace on the submesh
            Vs = FunctionSpace(submesh, "Lagrange", 1)

            # interpolate_nonmatching_mesh required in parallel, 
            # interpolate works in series
            #us = interpolate_nonmatching_mesh(u, Vs)
            us = interpolate(u_0, Vs)
            v[i] = assemble(us*dx)

        return [v, x]