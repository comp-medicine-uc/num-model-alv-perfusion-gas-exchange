'''Main file for alveolar perfusion and gas exchange simulations in a sheet
mesh.
'''

__author__ = 'pzuritas'
__email__ = 'pzurita@uc.cl'

import sys
import os
sys.path.append(os.getcwd()[:-6])
import dolfin
from src.model import PerfusionGasExchangeModel
from src.params import params


print("Starting...")
folder = "slab_job"
path = os.path.join("../raw-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params)
print("Model initialised")
model.generate_slab_mesh(
    dims=(200, 6, 6), elems=(200, 6, 6), save=True, periodic=True, refined=True
)
print("Mesh generated")
model.mesh = dolfin.refine(model.mesh)
model.mesh = dolfin.refine(model.mesh)
print("Mesh refined")
print("Setting u")
model.set_u(value=(800/3, 0, 0), save=True)
print("u set")
print("Starting transport simulation")
x = model.sim_sbst(hb=False, save=False)
solution = model.sim_sbst(hb=True, save=True, guess=x)
print("Airflow:", model.compute_airflow())
print("Conservation:", model.compute_blood_conservation())
print("Done")