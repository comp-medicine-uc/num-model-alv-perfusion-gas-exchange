'''Main file for alveolar perfusion and gas exchange simulations in spherical
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
folder = "sphere_job"
path = os.path.join("../raw-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params)
print("Model initialised")
model.import_mesh(
    os.path.join("../raw-data", "sphere_small.xml"), meshtype='sphere',
    type="xml"
)
print("Mesh imported")
model.mesh = dolfin.refine(model.mesh)
print("Mesh refined")
print("Starting (P) simulation")
model.sim_p(save=True, meshtype="sphere")
print("(P) simulation done")
print("Starting transport simulation")
x = model.sim_sbst(hb=False, save=False)
solution = model.sim_sbst(hb=True, save=True, guess=x)
print("Done")