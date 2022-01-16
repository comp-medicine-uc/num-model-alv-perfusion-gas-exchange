'''Main file for alveolar perfusion and gas exchange simulations in TKD
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
folder = "tkd_job"
path = os.path.join("../raw-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params)
print("Model initialised")
model.import_mesh(
    os.path.join("../raw-data", "TKD_new_smooth.xml"), meshtype='tkd', type="xml", 
    periodic=True
)
print("Mesh imported")
model.mesh = dolfin.refine(model.mesh)
model.mesh = dolfin.refine(model.mesh)
print("Mesh refined")
print("Starting (P) simulation")
model.sim_p(save=True, meshtype="tkd")
print("(P) simulation done")
print("Starting transport simulation")
x = model.sim_sbst(hb=False, save=False)
solution = model.sim_sbst(hb=True, save=True, guess=x)
print("Done")