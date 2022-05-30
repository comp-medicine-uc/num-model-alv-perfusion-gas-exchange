'''Main file for alveolar perfusion and gas exchange simulations in spherical
mesh.
'''

__author__ = 'pzuritas'
__email__ = 'pzurita@uc.cl'

import sys
import os
# The following line adds the directory to the path in order to cross-reference
# files in the repo
sys.path.append(os.getcwd()[:-14])
import dolfin
from src.model import PerfusionGasExchangeModel
from src.params import params


print("Starting...")
folder = "sphere-job"
path = os.path.join("../../raw-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params)
print("Model initialised")
model.import_mesh(
    os.path.join("../../raw-data/meshes", "sphere.xml"), meshtype='sphere',
    type="xml"
)
print("Mesh imported")
model.mesh = dolfin.refine(model.mesh)
print("Mesh refined")
print("Starting (P) simulation")
model.sim_p(save=True, meshtype="sphere")
print("(P) simulation done")
print("Starting (T)) simulation")
x = model.sim_sbst(hb=False, save=False)
solution = model.sim_sbst(hb=True, save=True, guess=x)
print("Done")