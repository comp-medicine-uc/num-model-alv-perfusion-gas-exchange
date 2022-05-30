'''Main file for alveolar perfusion and gas exchange simulations in a sheet
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
folder = "sheet-job"
path = os.path.join("../../raw-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params)
print("Model initialised")
model.generate_sheet_mesh(
    dims=(200, 6, 6), elems=(200, 6, 6), save=True, periodic=True, refined=True
)
print("Mesh generated")
model.mesh = dolfin.refine(model.mesh)
model.mesh = dolfin.refine(model.mesh)
print("Mesh refined")
print("Setting u")
model.set_u(value=(800/3, 0, 0), save=True)
print("u set")
print("Starting (T)) simulation")
x = model.sim_t(hb=False, save=False)
solution = model.sim_t(hb=True, save=True, guess=x)
print("Done")