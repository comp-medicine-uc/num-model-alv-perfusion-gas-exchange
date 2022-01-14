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


folder = "sphere_job"
path = os.path.join("../raw-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params)
model.import_mesh(
    os.path.join("../raw-data", "sphere_small.xml"), meshtype='sphere',
    type="xml"
)
model.mesh = dolfin.refine(model.mesh)
model.mesh = dolfin.refine(model.mesh)
model.sim_p(save=True, meshtype="sphere")
x = model.sim_sbst(hb=False, save=False)
solution = model.sim_sbst(hb=True, save=True, guess=x)