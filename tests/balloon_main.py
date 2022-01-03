import sys
import os
import dolfin
sys.path.append(os.getcwd()[:-6])
from src.model import PerfusionGasExchangeModel
from src.params import params

mesh_path = "../raw-data/sphere_fine.xml"
folder = "sphere_fine_sims"
path = os.path.join("../raw-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params)
model.import_mesh(mesh_path, type="xml")
model.sim_p(save=True, meshtype="sphere")
boundaries = dolfin.File(model.folder_path+'/bnd/bnd.pvd')
boundaries << model.boundaries
x = model.sim_sbst(hb=False, save=False)
solution = model.sim_sbst(hb=True, save=True, guess=x)