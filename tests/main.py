from src.model import PerfusionGasExchangeModel
from src.params import params
import os


folder = "test-folder"
path = os.path.join("raw-data", folder)

model = PerfusionGasExchangeModel(folder_path=path, params=params)
model.generate_slab_mesh(dims=(200, 6, 200), elems=(40, 3, 30), save=True)
#model.generate_cylinder_mesh(end=(100, 0, 0), r=3, save=True)
model.sim_p()
model.sim_bst(final_time=1, num_steps=30)