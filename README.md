# Computational modeling of capillary perfusion and gas exchange in alveolar tissue

Pablo Zurita Soler (_[@pzuritas](https://github.com/pzuritas)_), Daniel E. Hurtado (_[@dehurtado](https://github.com/dehurtado)_)

## Abstract

Gas exchange is an essential function of the respiratory system that couples fundamentally with perfusion in respiratory alveoli. Current mathematical formulations and computational models of these two phenomena rely on one-dimensional approximations that neglect the intricate volumetric microstructure of alveolar structures. In this work, we introduce a coupled three-dimensional computational model of pulmonary capillary perfusion and gas exchange that conforms to alveolar morphology. To this end, we derive non-linear partial differential equations and boundary conditions from physical principles that govern the behavior of blood and gases in arbitrary alveolar domains. We numerically solve the resulting formulation by proposing and implementing a non-linear finite-element scheme. Further, we carry out several numerical experiments to validate our model against one-dimensional simulations and demonstrate its applicability to morphologically-inspired geometries. Numerical simulations show that our model predicts blood pressure drops and blood velocities expected in the pulmonary capillaries. Moreover, we replicate partial pressure dynamics of oxygen and carbon dioxide reported in previous studies. This overall behavior is also observed in three-dimensional alveolar geometries, providing more detail associated with the spatial distribution of fields of interest and the influence of the shape of the domain. We envision that this model opens the door for enhanced _in silico_ studies of gas exchange and perfusion on realistic geometries, coupled models of respiratory mechanics and gas exchange, and multi-scale analysis of lung function; furthering our understanding of lung physiology and pathology.

## Directories

- `raw-data`: Data generated from direct simulations.
- `src`: Source files.
- `tests`: Use examples and tests of the model.

## Dependencies

Coded in Python 3.8.2 64-bit.

- [`FEniCS`](https://fenicsproject.org/) 2019.1.0
- `numpy`
= `pyvista` & `pythreejs` for Jupyter Notebook VTK visualization
- `datetime`
- `os`
- `sys`