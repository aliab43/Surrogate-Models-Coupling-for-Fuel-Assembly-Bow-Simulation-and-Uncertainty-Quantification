#!/bin/bash

ln -s /volatile/catC/aa270113/FinalSimulation/FSI_Coupling3/TGaussianProcess_Hydraulic
ln -s /volatile/catC/aa270113/FinalSimulation/FSI_Coupling3/TGaussianProcess_matern12_J7
ln -s /volatile/catC/aa270113/FinalSimulation/FSI_Coupling3/TGaussianProcess_matern12_J41
ln -s /volatile/catC/aa270113/FinalSimulation/FSI_Coupling3/TGaussianProcess_matern12_V


python3 MecaUQ.py

python3 HydrauUQ.py

python3 ~/FinalSimulation/FSI_Coupling3/SimulateOneCycleDef.py