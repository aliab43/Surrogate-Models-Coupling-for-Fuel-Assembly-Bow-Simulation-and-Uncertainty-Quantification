# FSI_Coupling3

`FSI_Coupling3/` contains the main coupled mechanical-hydraulic surrogate workflow for fuel assembly bow simulation over a reactor cycle.

## Main Script

- `SimulateOneCycleDef.py`: runs the three main phases of the coupled simulation:
  - before irradiation
  - irradiation
  - after irradiation

The driver combines:

- mechanical surrogate predictions for modal coefficients
- hydraulic surrogate predictions for force fields
- contact and penetration filtering
- grid-clamping evolution during irradiation

## Core Modules

- `ClassMechanicalSurrogates.py`: wrapper around the mechanical Gaussian-process surrogate models.
- `ClassHydraulicSurrogate.py`: loads the hydraulic surrogate models and reconstructs the force matrix.
- `PenetrationDetection.py`: computes filtered deformation fields and applies contact constraints.
- `GridClampingUpdate.py`: updates clamping forces during irradiation.
- `post_traitement.py`: post-processing helpers and plotting.
- `HydrauUQ.py`, `MecaUQ.py`, `uncertaintyProp.py`: uncertainty propagation utilities.

## Supporting Utilities

- `cycle.py`: generates the illustrative reactor-cycle displacement figure.
- `generate_hydrau_DOE.py`: builds hydraulic DOE samples in `.dat` and `.npy` formats.
- `generate_doe_meca.py`: builds the mechanical DOE used to train or refresh the surrogates.
- `plot_*` scripts: standalone visualization helpers for selected study outputs.

## Versioned Files in This Folder

This folder contains the source code plus a small set of lightweight reference files and figures, including:

- `DOE_for_mechanics/Doe_for_surrogate.txt`
- `Matrix_M_N.dat`
- `Matrix_M_N.root`
- selected `.png` figures used for illustration

These files document the workflow, but they are not enough to reproduce the full production run on their own.

## External Assets Required for a Full Run

The complete coupled simulation depends on large local assets that are not versioned in Git, such as:

- `Doe_for_mechanics.npy`
- `Doe_for_hydraulic.npy`
- kriging export directories such as `TGaussianProcess_*`

Without those files, this folder should be read as a documented research software snapshot rather than a fully runnable public package.

## Generated Outputs

The following files and folders are generated locally during runs and are intentionally ignored by Git:

- `nohup.out`
- `output.txt`
- `final.txt`
- `displacement_before_irradiation.txt`
- `displacement_irradiation.txt`
- `displacement_after_irradiation.txt`
- `result_of_simulations/`
- `iteration_deformation_init/`

The main driver now writes its default text outputs relative to `FSI_Coupling3/`, so running it from the project root does not pollute the repository root with generated files.
