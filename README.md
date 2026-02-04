# Surrogate Models Coupling for Fuel Assembly Bow Simulation and Uncertainty Quantification

This repository contains a coupled mechanical–hydraulic surrogate workflow to simulate fuel assembly bow evolution over a full cycle (before irradiation, during irradiation, after irradiation), plus tools for uncertainty quantification and post‑processing.

The core idea is to replace expensive high‑fidelity codes with Gaussian‑process (GP) surrogate models, then couple them in a fixed‑point loop with contact/penetration detection and grid‑clamping evolution. The workflow produces assembly bow fields and modal coefficients for each phase and supports Monte Carlo uncertainty propagation.

## Project Goals

- Build fast surrogate models for mechanical bow and hydraulic forces.
- Couple the two surrogates with a stable, convergent FSI‑like loop.
- Track bow evolution through three phases of a reactor cycle.
- Quantify uncertainty in outputs due to measurement and model variability.
- Provide reproducible post‑processing and validation plots.

## Main Entry Point

`FSI_Coupling3/SimulateOneCycleDef.py`

This script runs the full cycle:

1. **Pre‑irradiation** (fuel insertion and pump activation).
2. **Irradiation** (time stepping with grid clamping evolution).
3. **Post‑irradiation** (power‑off relaxation).

Outputs are modal coefficients and displacement text files for each phase.

## How the Coupled Simulation Works

The coupled loop follows a fixed‑point update:

1. **Mechanical state** starts from DOE data (`Doe_for_mechanics.npy`) with modal coefficients `C, S, W`.
2. **PenetrationDetection** computes a filtered deformation and enforces contact constraints with penalty forces.
3. **HydraulicSurrogateModel** predicts hydraulic forces from the current modal state.
4. Forces are projected onto the modal basis and fed back into the mechanical DOE.
5. A relaxation factor is applied to stabilize convergence.
6. Convergence is checked on the modal coefficients.

The irradiation phase also updates **grid clamping** to account for fluence‑dependent evolution.

## Key Classes and What They Do

`FSI_Coupling3/ClassMechanicalSurrogates.py`
- `MechanicalSurrogate`
- Wraps GP models for mechanical bow coefficients.
- Provides `callSurrogateModelsday7`, `callSurrogateModelCreep`, and `callSurrogateModelday41`.
- These correspond to different physical phases (early, creep/irradiation, late).

`FSI_Coupling3/ClassHydraulicSurrogate.py`
- `HydraulicSurrogateModel`
- Dynamically loads GP models for hydraulic forces at each grid/assembly.
- Produces a force matrix (grids × assemblies).

`FSI_Coupling3/PenetrationDetection.py`
- `PenetrationDetection`
- Builds assembly deformation fields, checks penetration/contact, and applies penalty forces.
- Updates mechanical DOE modal coefficients to enforce constraints.
- Can plot intermediate deformation (saved under `FSI_Coupling3/iteration_deformation_init/`).

`FSI_Coupling3/GridClampingUpdate.py`
- `GridClampingUpdater`
- Updates grid‑clamping force based on fluence evolution using interpolation.

## Data and DOE Inputs

The simulation expects DOE and surrogate data already in the repository:

- `FSI_Coupling3/DOE_for_mechanics/` contains:
  - `Doe_for_mechanics.npy`
  - `Doe_for_hydraulic.npy`
- GP model directories:
  - `FSI_Coupling3/TGaussianProcess_matern12_J7/`
  - `FSI_Coupling3/TGaussianProcess_matern12_V/`
  - `FSI_Coupling3/TGaussianProcess_matern12_J41/`
  - `FSI_Coupling3/TGaussianProcess_Hydraulic/`

If you regenerate DOE or surrogates, update these folders before running the cycle.

## Data Availability Note

The full fuel‑assembly case requires large covariance matrices and `.npy` datasets that are too heavy to store in this Git repository. If the cycle fails to run due to missing files, it is because those GP surrogate assets are not included here. This repository is intended to present the code, workflow, and development approach; the large surrogate data must be provided separately.

## Running the Cycle

From the project root:

```bash
python3 FSI_Coupling3/SimulateOneCycleDef.py
```

Main outputs:

- `FSI_Coupling3/displacement_before_irradiation.txt`
- `FSI_Coupling3/displacement_irradiation.txt`
- `FSI_Coupling3/displacement_after_irradiation.txt`

## Uncertainty Quantification

The repository includes scripts for uncertainty propagation and plotting:

- `FSI_Coupling3/uncertaintyProp.py`
- `FSI_Coupling3/PlotUncertaintySurrogate.py`
- `FSI_Coupling3/plot_initial_uncertainty.py`
- `FSI_Coupling3/post_traitement.py`

These scripts perform Monte Carlo sampling and generate plots for uncertainty bands and distributions.

## Benchmark Validation

A dedicated analytical benchmark validates the coupled GP uncertainty method:

- `FSI_Coupling3/Analetical_Benchmark_Coupled_GP_Validation.py`

See `README_BENCHMARK.md` for full details.

The analytical benchmark runs out‑of‑the‑box and does not require the large surrogate data files.

## Repository Structure

- `FSI_Coupling3/SimulateOneCycleDef.py`: full cycle driver.
- `FSI_Coupling3/PenetrationDetection.py`: contact detection and penalty update.
- `FSI_Coupling3/ClassMechanicalSurrogates.py`: mechanical GP surrogates.
- `FSI_Coupling3/ClassHydraulicSurrogate.py`: hydraulic GP surrogates.
- `FSI_Coupling3/GridClampingUpdate.py`: fluence‑based clamping evolution.
- `FSI_Coupling3/*UQ.py`: uncertainty workflow.
- `FSI_Coupling3/Analetical_Benchmark_Coupled_GP_Validation.py`: analytical validation.

## Requirements

Typical Python dependencies include:

- `numpy`
- `scipy`
- `matplotlib`
- `scikit-learn`

Create and activate your environment, then install as needed.

## Contact

For questions about the model assumptions, surrogate construction, or validation strategy, open an issue or contact the repository maintainer.

## Copyright

Copyright (c) 2026 Ali Abboud.  
Contact: ali.ib.abboud95@gmail.com, ali.abboud@polytechnique.edu
