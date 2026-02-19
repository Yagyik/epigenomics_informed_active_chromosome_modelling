# Chromosome Modelling Pipeline (Config-Driven)

## Setup and Environment

System dependencies used by simulation build/run paths:

- `mpicc` (OpenMPI)
- `gsl` (`-lgsl -lgslcblas`)
- `bc`
- C toolchain (`gcc/g++`, `make`)

Python dependencies are in `requirements.txt`.

Linux/macOS automated setup:

```bash
./scripts/setup_environment.sh
```

Python-only setup (skip system package install):

```bash
./scripts/setup_environment.sh --python-only
```

Venv-only manual setup:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Windows (PowerShell, pip/venv):

```powershell
.\scripts\setup_environment.ps1
```

Alternative virtualenv tool:

```bash
python3 -m pip install virtualenv
python3 -m virtualenv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Cross-platform notes:

- Best portability: use pip + virtualenv and run simulations where `mpicc` + `gsl` are available.
- On Windows, recommended simulation path is WSL2 (Ubuntu) plus `./scripts/setup_environment.sh`.
- Native Windows can use `scripts/setup_environment.ps1` for Python setup, but MPI/GSL builds are typically more reliable in WSL2 or Linux containers.

## To Test The Code

This test is a local non-Slurm test of simulations.

To test the code, run:

```bash
python3 run_test.py --config pipeline/pipeline_test_config.json
```

Short smoke test (single run, shorter runtime):

```bash
python3 run_test.py --config pipeline/pipeline_test_config.json --max-runs 1 --simlength 2000
```

Parameters for this simulation are in:

- `pipeline/params/dediff_simulation_test.json`

The stage wiring for the test run is in:

- `pipeline/pipeline_test_config.json` (`stage_params_files.dediff_simulation_test`)

Generated outputs are in newly created subfolders under:

- `dediffSimulation/Set_<set>_in<inSuffix>_<pref>_Tp<Tp>/<run_dir>/`

Examine:

- `dediffSimulation/Set_<set>_in<inSuffix>_<pref>_Tp<Tp>/<run_dir>/Conf-thermo.dat`
- Example: `dediffSimulation/Set_1_in3_ell_Tp0.0001/46ch_T0.0001_sga0.001_sgd0.001_epb8_epa6_epm4_lmn0.1_sfdf0.001_sgat2_as0/Conf-thermo.dat`

Use column 2 of `Conf-thermo.dat` to examine gross interchromosomal intermingling energy.

## Repository Description

This repository provides a staged, config-driven pipeline for:

1. `chrom_props` -> `ChromPartitionProps.dat`
2. `ipd_hic_basic` -> `IPD.dat`
3. `iad_rnaseq_hicmask` -> `IAD.dat`
4. `expected_mat` -> `ConstraintsVector.dat`
5. `init_positions_*` -> start/int/constraint files
6. `prep_surf_meshes*` -> prepped surface files
7. `dediff_simulation` / `anneal_simulation` -> Slurm submission script generation
8. `dediff_simulation_test` -> local non-Slurm dediff test run

Processed upstream outputs are centralized in:

- `processed_HiC_epigenomic_data/`

Canonical files:

- `ChromPartitionProps.dat`
- `IPD.dat`
- `IAD.dat`
- `ConstraintsVector.dat`

## Other Run Modes

Run test mode through the pipeline entrypoint:

```bash
python3 pipeline/run_pipeline.py --config pipeline/pipeline_test_config.json --mode test
```

Run full pipeline:

```bash
python3 pipeline/run_pipeline.py --config pipeline/pipeline_config.json
```

Dry-run:

```bash
python3 pipeline/run_pipeline.py --config pipeline/pipeline_config.json --dry-run
```

## Execution and Config Notes

- Stage orchestration is Python-first (`pipeline/run_pipeline.py`, stage wrappers, `run_test.py`).
- Simulation wrappers currently call shell launchers (`run_46ch_script_usePrepped.sh`, `run_46ch_script_usePrepped_reprog.sh`) per parameter combination.
- Compile toolchain is configured in stage params (`compiler`, set to `mpicc` for simulation params).

Config files:

- User-facing config: `pipeline/pipeline_config.json`
- Core workflow config: `pipeline/core_config.json`
- Stage parameter files: `pipeline/params/*.json`

`pipeline/pipeline_config.json` contains:

- `core_config`
- `paths`
- `files`
- `stage_params_files`
- `stages` (only `script`, `required_inputs`, `required_outputs`)

`pipeline/core_config.json` contains:

- `pipeline_modes`
- `naming`
- core chrom-props clustering settings
- IPD/IAD control settings

## Core Scripts

- `genChromProps/make_chrom_props.py`
- `genConstraintsInteractions/make_ipd.py`
- `genConstraintsInteractions/make_iad.py`
- `genInitPositions/vectorizeConstraints.py`
- `genInitPositions/run_genStart_ell.py`
- `genInitPositions/run_genStart_sph.py`
- `genSurfMeshes/run_prepSurf_ell.py`
- `genSurfMeshes/run_prepSurf_sph.py`
- `dediffSimulation/genSlurm_submission_script.py`
- `annealSimulation/genSlurm_submission_script.py`
- `run_test.py`
- `pipeline/run_pipeline.py`
