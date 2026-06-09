# VASP C-O Distance Sampling and Free Energy Analysis

This repository contains a compact workflow for VASP structure sampling, batch optimization/molecular dynamics calculations, and post-processing of C-O distance-dependent force and free energy profiles.

The core project consists of six Python scripts. Starting from an initial `POSCAR`, the workflow generates perturbed structures, submits VASP optimization and MD jobs, collects `OUTCAR` files, extracts C-O distances/energies/force projections, and builds confidence-interval-filtered force and relative free energy curves.

The calculation workflow is designed for a VASP cluster environment using an LSF job scheduler. The analysis and plotting steps can also be run in a local Python environment after the VASP output files are available.

## Project Layout

```text
.
|-- 01_data_gen.py              # Generate randomly perturbed POSCAR files
|-- 02_data_gen.py              # Submit optimization/MD jobs and monitor them
|-- POSCAR                      # Initial structure
|-- inopt/                      # VASP input templates for optimization
|   |-- INCAR
|   |-- KPOINTS
|   |-- POTCAR
|   `-- vasp.script
|-- inmd/                       # VASP input templates for MD
|   |-- INCAR
|   |-- KPOINTS
|   |-- POTCAR
|   `-- vasp.script
|-- opt/                        # Optimization output directories
`-- std/
    |-- 03_data_gen.py          # Collect OUTCAR files from MD outputs
    `-- OUfolder/
        |-- 04_data_analysis.py # Extract C-O distance, energy, and force projection
        |-- 05_data_analysis.py # Scan CI windows and evaluate stability
        `-- 06_data_figures.py  # Generate CI-filtered force/free-energy figures
```

Directories such as `opt/`, `std/CONTCAR_*`, `std/OUfolder/OUTCAR_*`, and large VASP-generated files such as `WAVECAR`, `CHGCAR`, and `vasprun.xml` are calculation outputs rather than core source files.

The `example/` directory included in this repository provides a representative case for running and testing the workflow. It is intended as sample data to illustrate the expected directory structure, intermediate files, and analysis outputs.

`POTCAR` files are not uploaded to this repository. They are excluded because VASP pseudopotential files are subject to the VASP license; users should provide their own licensed `POTCAR` files when running new calculations.

## Workflow

### 1. Generate perturbed POSCAR files

`01_data_gen.py` reads the root-level `POSCAR`, applies random displacements to movable atoms, and uses a coordinate hash to avoid duplicate structures.

```bash
python 01_data_gen.py
```

Default output directory:

```text
new_POSCAR_folder/
```

By default, the script generates files from `POSCAR_61` to `POSCAR_100`. To change the generated range or displacement magnitude, edit:

- `target_start`
- `target_end`
- `random_move_vector(max_move=0.005)`
- `C_index`, `O1_index`, and `O2_index`

### 2. Run VASP optimization and MD jobs

`02_data_gen.py` reads a set of `POSCAR_*` files, submits optimization jobs, converts the optimized `CONTCAR` files into MD-stage `POSCAR` files, and then submits MD jobs.

The default path configuration is:

```python
OPT_INPUT_DIR = "inopt"
MD_INPUT_DIR = "inmd"
POSCAR_DIR = "poscar"
OPT_ROOT = "opt"
STD_ROOT = "std"
```

Before running this step, make sure the generated structures are in `poscar/`:

```bash
mv new_POSCAR_folder poscar
python 02_data_gen.py
```

Alternatively, change `POSCAR_DIR` in `02_data_gen.py` to `new_POSCAR_folder`.

This script requires the following LSF commands:

- `bsub`
- `bjobs`
- `bkill`

During optimization, the script monitors energy changes in `OSZICAR`. During MD, it monitors the first C-O distance and stops a job if the distance leaves the default range of `1.3-3.8 Angstrom`. The MD summary is written to:

```text
std/all_md_steps.csv
```

### 3. Collect OUTCAR files

Run the collector from inside `std/`:

```bash
cd std
python 03_data_gen.py
```

The script scans calculation subdirectories and copies each `OUTCAR` file into:

```text
std/OUfolder/OUTCAR_<number>
```

### 4. Extract C-O distance, energy, and force projection

Run the analysis script from inside `std/OUfolder/`:

```bash
cd OUfolder
python 04_data_analysis.py
```

`04_data_analysis.py` uses ASE to read VASP `OUTCAR` trajectories and extracts the following values for each ionic step:

- C-O distance
- potential energy
- force on the C atom projected along the C-O direction
- force on the O atom projected along the C-O direction
- `force_o - force_c`

Default output:

```text
sorted_results1.csv
```

### 5. Scan confidence interval windows

`05_data_analysis.py` scans different confidence interval (CI) windows, compares adjacent CI-averaged force curves, and reports convergence/stability metrics.

```bash
python 05_data_analysis.py
```

Output directory:

```text
CI_diff_results/
```

Main outputs:

- `delta_force_diff_scan_lower_long.csv`
- `delta_force_diff_scan_lower_wide.csv`
- `delta_force_diff_scan_upper_long.csv`
- `delta_force_diff_scan_upper_wide.csv`
- `metrics_lower_adjacent.csv`
- `metrics_upper_adjacent.csv`
- `plateau_lower_summary.csv`
- `plateau_upper_summary.csv`

### 6. Generate force and free energy figures

`06_data_figures.py` filters `force_diff` using selected CI windows, plots the force data, integrates the CI-averaged force curve into a relative free energy curve, and also produces binned force/free energy profiles.

```bash
python 06_data_figures.py
```

The default CI setting is defined in the script:

```python
ci_settings = [
    (0.25, 0.95),
]
```

Example output directory:

```text
CI_0.25_0.95/
```

Main outputs:

- `force_scatter_and_CI_filtered_points.png`
- `CI_mean_free_energy_curve.png`
- `Binned_CI_mean_force_curve.png`
- `Binned_CI_mean_free_energy_curve.png`
- `ci_filtered_force_free_energy.csv`
- `binned_force_free_energy_CI.csv`

## Dependencies

Python dependencies for post-processing:

```bash
pip install numpy pandas matplotlib ase tqdm "scipy<1.14"
```

Notes:

- `02_data_gen.py` requires VASP, an LSF scheduler, and a cluster-specific `vasp.script`.
- `04_data_analysis.py` relies on ASE to parse VASP `OUTCAR` files.
- `06_data_figures.py` uses `scipy.integrate.cumtrapz`. Newer SciPy versions may remove this interface, so `scipy<1.14` is recommended unless the script is updated to use `cumulative_trapezoid`.

## Important Notes

1. `02_data_gen.py` reads structures from `poscar/` by default, while `01_data_gen.py` writes to `new_POSCAR_folder/`. Make sure these directory names are aligned before submitting jobs.
2. `04_data_analysis.py`, `05_data_analysis.py`, and `06_data_figures.py` use `sorted_results1.csv` as the shared intermediate data file.
3. The scripts use the first C atom and the first O atom in the structure as the C-O pair for distance and force projection analysis. If your system contains multiple C/O atoms, verify that this matches the intended reaction coordinate.
4. `POTCAR` is subject to the VASP license and is intentionally not included in this repository. Users must provide their own licensed pseudopotential files.
5. VASP output files such as `WAVECAR`, `CHGCAR`, `OUTCAR`, and `vasprun.xml` can be very large. For public release, keep only small example data in the repository and consider using Git LFS, GitHub Releases, or an external data repository for full calculation outputs.

## Citation

If this workflow is used in a publication or report, add the appropriate citations for VASP, ASE, and this project.

## License

No open-source license has been specified yet. If you plan to publish this repository and allow reuse, add a suitable `LICENSE` file before release.
