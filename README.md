# InfiniteLIOMs

This repository contains data and Julia scripts for the paper "Algorithm for finding local integrals of motion in quantum lattice models in the thermodynamic limit." The scripts compute local integrals of motion (LIOMs) for 1D XXZ model, using two types of basis operators: Pauli strings basis and symmetry resolved basis. See the text for details and definitions.

## Overview of Scripts

### Prerequisites
Julia environment can be initialized via the `Project.toml` file, containing information about dependencies. In the `scripts` directory, run

```bash
julia --project=.
```
and then

```julia
] instantiate
```
to set everything up.

### 1. `lioms_trans_sym_pauli.jl`
This script computes M-local LIOMs in the Pauli strings basis, for the XXZ model.

#### Usage:
Run the script with the following command-line arguments:
- `--delta` (`-d`): Anisotropy parameter Δ (default: 0.5).
- `--max-supp` (`-M`): Maximum support M for LIOMs (default: 3).

#### Example:
```julia
> julia -t10 --project=. lioms_trans_sym_pauli.jl -d 0.5 -M 2
Generating operators...
Generated 12 operators
Computing norms...
Precomputing commutators...
10 smallest eigenvalues:
[0.0, 1.0658141036401503e-14, 2.0, 2.0, 2.000000000000032, 7.999999999999998, 8.000000000000007, 17.999999999999996, 18.0, 18.0]
10 eigenvectors corrsponding to smallest eigenvalues:
∑_l X1     0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
∑_l Z1     1.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
∑_l Y1     0.00000  0.00000  1.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
∑_l XX     0.00000  0.66667  0.00000  0.00000 -0.00000  0.00000 -0.00000  0.11171  0.43304  0.59628
∑_l ZX     0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.70711  0.00000  0.00000  0.00000
∑_l YX     0.00000  0.00000  0.00000  0.00000  0.70711  0.00000 -0.00000 -0.68469  0.17662  0.00000
∑_l XZ     0.00000  0.00000  0.00000  0.00000  0.00000  0.00000 -0.70711  0.00000 -0.00000  0.00000
∑_l ZZ     0.00000  0.33333  0.00000  0.00000  0.00000 -0.00000  0.00000 -0.22341 -0.86608  0.29814
∑_l YZ     0.00000  0.00000  0.00000  0.00000 -0.00000 -0.70711 -0.00000  0.00000  0.00000  0.00000
∑_l XY     0.00000  0.00000  0.00000  0.00000 -0.70711 -0.00000 -0.00000 -0.68469  0.17662  0.00000
∑_l ZY     0.00000  0.00000  0.00000  0.00000 -0.00000  0.70711 -0.00000  0.00000  0.00000  0.00000
∑_l YY     0.00000  0.66667  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000 -0.74536
Size of basis for M=2 : 12
Number of LIOMs found: 2
Elapsed time: 0.762079875 s
```
In this case, the script correctly reproduces the total $S^z$ and the Hamiltonian of XXZ model.

### 1. `lioms_trans_sym.jl`
This script computes LIOMs in a symmetry-resolved basis, for the XXZ model. It allows for constraints based on time-reversal symmetry, parity, and conservation of total $S^z$.

#### Usage:
Run the script with the following command-line arguments:
- `--delta` (`-d`): Anisotropy parameter Δ (default: 0.5).
- `--max-supp` (`-M`): Maximum support M for LIOMs (default: 3).
- `--time-reversal` (`-T`): Behavior under time-reversal symmetry (`even`, `odd`, `both`; default: `odd`).
- `--parity` (`-P`): Behavior under parity (`even`, `odd`, `both`; default: `even`).
- `--conserve-Sz` (`-C`): Conserve total Sz? (`yes`, `no`, `both`; default: `yes`).

#### Example:
```julia
> julia -t10 --project=. lioms_trans_sym.jl -d 0.5 -M 3 -T odd -P even -C yes
Generating operators...
Generated 3 operators
Computing norms...
Precomputing commutators...
3 smallest eigenvalues:
[-1.1814295211437466e-16, 0.37499999999999983, 3.9999999999999942]
3 eigenvectors corrsponding to smallest eigenvalues:
∑_l i(321 - 123)    0.81650  0.57735 -0.00000
∑_l i(231 - 213)   -0.40825  0.57735 -0.70711
∑_l i(312 - 132)   -0.40825  0.57735  0.70711
Size of basis for M=3, restricted to time reversal = odd, parity = even, Sz conservation = yes: 3
Number of LIOMs found: 1
Elapsed time: 0.994267708 s
```
The following encoding is used: $\mathrm{I} \to 0,\, S^+ \to 1,\, S^z \to 2,\, S^- \to 3$.
For the above parameters, one LIOM is produced, corresponding to the energy current $Q_3$.

### Output 
Both scripts save results in the `results/` directory under subfolders for each basis:
- `results/1D_XXZ/symmetry_resolved/` for `lioms_trans_sym.jl`.
- `results/1D_XXZ/Pauli_strings/` for `lioms_trans_sym_pauli.jl`.

The outputs include:
- `eigenvalues_*.txt`: Eigenvalues of the $F$ matrix
- `eigenvectors_*.txt`: Eigenvectors of the $F$ matrix.
- `operators_*.txt`: Basis operatrs.
- `log_*.txt`: Log file with computation parameters and summary.

## Overview of data
All data shown in paper is hosted on [Zenodo](https://doi.org/10.5281/zenodo.15363681). In `plots/` there are Python scripts creating the plots, assuming that the data is downloaded and unpacked in `data/`. If you are interested in structure of the data files, please consult the plotting scripts for details.


## Acknowledgments
The algebra of operators in these scripts is implemented via the library [`PauliStrings.jl`](https://github.com/nicolasloizeau/PauliStrings.jl) developed by Nicolas Loizeau, J. Clayton Peacock and Dries Sels and published in: Loizeau, J. C. Peacock, and D. Sels, Quantum many-body simulations with PauliStrings.jl, [SciPost Physics Codebases , 054 (2025)](https://scipost.org/10.21468/SciPostPhysCodeb.54).

Marcin Mierzejewski acknowledges support by the National Science Centre (NCN), Poland via project 2020/37/B/ST3/00020. Jakub Pawłowski acknowledges support by the National Science Centre (NCN), Poland via project 2023/49/N/ST3/01033. 

## Contact
For question, please contact `jakub.pawlowski@pwr.edu.pl`.
