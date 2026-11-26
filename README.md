# Hybrid ion–atomic-ensemble entanglement distribution protocol comparison

Code accompanying the paper/preprint

> "Comparison of long-distance entanglement generation protocols with trapped ions and multiplexed atomic ensembles,"

The code computes the figures of merit for the protocols studied in the paper. These are the average entanglement distribution time $T_{\rm ent}$, the worst-case memory and ion waiting times, the crosstalk figure $N_{\rm ion}$, and the optimal photon probabilities. It covers two ways of bridging the bandwidth mismatch between trapped-ion and atomic-ensemble photons:

- **Memory matching (M-M):** the ensemble memory's readout is shaped to match the ion photon.
- **SPDC matching (SPDC-M):** a time-binned SPDC source at the edge node is matched to the ion photon.

Each scheme is covered for single-click (SC) and double-click (DC) protocols, with and without a central memory repeater, and with or without turning off the ion (and SPDC) drive after a click.

---

## Contents

```
.
├── Protocol_Fns.jl              # Memory_Match, SPDC_Match and helper functions
├── Distance_Plots.ipynb         # T_ent, waiting times and N_ion vs distance
├── Distance_Plots_Cap.ipynb     # same as Distance_Plots, with a cap on N_ion for SPDC-M
├── Dephasing_Plots.ipynb        # T_ent and optimal P_i vs dephasing parameter r
├── Memory_Plots.ipynb           # T_ent vs memory efficiency η_m
├── MM_vs_SPDCM.ipynb            # direct M-M vs SPDC-M comparison vs η_AFC and T_EN
├── LICENSE                      
└── README.md
```


## Requirements

- Julia ≥ 1.12. 
- Julia packages:
  - [`Optim.jl`](https://github.com/JuliaNLSolvers/Optim.jl), needed by `Protocol_Fns.jl`.
  - [`CairoMakie.jl`](https://docs.makie.org/) and [`LaTeXStrings.jl`](https://github.com/JuliaStrings/LaTeXStrings.jl), needed only for the notebooks.
  - [`IJulia.jl`](https://github.com/JuliaLang/IJulia.jl), to run the notebooks in Jupyter.

### Installation

```julia
using Pkg
Pkg.add(["Optim", "CairoMakie", "LaTeXStrings", "IJulia"])
```

## Usage

I recommend reading the function documentation in `Protocol_Fns.jl` (there is one function for each matching method, `Memory_Match` and `SPDC_Match`), and then looking at the attached notebooks (which mostly create figures from the paper), and going from there. Each one sets up the parameters, calls the functions and plots the results for one or more figures in the paper. One can modify the input parameters and create similar graphs, or use the functions directly to calculate whatever one wants.

| Notebook | Figures |
|---|---|
| `Distance_Plots.ipynb` | Fig. 6, Fig. 10, Figs. S2, S3, S6 |
| `Dephasing_Plots.ipynb` | Fig. 7, Fig. S4 |
| `Memory_Plots.ipynb` | Fig. 8 |
| `MM_vs_SPDCM.ipynb` | Fig. 9, Fig. S5 |
| `Distance_Plots_Cap.ipynb` |  n/a  |

Some things to note:
In general the code does not check if the optimized photon probabilities are sufficiently small such that the perturbation theory assumed in the paper still holds; generally this will be true for reasonable inputs, but you should check yourself that they are sufficiently small (to whatever you deem as a valid regime for the various first order approximations to hold). The notebooks do this by plotting all (or most in some cases) of the various photon probabilities.

The optimizer is not guaranteed to work every time of course, but it seems to be quite robust for all regimes studied in this work. One can of course tweak the optimization procedure in the code functions if it is deemed necessary. One can also avoid optimization and just input specific photon probabilities for both of the matching functions.

We have tried to make the error messages decriptive if some parameter exceeds a physical regime; this may not be perfect, and you may need to look into the functions a bit if this happens.

## Citation

If you use this code, please cite our paper.

## License

This code is released under the MIT License; see LICENSE for details.

## Acknowledgments

We thank T. E. Northup, H. de Riedmatten, S. D. C.
Wehner, M. van Hooft, N. Sangouard, B. Lanyon, S. Maiti, V. Krutyanskiy, and M. Teller for fruitful discussions. 
This work was funded by the European Union's Horizon Europe research and innovation programme under grant agreement No. 101102140 -- QIA Phase 1 and No. 101297140 -- QIA Phase 2. 

Funded by the European Union. Views and opinions expressed are however those of the authors only and do not necessarily reflect those of the European Union or European Commission. Neither the European Union nor the granting authority can be held responsible for them.

## Contact

Chris Gustin; chris.gustin@nbi.ku.dk