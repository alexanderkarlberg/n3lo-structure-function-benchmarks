# n3lo-structure-function-benchmarks

This repository contains the C++ codes used to benchmark the inclusive deep-inelastic-scattering (DIS) structure functions at $\mathcal{O}(\alpha_S^3)$, _i.e._ N<sup>3</sup>LO, as implemented in `apfel++` and `Hoppet` and presented in:

* V. Bertone and A. Karlberg, _Benchmark of deep-inelastic-scattering structure functions at_ $\mathcal{O}(\alpha_S^3)$, arXiv:24XX.XXXXX,

along with the suite of `Matplotlib` scripts used to produce the plots shown in that paper and many others.

In order to compile the code, both `apfel++` and `Hoppet` must be installed. See:

* [`apfel++`](https://github.com/vbertone/apfelxx),
* [`Hoppet`](https://github.com/hoppet-code/hoppet),

for the respective installation instructions.

# Installation

To compile the benchmark codes, simply go into the [`code`](code/) folder and type:

    make [-j]

This will produce three executables. They all use the same PDF initial condition of section 1.3 in http://arxiv.org/abs/hep-ph/0204316.

# Execution

In the [`code`](code/) folder, running the code:

    ./StructureFunctionsJoint

will print on the screen comparison tables between `apfel++` and `Hoppet` for different values of $x_{\rm B}$ and $Q$ and for all perturbative orders from LO to N<sup>3</sup>LO. On top of the structure functions, the tables also show a few PDF combinations. Results on screen are presented as |`apfel++` / `Hoppet` - 1|. This code will also dump into the [`plots`](plots/) folder the files: `StructureFunctions_N0LO.dat`, `StructureFunctions_N1LO.dat`, `StructureFunctions_N2LO.dat`, and `StructureFunctions_N3LO.dat` that can be used to produce a number of plots (see below). Moreover, also the files: `table_N0LO_APFEL.tex`, `table_N1LO_APFEL.tex`, `table_N2LO_APFEL.tex`, `table_N3LO_APFEL.tex`, `table_N0LO_HOPPET.tex`, `table_N1LO_HOPPET.tex`, `table_N2LO_HOPPET.tex`, and `table_N3LO_HOPPET.tex` will be created in the [`tables`](\tables) folder.

Similarly to the code above, running the code:

    ./TabulateStructureFunctions

will print on screen comparison tables for both `apfel++` and `Hoppet` this time showing the values of structure functions and, again, few PDF combinations.

Finally, running the code: 

    ./ScaleVariations

will produce files in the (`plots`)[plots/] folder with the structure function $F_2^{\rm NC}$ evaluated at multiple values of the renormalisation and factorisation scales $\mu_{\rm R}$ and $\mu_{\rm F}$, respectively. The files produced are: `F2NC_Scale_Variations_N0LO.dat`, `F2NC_Scale_Variations_N1LO.dat`, `F2NC_Scale_Variations_N2LO.dat`, and `F2NC_Scale_Variations_N3LO.dat` and can be used to produce plots.

# Producing the plots

The [`plots`](plots/) folder, where some of the codes described above dump part of their output, also contains some `python` scripts that utilise `matplotlib` to produce plots in `.pdf` format. The scripts are:

* `PerturbativeConvergence.py`
* `ScaleVariations.py`
* `StructureFunctions_N0LO.py`
* `StructureFunctions_N1LO.py`
* `StructureFunctions_N2LO.py`
* `StructureFunctions_N3LO.py`

and are responsible for all plots contained in the [`plots`](plots/) folder.

# Tables

The [`tables`](tables/) folder contains the benchmark tables presented in the paper in `LaTeX` format. The file `tables.tex` collects the table file produced by the `StructureFunctionsJoint` code and can be compiled with `LaTeX` to produce [this](tables/tables.pdf) file.

# Contacts

- Valerio Bertone: valerio.bertone@cern.ch
- Alexander Karlberg: alexander.karlberg@cern.ch
