# n3lo-structure-function-benchmarks

This repository contains the C++ files that were used to generate all the benchmark results in

* Valerio Bertone and Alexander Karlberg: arXiv:24XX.XXXXX

In particular the files benchmark the deep inelastic scattering structure functions at $\mathcal{O}(\alpha_S^3)$ as implemented in <span style="font-variant:small-caps;">apfel</span>++ and <span style="font-variant:small-caps;">Hoppet</span>.

In order to compile the code both programs must be installed

* <span style="font-variant:small-caps;">apfel</span>++ (https://github.com/vbertone/apfelxx)
* <span style="font-variant:small-caps;">Hoppet</span> (https://github.com/hoppet-code/hoppet)

# Installation
To compile the benchmarks simply type

    make [-j]

There are three executables that can be run. The all use as the PDF initial condition section 1.3 of http://arxiv.org/abs/hep-ph/0204316. Typing

    ./TabulateStructureFunctions

will print a comparison of <span style="font-variant:small-caps;">apfel</span>++ and <span style="font-variant:small-caps;">Hoppet</span> on the screen for Q =  91.1876 GeV. The benchmark shows a few PDF combinations along with structure functions at N3LO.

To reproduce the tables from arXiv:24XX.XXXXX one can run

    ./StructureFunctionsJoint

This will produce two files in addition to a print on the screen that should show that the ratio between the output of the two evolution codes is 1. `table_N3LO.tex` contains the tables that are printed in the paper. `StructureFunctions_N3LO.dat` contains all the structure functions and reduced cross sections at N3LO for a number of Q-values and finely spaced in x.

Finely one can run 

    ./ScaleVariations

This will produce a file with the structure functions evaluated at multiple values of the renormalisation and factorisation scale. The output can be found in `F2NC_Scale_Variations_N1LO.dat`, `F2NC_Scale_Variations_N2LO.dat`, and `F2NC_Scale_Variations_N3LO.dat`.