# HeavyMesonLifetimes
Matrix Elements for Heavy Meson Lifetimes from Lattice QCD

Software source code based on [Grid](https://github.com/paboyle/Grid) and [Hadrons](https://github.com/aportelli/Hadrons), which are both free software under GPLv2. Both Grid & Hadrons need to be installed as dependencies, and production binaries were using Grid up to commit [5a4f9bf2](https://github.com/paboyle/Grid/tree/5a4f9bf2e35787e39e4f87c37d8acd0c56fa49c9) and a fork of Hadrons up to [579427b](https://github.com/mbr-phys/Hadrons/tree/feature/GradientFlow); this work has been submitted as a [PR to the main Hadrons branch](https://github.com/aportelli/Hadrons/pull/137).

Round 2 (HQE operators, `HQE_plan.md`) needs the gradient-flow Hadrons of branch `feature/HQE` (fork [mbr-phys/Hadrons](https://github.com/mbr-phys/Hadrons)), which provides `MGradientFlow::FermionFlow`, `MContraction::HeavyMeson3pt` and `MContraction::WeakNonEye3ptDeriv`.
It builds two executables: `FermionFlow_Lifetimes` (everything; with `runPar.eyes = 1` also the dim-6 $\Delta B=0$ eyes and the stochastic condensates) and `FermionFlow_Eyes` (eyes and condensates only, from wall propagators saved by `FermionFlow_Lifetimes` with `runPar.saveWalls = 1`, using the same input file and `runId`). See `HQE_plan.md` §7.2.

The code can be compiled using the sequence of commands below
```
./bootstrap.sh
mkdir build
cd build
../configure --with-grid=<Grid prefix> --with-hadrons=<Hadrons prefix> CXX=<compiler used for Grid>
make
```
