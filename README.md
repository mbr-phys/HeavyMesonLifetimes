# HeavyMesonLifetimes
Matrix Elements for Heavy Meson Lifetimes from Lattice QCD

Software source code based on [Grid](https://github.com/paboyle/Grid) and [Hadrons](https://github.com/aportelli/Hadrons), which are both free software under GPLv2. Both Grid & Hadrons need to be installed as dependencies, and production binaries were using Grid up to commit [5a4f9bf2](https://github.com/paboyle/Grid/tree/5a4f9bf2e35787e39e4f87c37d8acd0c56fa49c9) and a fork of Hadrons up to [579427b](https://github.com/mbr-phys/Hadrons/tree/feature/GradientFlow); this work has been submitted as a [PR to the main Hadrons branch](https://github.com/aportelli/Hadrons/pull/137).

Round 2 (HQE operators, `HQE_plan.md`) needs the gradient-flow Hadrons of branch `feature/HQE` (fork [mbr-phys/Hadrons](https://github.com/mbr-phys/Hadrons)), which provides `MGradientFlow::FermionFlow`, `MContraction::HeavyMeson3pt` and `MContraction::WeakNonEye3ptDeriv`.
It builds two executables: `FermionFlow_Lifetimes` (everything; with `runPar.eyes = 1` also the dim-6 $\Delta B=0$ eyes and the stochastic condensates) and `FermionFlow_Eyes` (eyes and condensates only, from wall propagators and wall sources saved by `FermionFlow_Lifetimes` with `runPar.saveWalls = 1`; its `runId` only seeds the loop noise, so a new `runId` gives new noise hits). See `HQE_plan.md` §7.2. With `hqePar.kaonMixing = true` (needs `mesonPar.light = 1`) it also writes the dim-6 $\Delta S=2$ kaon-mixing 3pts `ls_DeltaHs2_…`, the `DeltaHs2` contraction with the strange quark in the heavy role. All 3pts (connected four-quark, eyes, two-quark) have the heavy quark running from the wall at `t_i` to the wall at `t_i + dT`, with `adj(Gi)` at `t_i` and `Gf` at `t_f` for a result tagged `_<Gi>_<Gf>_`; the contraction modules write raw traces without Wick signs (see `four-quark-combination.tex` and `HQE_Hadrons.tex` in the gradflow notes).
Example inputs: `run/input_fermion.xml` (`FermionFlow_Lifetimes`, saving the walls) and `run/input_eyes.xml` (`FermionFlow_Eyes`). The eyes input must repeat the saving run's `configPar`, quark names, `flowPar`, `mesonPar.light` and `runPar.wallStem`/`wallIo32`. Using the same file for both also works.

`flowPar.heavyLinks = stout` flows the heavy walls from the stout links of the heavy action (`Vgauge`, a second flow chain whose modules and gauge observables are `WilsonFlowV_t<τ>`), and the two-quark vertices use those links, at τ = 0 too; light/strange walls, eye noise, four-quark modules and eyes stay on the thin links, except that the dim-7 derivatives on the heavy legs use the stout chain (`WeakNonEye3ptDeriv.legGauges`, Hadrons `feature/HQE`); total-derivative relations between the dim-7 specs then hold only up to O(a²). Default `thin`: one chain from `trans_gauge`. `FermionFlow_Eyes` then also reads `stoutPar`, which must match the saving run.

Propagators from an earlier run (e.g. `HeavyMassScan` of SmearingSearch with `savePar.save = 1`) can be read in instead of solved, flavour by flavour: `hdwfPar.prop`, `sdwfPar.prop` and `ldwfPar.prop` are `none` (solve, or in `FermionFlow_Eyes` read `<wallStem>/<quark>_<t>`) or a file stem, and the propagator of wall `t` is then read from `<prop>_<t>.<traj>.bin`, e.g. `props/b0.45` for the heavy mass 0.45 of HeavyMassScan. Loaded propagators have no Ward identity check (no 5D solution). With `runPar.loadSrcs = 1` the wall sources are read from `runPar.srcStem` too (`src_<t>`, and `z2_<t>` if `smearPar.smear = 1`), so solved and loaded propagators share the same sources whatever the runId; with `loadSrcs = 0` they are generated, which needs the runId and smearing of the saving run. `runPar.propIo32` gives the precision of both. The loading run must repeat the saving run's `configPar` (gauge files, trajectories, `source_loc`, sources), quark masses and actions, `stoutPar`, and the light/strange smearing (HeavyMassScan on thin links, `sourcePar.jacobi.links = thin`); none of this is checked.

The code can be compiled using the sequence of commands below
```
./bootstrap.sh
mkdir build
cd build
../configure --with-grid=<Grid prefix> --with-hadrons=<Hadrons prefix> CXX=<compiler used for Grid>
make
```
