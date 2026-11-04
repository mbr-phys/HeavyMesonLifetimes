# HeavyMesonLifetimes round 2: HQE operators and dim-7 mixing

This is the plan for updating `FermionFlow_Lifetimes` so that one run produces everything needed for:
- the HQE matrix elements $\langle\bar QQ\rangle$, $\mu_\pi^2$, $\mu_G^2$, $\rho_D^3$, $\rho_{LS}^3$, with their remainders and full-$D$ forms;
- the dim-6 $\Delta B=0,2$ four-quark operators (as in round 1);
- the dim-7 $\Delta B=2$ operators $R_0,R_1,\tilde R_1,R_2,\tilde R_2,R_3,\tilde R_3$.

Everything is done at every measured flow time $\tau$, including $\tau=0$.

**Sources.**
- Physics and operator definitions: `Hadrons-NewGF/lifetimes_round2.md` (§2 HQE, §4 dim-7) and `~/library/LQCD/gradflow/notes/hqe_notes.tex`.
- Hadrons side: `Hadrons-NewGF/newCode.md`, which documents `MContraction::HeavyMeson3pt` (two-quark) and
  `MContraction::WeakNonEye3ptDeriv` (four-quark with gamma pairs and derivatives) on branch `feature/HQE`.
- App structure: `../SParameterGF`, which started from this app and has since been updated to the current gradient-flow
  Hadrons (§1.1). P0 follows its patterns wherever they apply.

**Status (2026-10-29): P0, P1 and P2 code written, not compiled or run, uncommitted.** The dim-7 operators now use
`MContraction::WeakNonEye3ptDeriv` instead of a four-quark part of `HeavyMeson3pt` (§9). Changed: `FermionFlow_Lifetimes.cpp` (rewritten),
`Functions.hpp`, `Inputs.hpp`, `run/input_fermion.xml`, `README.md`. P3: dim-6 eyes and condensates written (2026-11-02, item 15,
§7.2), the rest not started. See the step list below and §9.

---

## Implementation overview

Steps in order. Each phase ends with its validation before the next starts; details are in the sections given.
✅ = code written (not compiled), ⏳ = waiting for a build/run, no mark = not started.

**Prerequisite.** Build and test `feature/HQE` Hadrons (`./bootstrap.sh`, build, `Test_heavy_meson_3pt`; see `Hadrons-NewGF/newCode.md`).

**P0: port to the current gradient-flow Hadrons (no new physics)**

1. ✅ `Inputs.hpp`: add `runPar.flowChaining` (`step` | `interval`) (§4.1).
2. ✅ `Functions.hpp`: add `make_fermionflow`, `measuredSteps`, `tauLabel` (§2.1, §4.2).
3. ✅ `FermionFlow_Lifetimes.cpp`:
   - per-wall struct `{t, b, s, l, z2, gsrc}` and wall-pair lookup $t_f=(t_i+\Delta t)\bmod T$ (§2.2, bugs 3–4);
   - `WilsonFermionFlow` → chained `FermionFlow` in both chaining modes, with the `SParameterGF` measurement schedule (§2.1, §4.3);
   - fix bugs 1, 2 and 5, and unify the 2pt names (§2.2);
   - leave the dim-6 `WeakNonEye3pt` calls as they are for now.
4. ⏳ Validate (needs a build): regression against the round-1 binary on one configuration (§6.1); chaining benchmark `step` vs `interval` (§6.6).

**P1: probe gammas, dim 6 and the two-quark gamma vertices**

5. ✅ `Inputs.hpp`: add the `HQEPar` block and its input checks (§4.1).
6. ✅ Dim-6 `WeakNonEye3pt` for every (`gammaIn`, `gammaOut`) in `probeGammas`, with round-1-style names (§3.2, §4.3).
7. ✅ `HeavyMeson3pt` two-quark instance `<q>h_HQE_…` with the `all` vertex only (§3.1, §4.2–4.3).
8. ⏳ Validate (needs a build): $\gamma^0$ normalisation plateau and $\langle\bar QQ\rangle$ (first items of §6.3).

**P2: HQE operators and dim-7 $\Delta B=2$**

9. ✅ Add the `ops` vertices to the two-quark instance (§3.1).
10. ✅ `WeakNonEye3ptDeriv` dim-7 modules `<q>h_DeltaHs2_dim7_<G>_<G>_…` (§3.2, §4.2–4.3):
    - derivative specs from `dim7Specs()` (Lenz–Nierste);
    - one module per probe pair.
11. ✅ Example input with `<hqePar>` and `flowChaining` (§4.4).
12. ⏳ Validate (needs a build):
    - dim-7 `none` entries equal `WeakNonEye3pt` (§6.2);
    - HQE sanity checks (§6.3);
    - one-off parity check with `gammaPairs = all` (§6.4);
    - timings (§6.5).

**P3 (optional, later; details in §7)**

Only once P0–P2 run in production and the §6.5 timings show where the time goes.

13. **Share the heavy-quark work between spectators** (Hadrons, `HeavyMeson3pt`).
    - Today each spectator ($s$, $l$) needs its own two-quark instance, so the heavy-propagator derivatives and $F_{\mu\nu}$
      insertions for a wall pair are computed twice.
    - Letting `spectator`/`sinkNoise` be lists would do that work once per wall pair and τ.
14. **Stop keeping the unflowed spectators in memory** (Hadrons, `MSink::SourceSink` + `HeavyMeson3pt`).
    - The two-quark 3pts need the unflowed $s$/$l$ propagators at every τ, but only through one time-slice sum
      ($q_1^{\rm snk}$, the noise-projected spectator at the sink).
    - Computing that sum once at τ = 0 (`SourceSink` with `sourceOp = adj`, now available) and passing it to `HeavyMeson3pt` would let
      the full spectator propagators be freed.
    - Saves $n_{\rm walls}\times n_{\rm spec}$ propagators held for the whole run (§5).
15. ✅ **Dim-6 $\Delta B=0$ eye diagrams and condensates** (§7.1, §7.2; one Hadrons change: `MSink::SourceSink` `sourceOp`).
    - New stochastic inversions on diluted Z2 noise for $s$/$l$; this is the only P3 item with extra inversions.
    - Noise and solutions flowed in the same `FermionFlow` as the propagators, then `Loop` + `WeakEye3pt` with the flowed heavy
      legs and the unflowed spectator, plus the `SParameterGF` condensates on the same fields.
    - Two ways to run, to be benchmarked (§7.2): in the main run (`runPar.eyes = 1`), or in `FermionFlow_Eyes` from wall propagators
      saved by the main run (`runPar.saveWalls = 1`).
    - ⏳ Validate (needs a build): the two ways must agree (§7.2), and benchmark them.
16. **Other probes** (app only).
    - $B^*$ HQE parameters via vector probe gammas (`GammaX/Y/Z`) in the triads and 2pts, if needed later.

Far future, not planned: $\Delta B=0$ dim-7 (connected and eye) and the Darwin four-quark cross-check, which need new Hadrons eye code.

**Still open** (§8): gauge field for the τ = 0 operators (`trans_gauge` for now, §8.1); probes $P$ only or $P$ and $A_0$ (§8.4);
production chaining mode, settled by the benchmark in step 4 (§8.7).

---

## 0. Summary

| Step | Change | New Hadrons code needed? |
|---|---|---|
| P0 | Port to the `feature/HQE` Hadrons following `SParameterGF`: `WilsonFermionFlow` → `FermionFlow`, flow-time bookkeeping, input validation; known bugs | no |
| P1 | Keep the dim-6 `WeakNonEye3pt` calls (new lattices), extended to the probe gammas; add the two-quark gamma-vertex 3pts (`HeavyMeson3pt`) | no |
| P2 | Add the HQE operator vertices (`ops`) and the dim-7 $\Delta B=2$ `WeakNonEye3ptDeriv` modules (Lenz–Nierste specs); XML switches | no |
| P3 (opt.) | Hadrons-side efficiency (several spectators per instance, sliced spectator); stochastic part for dim-6 $\Delta B=0$ eyes (app-side, modelled on `SParameterGF`'s `FermionFlow_Condensates`, §7.1). Far future: $\Delta B=0$ dim-7 | efficiency items only |

P0–P2 need no further Hadrons changes; `feature/HQE` already has everything they use. Each step is validated before the next (§6).

---

## 1. Current app (round 1), for reference

- **Gauge field.** `gauge` (NERSC) → `trans_gauge` (`CShiftGaugeField` to `source_loc`) → `Vgauge` (stout smearing, used **only** by the
  heavy Möbius action), plus `trans_gaugeF` (single precision, used by the light/strange mixed solver).
- **Walls.** One Z2 wall per entry of `sources`: `z2_p000_t<t>`, Jacobi-smeared to `sm_z2_p000_t<t>_<t>` if `smear == 1`. Propagators:
  - heavy `b_<t>` (point Z2);
  - strange `s_<t>` and, if `light == 1`, light `l_<t>` (from the smeared Z2 if `smear == 1`).
- **Flow.** One `WilsonFermionFlow` module **per step**. Flowed propagators are `<q>_t<τ>` for measured steps.
- **Per wall and τ:**
  - 2pts `Meson`: sh, ss, hh (+ lh, ls, ll), with the `gamma_pairs` from `mesonPar`;
  - smeared–smeared 2pts at τ = 0 only;
  - `WardIdentity` at τ = 0.
- **Per wall pair and τ.** Pairs are formed as $(t_i,t_f=t_i+\Delta t \bmod T)$ for each `deltaTs`. Three `WeakNonEye3pt` (G5, G5) calls each for sh and lh:
  - `DeltaHs2` with legs $(h_i, s_i, s_f, h_f)$;
  - `DeltaHs0_T1` with $(s_i, h_i, s_f, h_f)$;
  - `DeltaHs0_T2` with $(s_i, s_i, h_f, h_f)$.
- **Output.** One `ResultGroup` per τ, all written to one file by `WriteResultGroup`.

### 1.1 What `../SParameterGF` already does (reviewed 2026-10-28)

`SParameterGF` was built from this app and updated for the current gradient-flow Hadrons. Its main driver now targets
CLS open-boundary ensembles (Wilson-clover); a periodic DWF variant exists but is disabled. What carries over:

| `SParameterGF` pattern | Where | Use here |
|---|---|---|
| `FermionFlow` chained **one step per module**: `steps = 1`, `meas_interval = 1`, `defaultType = "PropagatorField"`; step $t$ takes `props` = the previous step's `outProps` and `gauge` = previous module + `"_U"`; `outProps` = `q + "_t<τ>"` always set | `FermionFlow_SParameter.cpp`, `FermionFlow_Condensates.cpp` | P0 flow chain (§2.1) |
| Measurement schedule on step indices: `coarsen_step = ceil(Tcoarsen/step_size − 1e-12)`, interval ×4 from there, `last_contraction_step`, always measure the last step | same | replaces the round-1 `ti == Tcoarsen` double compare (bug 2) |
| Flow module (`gauge_obs`) put in every τ's `ResultGroup`, measured or not; `WilsonFlow` with `steps = 0` for τ = 0 | same | unchanged from round 1 |
| Consistent 2pt names `<flav>_LL_t<t>_<τ>` at τ = 0 and τ > 0 | `FermionFlow_SParameter.cpp` | fixes bug 6 |
| Input checks that fail the run (`gamma_src`/`gamma_snk` length mismatch, `nSrcs`/`nsparse` > 0; `SParameterGF` also rejects an empty `gamma_snk`) | same | same checks for the new inputs (§4.1); an empty `gamma_snk` keeps meaning `all` here (§2.2) |
| `deltaTs`, `stp`, `Tmax` index arithmetic removed (commit "remove unnecessary deltaTs") | git history | here the wall pairing still needs Δt: use the lookup of bug 3 |
| Stochastic sidecar: `SparseSpinColorDiagonal → Z2Diluted → PropagatorVectorUnpack`, solve per component, flow **noise and solution together in one batched `FermionFlow`**, then `StochasticCondensate` and `DslashField` | `FermionFlow_Condensates.cpp`, helpers in `Functions.hpp` | P3: eyes (`Loop` + `WeakEye3pt` on the same flowed pairs, §7.1); the condensates themselves are not needed here |
| Helpers `make_stochastic_noise_diluted`, `make_stochastic_condensate`, `make_dslash_field`, `make_stochastic_meson` | `Functions.hpp` | copy when P3 starts; not needed for P0–P2 |

What does **not** carry over:
- **Open boundaries.** `LoadOpenQcd`, `WilsonFlow_OBC`, `FermionFlow_OBC` (`cG`), fermion `bc = 0` and `ContractConservedSink` are on
  `feature/contractConservedCurrents`, not on `test/AdjointFlow`/`feature/HQE`. This app stays on periodic DWF ensembles (NERSC, `bc = −1`;
  decided, §8.6), so none of this is needed.
- **DWF conserved currents.** `FermionFlow_SParameter_DWF.cpp` exits early ("fermion flow for 5D propagators not implemented"):
  conserved currents need the 5D propagator, and only the 4D one is flowed. This app only contracts local operators on 4D
  propagators, so it is unaffected. Do not add flowed conserved-current contractions here.
- **Flowed wall sources and `StochasticMeson`.** `SParameterGF` flows its point sources for the mixed-flow `stoch_*` 2pts. The HQE
  double ratios only need `Meson` on flowed propagators (unflowed wall, flowed sink bilinear) and `HeavyMeson3pt` takes the
  **unflowed** `sinkNoise`, so the wall sources are not flowed here.

Issues seen in `SParameterGF` that this port should avoid (not fixed there; that repo was not changed):
- τ labels use `setprecision(2)`, so a `step_size` below 0.01 (or not a multiple of it) gives colliding `_t<τ>` field and module
  names. Build labels from the step index or with enough digits.
- `source_loc` is still overwritten in the source loop (dead code there too).
- `CL_gamma_pairs`/`snkCL_gamma_pairs` are built but unused.
- Stale comments and docs: `FermionFlow_Condensates.cpp` says "antiperiodic time" with `boundary = "1 1 1 0"`; `Spar.md` still
  describes translating the gauge field and cites a `stochasticMeson.md` that is not in the repo.
- `FermionFlow_SParameter_Z2.cpp` (untracked) does not compile (undefined `stochSrcs`, malformed `for`, redeclared names, missing
  `n_sources` field), and the DWF variant references `MesonPar` fields (`moms`, `ss`) that no longer exist. Neither is in
  `Makefile.am`.

---

## 2. P0: port to the `feature/HQE` Hadrons (no new physics)

### 2.1 Flow modules
- **Module.** `MGradientFlow::WilsonFermionFlow` no longer exists; use `MGradientFlow::FermionFlow`. Its parameters are `props`,
  `propTypes`/`defaultType` (set `defaultType = "PropagatorField"`), `outProps`, `gauge`, `bc`, `steps`, `step_size`,
  `meas_interval` and `output`.
  - Outputs: `<name>`, `<name>_U` (the flowed gauge field at the module's **last** step), and the flowed fields.
  - Without `outProps`, the flowed fields are named `<q>_t<step_size*i>`, counting from the module's **own** start. That is wrong
    for chained modules, so always pass `outProps`.
- **Chain one module per step, as `SParameterGF` does** (§1.1). For step $n$:
    ```
    gauge         = (n == 1) ? "trans_gauge" : flowName(n-1) + "_U"
    props         = (n == 1) ? unflowed b_t, s_t, l_t : outProps of step n-1
    outProps      = q + "_t" + label(n)          (same names as round 1)
    steps         = 1,  meas_interval = 1,  defaultType = "PropagatorField",  bc = flowPar.bc
    ```
  - `<flowName(n)>_U` is then exactly the flowed gauge field at $\tau_n$ that the insertions need.
  - Every flow module goes into its τ's `ResultGroup` (for `gauge_obs`), measured or not, as round 1 and `SParameterGF` do.
  - Measured steps follow the `SParameterGF` schedule: base `meas_interval`, ×4 from `coarsen_step = ceil(Tcoarsen/step_size − 1e-12)`,
    counted from `last_contraction_step`, and always the last step.
  - `label(n)` must be unique per step (see the `setprecision(2)` issue in §1.1).
- **Alternative, selectable at run time:** one module per measurement interval (`steps = n_k − n_{k−1}`, `meas_interval = steps`), which
  cuts the module count. The mode is chosen by `runPar.flowChaining` = `step` | `interval` (§4.1) so both can be benchmarked on the
  same input (§6.5). Everything downstream (names, groups, contractions) is identical in the two modes; only the flow modules differ:
  - `step`: one `FermionFlow` per step, each in its own `FlowTime_t<τ>` group (as `SParameterGF`);
  - `interval`: one `FermionFlow` per measured step $n_k$, covering $n_{k-1}\to n_k$, in the `FlowTime_t<τ_k>` group. Its `gauge_obs`
    holds all the steps of the interval, so unmeasured steps get no group of their own.
- **`gauge_obs` bookkeeping (both options).** `FermionFlow::execute` calls `gauge_status` once on the **input** field (flow time 0) and
  then after every step, with the flow time counted from the module's own start. So every step is measured either way, but:
  - each module's first entry repeats the previous module's last state (and the first module's repeats the τ = 0 `WilsonFlow`);
  - the recorded times are relative ($0,\epsilon,\dots$), not absolute τ.

  The analysis must add the module's start τ and drop each module's first entry. (A small optional Hadrons change, a start-time
  parameter or skipping the initial measurement, would remove this; not needed for P0.)
- **Per-interval constraint.** With `outProps` set, `meas_interval` must equal `steps`: `setup` creates the `outProps` fields once per
  measured step, so more than one measured step per module would create the same names twice and fail.
- Keep `MGradientFlow::WilsonFlow` with `steps = 0` for the τ = 0 gauge observables.
- New helper `make_fermionflow(...)` in `Functions.hpp`.

### 2.2 Bugs and fragilities to fix while porting
1. **Light-quark `WardIdentity` is skipped when `smear == 1`.** Its `createModule`/`push_back` sit inside the `else` branch (misplaced brace).
   Restructure the per-wall block like `SParameterGF`'s: one `WardIdentity` per flavour right after its propagator.
2. **`ti == Tcoarsen` compares doubles.** Replace it with the step-index schedule above (as in `SParameterGF`).
3. **Wall pairing assumes evenly spaced `sources`.** The index `r = (stp*q + dT < Tmax) ? …` is only valid then. Instead look up $t_f=(t_i+\Delta t)\bmod T$ in `sources`, and error if it is missing.
4. **Flowed propagators are found by position.** `qs[y]`, `qs[z]`, `qs[zl]` index by position with `split`. Use a small struct per wall instead: `{t, b, s, l, gsrc}`.
5. **`source_loc` is overwritten** inside the source loop after `CShiftGaugeField` was already created, so it has no effect. Remove that line, or decide what it was meant to do.
6. **Inconsistent 2pt names.** They are `sh_SP_t0_0.00` at τ = 0 but `sh_t0_0.01` at τ > 0. Unify to `sh_<SM>_t<t>_<τ>` (as `SParameterGF`'s
   `<flav>_LL_t<t>_<τ>`), and document the change for the analysis scripts.
7. **`flowPar.Tcoarsen` was an `int`**, so the example's `9.01` could not be read as intended. Now `double` (as in `SParameterGF`).
8. **The example input had no `<ldwfPar>`** although `light = 1`, so the light quark was read with zero parameters. Added a block
   with placeholder values.
- *Kept on purpose:* an empty `gamma_snk` means `all` for the 2pts (unlike `SParameterGF`, which errors).

### 2.3 Build
- Build against a Hadrons install of `feature/HQE`, which in turn needs a Grid with `traceColour`/`traceSpin` (see `newCode.md`).
  `feature/HQE` sits on `test/AdjointFlow`; the open-boundary modules `SParameterGF` uses are **not** in it and are not needed (periodic ensembles, §8.6).
- Update the README's version line.

---

## 3. What has to be computed

### 3.1 Two-quark 3pts (HQE), per (spectator $q$, wall pair $t_i\to t_f$, τ)

$$C[g_{\rm Snk},X,g_{\rm Src}](t)=\sum_{\mathbf x}\mathrm{tr}\big[g_{\rm Snk}\,q_1^{\rm snk}\,(g_{\rm Src}^\dagger\gamma_5)\,S_{h,\tau}(x\leftarrow t_i)^\dagger(\gamma_5X)\,S_{h,\tau}(x\leftarrow t_f)\big],\quad
q_1^{\rm snk}=\sum_{\mathbf z}\eta_f^\dagger(\mathbf z,t_f)\,S_q(\mathbf z,t_f\leftarrow t_i)$$

`HeavyMeson3pt` parameters:

| Parameter | Value | Notes |
|---|---|---|
| `qSrc` | `b_<t_i>_t<τ>` (`b_<t_i>` at τ = 0) | |
| `qSnk` | `b_<t_f>_t<τ>` | |
| `spectator` | `s_<t_i>` or `l_<t_i>` | **unflowed** at every τ: only the operator is flowed |
| `sinkNoise` | `sm_z2_p000_t<t_f>_<t_f>` if `smear == 1`, else `z2_p000_t<t_f>` | the conjugate is taken inside the module |
| `tSnk` | $t_f$ | |
| `gamma` | see below | |

The `gauge` is `<flowName>_U` at τ > 0 and `trans_gauge` at τ = 0. Use `bc = flowPar.bc` and `flowTime = τ`.

**Gamma triads** (`"sink vertex source"`). Source and sink are the **probe gammas** only: at most $P$ = `Gamma5` and $A_0$ = `GammaTGamma5`,
possibly only $P$ (to be decided, input `probeGammas`). The vertex sets are `all` (the 16 Γ, default) and `ops` (the HQE operators). So:
- probes = {P}: `Gamma5 all Gamma5`, `Gamma5 ops Gamma5`, i.e. 16 + 13 = 29 correlators per instance;
- probes = {P, A₀}: those for all four sink/source combinations, 4 × 29 = 116 correlators.

They are only extra `traceSpin` calls, which are cheap.

| Vertex | Notes quantity |
|---|---|
| `all` (16 Γ) | $\langle\bar QQ\rangle$ (`Identity`), $\gamma^0$ normalisation (`GammaT`), $\mathbb D_2[1-\gamma^0]$, mixing subtractions |
| `D2` | $\mu_\pi^2$; $\mathcal R_D$ |
| `D4`, `D4L` | mass check $C[iD_0]/C[\gamma^0]\to m_Q$; Re parts |
| `D4D4`, `D4D4D4`, `D2D4` | full-$D$ appendix, $\mathcal R_D$ |
| `sigmaB`, `sigmaE` | $\mu_G^2$, $\bar Q\mathcal EQ$ |
| `sigmaBD4`, `sigmaED4` | $\mathcal R_{LS}$, full-$D$ $\rho_{LS}$ |
| `DiE`, `divE` | $\rho_D^3$ invariant part (two discretisations: the difference is a cross-check) |
| `DiSigmaE` | $\rho_{LS}^3$ invariant part |

**Spectators.** `s` always, and `l` if `light == 1`. One `HeavyMeson3pt` instance takes one spectator, so this means one instance per spectator.
The heavy derivatives are then built twice per wall pair; this is the P3 optimisation.

**Probe 2pts.** These already exist: the flowed `sh`/`lh` `Meson` outputs from walls $t_i$ and $t_f$, and `hh`/`ss`/`ll` for PCAC.
The `mesonPar` `gamma_src`/`gamma_snk` must contain every (sink, source) combination of the probe gammas: (G5 G5), plus (GTG5 G5), (G5 GTG5)
and (GTG5 GTG5) if $A_0$ is used. The example input does.

### 3.2 Four-quark 3pts, per (spectator $q$, wall pair, τ)
- **Dim 6 ($\Delta B=0,2$): `WeakNonEye3pt`, as in round 1** (round 2 uses new lattices, so these are production data, not only a
  cross-check). One call per operator and per (`gammaIn`, `gammaOut`) in the probe gammas.
- **Dim 7 ($\Delta B=2$): `WeakNonEye3ptDeriv`**, the generalised `WeakNonEye3pt` (same legs and gammas, plus gamma pairs and
  derivative legs), one module per probe combination on the `DeltaHs2` legs. It is a separate module from `HeavyMeson3pt`
  because the two share no work ($F$ is only used by the two-quark part, and each differentiates its own legs).

Leg conventions are as in `WeakNonEye3pt`: $\bar b$ on the adjointed legs (`qBarLeft`, `qRight`), and $q$ on `qLeft`, `qBarRight`.

| `name` | `qBarLeft` | `qLeft` | `qBarRight` | `qRight` | round-1 equivalent |
|---|---|---|---|---|---|
| `DeltaHs2` | $h_i$ | $q_i$ | $q_f$ | $h_f$ | `sh_DeltaHs2_G5_G5` |
| `DeltaHs0_T1` | $q_i$ | $h_i$ | $q_f$ | $h_f$ | `sh_DeltaHs0_T1_G5_G5` |
| `DeltaHs0_T2` | $q_i$ | $q_i$ | $h_f$ | $h_f$ | `sh_DeltaHs0_T2_G5_G5` |

Each `WeakNonEye3ptDeriv` module uses `gammaPairs = diag` (the 16 Γ⊗Γ, default) and one (`gammaIn`, `gammaOut`) from the probes. With the
derivative spec `none` it equals the matching `WeakNonEye3pt` output (both traces, 16 G); keep it as a built-in check on every configuration.
The `DeltaHs0` rows are listed for the leg conventions only: at dim 6 they come from `WeakNonEye3pt`, and dim 7 $\Delta B=0$ is far future.
Here $h=$ flowed `b`, $q=$ flowed `s`/`l`, and all legs are at the same τ.

**Dim-7 $\Delta B=2$ derivative specifications: Lenz–Nierste basis** (on `DeltaHs2` only). The labels are those written to the output.

| spec (`legs @ dirs [:couple]`) | purpose |
|---|---|
| `none` | dim-6 $Q_i$; $R_0$; $R_1,\tilde R_1$ (× $m_s/m_b$ from PCAC) |
| `qLeft@3`, `qBarRight@3` | leading-power $R_{2,3}$: $v\!\cdot\!D$ on the light leg, both Wick assignments |
| `qBarLeft@3`, `qRight@3` | $\nabla_4$ on each $\bar b$ leg alone: the $\mp M$ check (`lifetimes_round2.md` §4.3, §7.6) |
| `qBarLeft qLeft@0123` (A), `qRight qLeft@0123` (B), `qBarRight qRight@0123` (C), `qBarLeft qBarRight@0123` (D) | full $R_2,R_3$ and $\tilde R_2,\tilde R_3$ |
| *(future, not in the default list)* `qLeft@0123:G2`, `qBarRight@0123:G1` | index-coupled operators, for a basis beyond Lenz–Nierste (§8.2) |

- **Why A–D are complete.** They are the four ways of assigning the derivative-carrying bilinear $(\bar b\,\overleftarrow D,\,D s)$ to the legs.
  The colour-rearranged $\tilde R$ and the Fierz partners only need the *other* trace of the same specification, and the module writes both.
- **Why `diag` is enough.** The $R_{2,3}$ Dirac structures ($\gamma^\mu P_L\otimes\gamma_\mu P_L$ and $P_L\otimes P_L$) reduce to $VV$+$AA$ and $SS$+$PP$ in the
  parity-even pseudoscalar-to-pseudoscalar matrix element, and $D_\rho\cdots D^\rho$ is parity-even. Index-coupled specs create $\Gamma_1\neq\Gamma_2$ through
  `couple` anyway. **Check this parity argument once, with `gammaPairs = all` on one configuration.**
- **Cost per entry.** It is (1 + 4 + 4×4 directions) × 16 pairs × 2 traces ≈ 670 traces. Each distinct leg is differentiated once per ρ for all specs together.
- **$\Delta B=0$ dim-7** (spectator effects at $1/m_b$): far future, not for now. It would be the same spec list on `DeltaHs0_T1`/`T2`.

### 3.3 Mapping to the analysis
Nothing new is needed on the Hadrons side for the double ratios.
- Two-quark: $\mathbb D_k$ from the `hqe` correlators and the flowed probe 2pts. $r_Q$ comes from the hh/hl PCAC ratios.
- Four-quark: $\frac{1}{m_b^2}\frac{\langle R_2\rangle}{\langle Q_1\rangle}=\frac{1}{r_Q^2M^2}\frac{\mathcal S_{AA}[\mathcal O_{R_2}]}{\mathcal S_{PP}[Q_1]}$. The probe swap uses the same
  four-quark correlators with the (GTG5) vs (G5) probe 2pts. $m_s/m_b$ comes from the ss/hh flowed PCAC ratios.
- The Wick assignments of A–D and the leading-power sign $\pm$ (`lifetimes_round2.md` §4.3–4.4) are combined in the analysis.

---

## 4. App code changes (P1, P2)

### 4.1 `Inputs.hpp`: a new `HQEPar` block

```cpp
class HQEPar : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(HQEPar,
        std::string, probeGammas,   // source/sink of all 3pts: "Gamma5" or "Gamma5 GammaTGamma5" (required)
        bool,        twoQuark,      // HeavyMeson3pt two-quark instances (HQE)
        std::string, vertices,      // two-quark vertex sets; empty = "all ops" (16 Γ + HQE operators)
        bool,        dim7,          // WeakNonEye3ptDeriv DeltaHs2 dim-7 modules (Lenz–Nierste specs, §3.2)
        std::string, gammaPairs,    // four-quark Γ1⊗Γ2 for dim 7; empty = "diag" (16 Γ⊗Γ)
        bool,        atTau0);       // also measure the new 3pts at tau = 0
};
```
Also add to `RunPar`:
```cpp
        std::string, flowChaining   // "step" | "interval" (required; anything else stops the run)
```

Read the HQE block as `<hqePar>` in `main`. The dim-6 `WeakNonEye3pt` calls always run and use `probeGammas` for `gammaIn`/`gammaOut`.
App-level defaults follow round 1 (an empty `gamma_snk` means `all`): an empty `vertices` means `all ops` (16 Γ + the HQE
operators), an empty `gammaPairs` means `diag` (16 Γ⊗Γ).
`probeGammas` has no default; an empty value, or a gamma other than `Gamma5`/`GammaTGamma5`, stops the run.
The Hadrons modules themselves always receive explicit values.
(Dropped: `dim7DeltaB0` (far future), `dim7Coupled` (future basis option, §8.2), `legacyWeakNonEye` (`WeakNonEye3pt` is now always on).)

### 4.2 `Functions.hpp`: new helpers
- `make_fermionflow(app, name, gauge, props, outProps, nSteps, flowPar)` (§2.1): `steps = meas_interval = nSteps` (1 in `step` mode),
  `defaultType = "PropagatorField"`.
- `measuredSteps(flowPar)`: the `SParameterGF` schedule (§1.1) as a list of step indices, and `tauLabel(n)` (unique per step).
- `make_hm3pt(app, name, gauge, bc, tau, qSrc, qSnk, spectator, sinkNoise, tSnk, triads)` wraps
  `createModule<MContraction::HeavyMeson3pt>` and returns the name.
- `deriv(legs, dirs, couple)` returns a `MContraction::WeakDerivSpec`; `none` is `{legs = "none", dirs = "", couple = ""}`, and
  every other spec sets `couple` explicitly (`none`, `G1` or `G2`).
- `dim7Specs()` returns the Lenz–Nierste list of §3.2 (the coupled specs can be added later behind a switch).
- `make_weaknoneye_deriv(app, qbl, ql, qbr, qr, gIn, gOut, gammaPairs, specs, gauge, bc, name)` wraps
  `createModule<MContraction::WeakNonEye3ptDeriv>` (argument order as `make_weaknoneye`).
- `probeTag(g)`: `Gamma5` → `G5`, `GammaTGamma5` → `GTG5` (round-1 tokens).
- Both modules have no defaults, so every parameter they use must be set. `output` and `momProjector` are left empty:
  results go through `ResultGroup`, at zero momentum only.

### 4.3 `FermionFlow_Lifetimes.cpp`: main loop

**Restructure as:**
```
walls    = [{t, b, s, l, z2, gsrc}]          // from configPar.sources
pairs    = [(wi, wf, dt)]                    // t_f = (t_i + dt) mod T looked up in walls
measured = measuredSteps(flowPar)            // SParameterGF schedule (§1.1)

τ = 0: gauge_τ = "trans_gauge"; group = [WilsonFlow(steps=0), WardIdentities, 2pts, 3pts below]
prev = 0
for n in 1..flowPar.steps:
    if runPar.flowChaining == "step":
        flow = make_fermionflow(n-1 → n, nSteps = 1)
        if n not in measured: corrgroups += ResultGroup("FlowTime_t" + tauLabel(n), [flow]); continue
    else:  // "interval"
        if n not in measured: continue
        flow = make_fermionflow(prev → n, nSteps = n - prev)
    prev = n
    group = [flow]
    gauge_τ = flow + "_U"
    2pts as now (names unified, §2.2.6)
    for (wi, wf, dt) in pairs:
        for q in {s, l if light}:
            for (gIn, gOut) in probes × probes:
                group += WeakNonEye3pt DeltaHs2, DeltaHs0_T1, DeltaHs0_T2   // names below
            (twoQuark and (τ > 0 or atTau0))
                group += make_hm3pt(q+"h_HQE_t"+t_i+"_dt"+dt+"_"+τ,
                    qSrc = h(wi,τ), qSnk = h(wf,τ), spectator = q(wi) unflowed,
                    sinkNoise = noise(wf), tSnk = t_f, gamma = triads(probes, vertices))
            (dim7 and (τ > 0 or atTau0))
                for (gIn, gOut) in probes × probes:
                    group += make_weaknoneye_deriv(h(wi,τ), q(wi,τ), q(wf,τ), h(wf,τ), gIn, gOut, gammaPairs, dim7Specs(),
                                                   q+"h_DeltaHs2_dim7_"+gIn+"_"+gOut+"_t"+t_i+"_dt"+dt+"_"+τ)
    corrgroups += ResultGroup("FlowTime_t" + τ, group)
WriteResultGroup(corrgroups)
```

- `h(w,τ)` is `b_<t>` at τ = 0 and `b_<t>_t<τ>` otherwise; likewise `q(w,τ)`.
- The spectator is always the unflowed `s_<t_i>`/`l_<t_i>`.

**Module and correlator names** (round-1 style: `<spectator><heavy>_<operator>_<probes>_t<t_i>_dt<Δt>_<τ>`, with $\tau$ from `tauLabel`):

| What | Module / result name | Inside the result |
|---|---|---|
| 2pts | `<flav>_<SM>_t<t>_<τ>` (§2.2.6) | as round 1 |
| dim-6 four-quark | `sh_DeltaHs2_G5_G5_t<t>_dt<Δt>_<τ>`, `sh_DeltaHs0_T1_…`, `sh_DeltaHs0_T2_…`; with $A_0$ also `…_G5_GTG5_…`, `…_GTG5_G5_…`, `…_GTG5_GTG5_…` | as round 1 (16 G, both traces) |
| two-quark HQE (`HeavyMeson3pt`) | `sh_HQE_t<t>_dt<Δt>_<τ>` | `gamma[k] = {sink, vertex, source}` and `corr[k][0][t]`; vertex = gamma name or operator name (`D2`, `sigmaB`, …) |
| dim-7 $\Delta B=2$ (`WeakNonEye3ptDeriv`) | `sh_DeltaHs2_dim7_G5_G5_t<t>_dt<Δt>_<τ>`; with $A_0$ also `…_G5_GTG5_…` etc. | `in`, `out`, and `contractions[]` entries `{deriv = none \| qLeft@3 \| … \| qBarLeft qLeft@0123, gamma1, gamma2, trace, corr[0][t]}` |

`lh_…` likewise for the light spectator. With probes = {P}, the dim-6 names are exactly round 1's.
The analysis reads `HeavyMeson3pt` entries by their `gamma` triad and `WeakNonEye3ptDeriv` entries by `deriv`/`gamma1`/`gamma2`/`trace`,
not by position.

### 4.4 Example input
- Add `<hqePar>` to `run/input_fermion.xml`: everything on, empty `vertices`/`gammaPairs`. The example uses `probeGammas = Gamma5`
  until §8.4 is settled; set `Gamma5 GammaTGamma5` to add $A_0$.
- Keep `<mesonPar>` gamma pairs as they are; they already include the PCAC/probe pairs.

---

## 5. Cost and memory (32³×64, double precision, one propagator ≈ 4.8 GB in total)

- **Inversions: none extra.** Everything reuses the round-1 propagators.
- **Per module:**
  - `HeavyMeson3pt`: about 5 propagators, one per source gamma (≤ 2), plus 6 colour-matrix fields (= 6/16 of a propagator in total);
  - `WeakNonEye3ptDeriv`: `L`, `R`, plus one propagator per distinct differentiated leg (4), i.e. ≈ 6;
  - so a peak of ≈ 7 propagators ≈ 35 GB of temporaries, freed after each module.

  With the naive scheduler, only one instance's temporaries are alive at a time.
- **Persistent fields:**
  - the flowed props of all walls at the current τ (as now);
  - **new:** the unflowed spectator `s_<t>`/`l_<t>` of every wall now stays alive for the whole run, because every τ uses it.
    That is $n_{\rm walls}\times n_{\rm spec}$ extra propagators. If this is too much, use the P3 sliced-spectator option.
  - **new:** the `sinkNoise` field of every wall (smeared or plain Z2 wall, propagator-sized) also stays alive for the whole
    run, for the same reason: $n_{\rm walls}$ extra propagators. The P3 sliced-spectator option removes these too.
- **Derivative work per instance:**
  - two-quark: about 3 spatial + 1 temporal + nested heavy derivatives;
  - four-quark: 4 legs × 4 directions per `WeakNonEye3ptDeriv` module, i.e. repeated for each probe pair (small next to the traces).

  This is all nearest-neighbour shifts, and small compared with one Möbius solve.
- **Duplication:** with `light == 1`, the heavy derivatives are computed twice per (pair, τ). P3 removes this.
- **Module count per τ:** walls × Δt × spectators × (3 × n_probes² `WeakNonEye3pt` + n_probes² `WeakNonEye3ptDeriv` + 1 `HeavyMeson3pt`);
  for the example input with two probes, 2 × 1 × 2 × (12 + 4 + 1) = 68, compared with 12 in round 1. The `WeakNonEye3pt` calls dominate the count but are cheap.

---

## 6. Validation of the app (in order)

1. **P0 regression.** Run on one configuration with the old binary (round-1 Hadrons) and the ported one, both without HQE.
   - The 2pts and `WeakNonEye3pt` outputs at every τ must agree to rounding.
   - This checks `FermionFlow` chaining, `outProps` naming and the measurement schedule.
2. **P1/P2.** Every dim-7 `WeakNonEye3ptDeriv` entry with `deriv = none` must equal the matching `DeltaHs2` `WeakNonEye3pt` output (both traces, 16 G).
   Both always run, so this check is available on every configuration.
3. **P2 sanity, on one configuration** (the exact identities are already in `Test_heavy_meson_3pt`):
   - `GammaT` vertex 3pt / 2pts → plateau at the vector-current normalisation;
   - $C[D4]/C[\gamma^0]$ consistent with the PCAC mass;
   - $\mu_\pi^2>0$ and $\langle\bar QQ\rangle<1$ (`lifetimes_round2.md` §7.4);
   - `DiE` vs `divE` agree up to discretisation;
   - `qBarLeft@3` and `qRight@3` show the $\mp M$ behaviour.
4. **Parity check.** Run once with `gammaPairs = all` for the dim-7 specs. The off-diagonal pairs must vanish within noise; if they don't, keep `all` (§3.2).
5. **Timing.** Record the `HeavyMeson3pt` (Sink, Trace, Insertion, FT, I/O) and `WeakNonEye3ptDeriv` (Insertion, Trace, FT, I/O) timers against the solves. Decide whether P3 is worth doing.
6. **Flow chaining benchmark.** Run the same input with `flowChaining = step` and `interval`. All correlators must agree to rounding,
   and the `gauge_obs` must agree once each module's first entry is dropped and its times are shifted by its start τ (§2.1).
   Compare total wall time and peak memory; keep the faster mode as the production setting.

---

## 7. P3 (optional, later)

| Item | Where | Why |
|---|---|---|
| `spectator`/`sinkNoise` as lists (one $q_1^{\rm snk}$ each; heavy derivatives shared) | Hadrons `HeavyMeson3pt` | removes the s/l duplication of heavy derivatives |
| Accept a precomputed `SlicedPropagator` spectator ($q_1^{\rm snk}$ via `MSink::Smear` + adjoint noise sink) | Hadrons (`SourceSink` adjoint option + `HeavyMeson3pt`) | frees the unflowed spectator fields after τ = 0 |
| ✅ Stochastic part: dim-6 $\Delta B=0$ eyes and condensates (§7.1, §7.2) | app + extra inversions + `SourceSink` `sourceOp` | SU(3)-breaking eyes |
| **Far future:** $\Delta B=0$ dim-7 (connected and eye), Darwin four-quark cross-check | app + new Hadrons eye code | lifetimes at $1/m_b$ |
| $B^*$ / other probes | app (`HeavyMeson3pt` `gamma` triads with `GammaX/Y/Z`) | vector-meson HQE parameters |

### 7.1 Stochastic part: eyes

**Only the eye is needed for the lifetimes.** The condensates are measured anyway, as diagnostics (decided 2026-11-02, see below).
- The flowed $\Delta B=0$ operator's mixing into $\bar b\Gamma b$ (dim 3, power-divergent in $1/\tau$) is handled by the small-flow-time
  matching: perturbative coefficients times the standard $\langle B|\bar b\Gamma b|B\rangle$ 3pt, which `HeavyMeson3pt` already computes
  (`all` vertex). The vacuum part of the light-quark loop is part of the four-quark matrix element itself, not something to subtract
  with a separately measured condensate. (For $\gamma^\mu P_L\otimes\gamma_\mu P_L$ the scalar part of $\langle L\rangle$ drops out of the contraction anyway.)
- The kinetic trace $\langle\eta^\dagger\not\!D[V_\tau]\phi\rangle$ for $Z_\chi(\tau)$ cancels in the double ratios.

**Diagnostics, almost free (included).** The flowed condensates $\langle\bar q\Gamma q\rangle_\tau$ and the $Z_\chi$ kinetic trace can be
measured on the eye's flowed pairs at no extra inversion cost. They are not inputs to the lifetime analysis.

**Cost of the optional diagnostics.** They reuse the eye's noise, solves and flowed pairs; per (noise component, flavour, measured τ) it is
one `StochasticCondensatePropagator` (`gammas = "all"`, `c_fl` from input) and one `DslashFieldPropagator` + `StochasticCondensatePropagator`,
exactly as in `SParameterGF/FermionFlow_Condensates.cpp`. Both modules are on `test/AdjointFlow`. No new Hadrons code.

**Four-quark vacuum condensates: not proposed here.**
- They need $\sum_x$ of products of two loops at the same site, $\mathrm{tr}[\Gamma L_A(x)]\,\mathrm{tr}[\Gamma L_B(x)]$ and $\mathrm{tr}[\Gamma L_A(x)\Gamma L_B(x)]$,
  with **independent** noises $A\neq B$ (otherwise biased). That doubles the inversions, and no existing module does the contraction
  (`StochasticCondensate` sums over $x$ first). It would need a small new Hadrons module taking two `Loop` fields.
- The signal is a vacuum disconnected quantity and noisy.
- For lifetimes they only enter through mixing of light-quark four-quark operators into the identity, which these $b$-containing operators
  do not have. If wanted for other physics, put them in a generic condensates program (`SParameterGF`-style), not in this app.

**Where it runs.** The eye needs the main run's heavy and spectator wall propagators. Either the stochastic part runs in the main
app, or the walls are saved to disk and a second executable loads them. Both are implemented, to be benchmarked (§7.2). In both,
the noise and solution fields are flowed in the **same** `FermionFlow` as the heavy propagators, so all fields share the gauge RK stages.

### 7.2 Eye implementation (2026-11-02): written, not compiled

**Contraction.** `WeakEye3pt` is the `Gamma3pt` contraction with the vertex $g_V$ replaced by $G_1\,L\,G_2$ (trace 1) or $G_1\,\mathrm{tr}(LG_2)$ (trace 2), for the operator gamma pairs `stochasticPar.eyeGammaPairs` (`diag`, `all` or a list; required, 2026-11-03; `all` is needed for the penguin-type operators via Fierz, `hqe_hadrons_map.tex` item 19):
- `qBarLeft` = heavy from $t_i$ and `qBarRight` = heavy from $t_f$, both at τ;
- `qSpectator` = the unflowed spectator from $t_i$, closed at $t_f$ against the adjoint of $t_f$'s sink noise:
  `MSink::SourceSink` (`source` = `gsrc` of $t_f$, `sourceOp = adj`) + `MSink::Smear`. This is the closure `HeavyMeson3pt` uses
  internally. The `adj` is essential: `MSource::Z2` noise is complex, and without it the eye averages to zero
  (RareK avoided this with Coulomb-gauge-fixed deterministic walls);
- `loop` $L=\frac1{n_{\rm src}}\sum_k\phi_k\eta_k^\dagger$ at τ (`MContraction::Loop` per component, `MUtilities::LinearComb` for the sum),
  with $\phi_k=D^{-1}\eta_k$ and both $\phi_k$ and $\eta_k$ flowed, so $L\to K S K^\dagger$.
- One module per (wall pair, spectator $f$, loop flavour $g\in\{s, l\}$, probe pair), at τ = 0 and every measured τ.
  No Wick signs are applied (as for `WeakNonEye3pt`).

**Condensates** (exactly `SParameterGF/FermionFlow_Condensates.cpp`), per noise component $k$, loop flavour and τ:
`StochasticCondensatePropagator` with `gammas = all` (`c_fl` from input at τ > 0, 0 at τ = 0) and the `DslashField` trace.

**Inputs** (`Inputs.hpp`, example in `run/input_fermion.xml`):
- `runPar.eyes` (0/1): eyes and condensates in this run; needs `<stochasticPar>` (`nSrcs`, `nsparse`, `c_fl`, as in `SParameterGF`, and `eyeGammaPairs`).
- `runPar.saveWalls` (0/1): save the unflowed wall propagators (`b`, `s`, and `l` if `light`) to `<wallStem>/<prop>.<traj>.bin`
  (`MIO::SavePropagator`, or `SavePropagatorIo32` with `runPar.wallIo32 = 1`). The directory is created if needed.
- `FermionFlow_Eyes <same input>`: loads those walls (`MIO::LoadPropagator[Io32]`), regenerates the Z2 wall sources, then does the
  noise, solves, flow and eye/condensate contractions only. It needs the **same `runId`** as the saving run, because the Z2 wall noise
  (used for the spectator sink) is seeded from `runId`, module name and trajectory. Output `<mesonPar.file_name>_eyes`; the
  application database, if set, becomes `<stem>_eyes.db`.

**Names** (suffix `_t<t_i>_dt<Δt>_<τ>` as for the other 3pts):
- eyes: `<f>h_DeltaHs0_eye_<g>_<G>_<G>_…`, e.g. `sh_DeltaHs0_eye_l_G5_G5_t0_dt32_0.00` (loop flavour $l$, spectator $s$; the `<G>` are the probe gammas). Entries carry metadata `{in, out, gamma1, gamma2, trace}`;
- condensates: `condensate_<g>_stoch_<k>_t<τ>` and `dslash_condensate_<g>_stoch_<k>_t<τ>`, as in `SParameterGF`
  (estimators: sum over $k$, divide by `nSrcs`).

**Code** (`Functions.hpp`): `Wall`/`WallPair` (moved from the main file), `make_wall`, `make_wall_pairs`, `parse_probes`,
`make_flow_chain` (the flow loop, shared by both executables), `make_save_wall`/`make_load_wall`, and the eye section
(`make_stochastic_noise_diluted`, `make_eye_noise`, `eye_flow_fields`, `make_eye_spectators`, `make_eyes`, and the
`SParameterGF` helpers `make_dslash_field`, `make_stochastic_condensate`).

**Cost** (32³×64, double, one propagator ≈ 4.8 GB):
- inversions: $n_{\rm src}\,n_{\rm sparse}^4$ components × 12 columns × $n_{\rm loop}$ flavours; with the example ($n_{\rm src}=1$,
  $n_{\rm sparse}=2$, $s$ and $l$) 16 × 12 × 2 = 384 solves, against 12 × 3 per wall for the walls;
- memory: the flowed noise and solutions are $(1+n_{\rm loop})\,n_{\rm src}\,n_{\rm sparse}^4$ propagators, 48 ≈ 230 GB for the example,
  held for the whole flow (and twice during each `FermionFlow` step). This, not the contractions, is likely to set the node count.
  `nsparse = 1` gives 3 propagators but noisier loops;
- disk for the saved walls: $n_{\rm walls}\times(1+n_{\rm spec})$ propagators per configuration, 6 × 4.8 GB (2.4 GB with `wallIo32`) for the example;
- `FermionFlow_Eyes` repeats the gauge flow and flows only the heavy walls (not the spectators), but no wall solves.

**Checks (needs a build).**
1. Same input, same `runId`: `runPar.eyes = 1` in the main run vs `FermionFlow_Eyes` from walls saved with `wallIo32 = 0`. The noise is
   seeded identically, so every eye and condensate must agree to rounding (to single precision with `wallIo32 = 1`).
2. The condensates at τ = 0 vs `SParameterGF` on the same configuration and quark parameters (statistically).
3. Benchmark: total time and memory of (main with eyes) vs (main with `saveWalls`) + `FermionFlow_Eyes`.

---

## 8. Decisions needed before coding

1. **Operator gauge field at τ = 0.** *For now: `trans_gauge`.* **Open, to discuss later.**
   - `trans_gauge`: continuous with τ > 0 (the flow starts from it), matches the small-flow-time matching and the light-quark actions.
   - `Vgauge` (stout, used only by the heavy Möbius action): the heavy propagator obeys its equation of motion with these links, so
     EOM-based checks ($C[iD_0]/C[\gamma^0]\to m_Q$, the `D4`/`D4L` remainders) have smaller cut-off effects.
   - Middle option: `trans_gauge` everywhere plus, at τ = 0 only, one extra two-quark instance per wall pair on `Vgauge` with the `ops`
     vertices. No extra solves.
   - τ = 0 is a reference point, not used in the flow analysis, so this only affects how the τ = 0 checks are read.
2. **Dim-7 basis.** *Decided for now: Lenz–Nierste only* ($R_0..R_3$, $\tilde R_{1,2,3}$; specs A–D, leading power, $\mp M$ check, `none`).
   **Future option:** index-coupled operators, which contract the derivative index with a gamma index (e.g. $(\bar b\gamma^\mu P_L iD^\nu s)(\bar b\gamma_\nu P_L s)$).
   They appear in other bases and EOM reductions. In `WeakNonEye3ptDeriv` they are `couple = G1/G2` specs (`qLeft@0123:G2`, `qBarRight@0123:G1`),
   about 128 extra traces per spec; adding them later is a change to `dim7Specs()` plus the analysis and matching.
3. ~~**$\Delta B=0$ dim-7** in the first production, or later?~~ Decided: far future. Beyond dim-6 only $\Delta B=2$ for now.
4. ~~**Gamma triads.**~~ *Decided:* source/sink only the probes, at most $P$ and $A_0$ (possibly only $P$, **still to decide**, input
   `probeGammas`); the 16 Γ (`all`) plus `ops` at the two-quark vertex and `diag` (16 Γ⊗Γ) for the four-quark insertion by default.
5. ~~**Output names.**~~ *Decided:* round 2 uses new lattices, so `WeakNonEye3pt` stays in production for dim 6. Names follow round 1,
   extended for the probes and the new modules (`sh_HQE_…`, `sh_DeltaHs2_dim7_<G>_<G>_…`), see §4.3.
6. ~~**Ensembles.** Periodic DWF or also CLS open-boundary?~~ Decided: this app stays on periodic DWF ensembles. No rebase onto
   `feature/contractConservedCurrents` and no `bc = 0` support in the HQE modules are needed.
7. **Flow chaining.** One `FermionFlow` per step as in `SParameterGF` (proposed for P0), or per measurement interval later (§2.1)?
   Both measure the gauge observables at every step (`gauge_status` runs after each step inside a module), with the same caveats:
   per-module relative flow times and a repeated first entry (§2.1). The flow work is identical, so the choice is only module overhead
   vs. a slightly more complex schedule (`meas_interval = steps` per module). Both are implemented behind `runPar.flowChaining`
   and benchmarked (§6.6); `step` is the reference for the P0 regression.

---

## 9. Implementation log

### P0 + P1 (2026-10-28): written, not compiled

What the code does now (`FermionFlow_Lifetimes.cpp`):
- **Input checks before `Grid_init`**, each stopping the run with a message:
  - `runPar.flowChaining` ∈ {`step`, `interval`};
  - `flowPar.steps` ≥ 0, `step_size` > 0, `meas_interval` ≥ 1 (when flowing), `bc` = ±1;
  - `sources` non-empty, all in $[0,T)$;
  - `hqePar.probeGammas` non-empty, only `Gamma5`/`GammaTGamma5`, no repeats;
  - `hqePar.dim7 = true` (P2 not implemented yet; check removed in P2);
  - unless `gamma_snk` is empty (= `all`), the 2pt gamma pairs must contain every probe pair;
  - every $t_f=(t_i+\Delta t)\bmod T$ must be a source wall.
- **Walls.** `struct Wall {t, tag, z2, gsrc, q["b"/"s"/"l"]}`; `struct WallPair {i, f, dT}` from the lookup above (bugs 3, 4).
  `WardIdentity` for every flavour right after its propagator (bug 1). The `source_loc` overwrite is gone (bug 5).
- **Contractions** in one lambda `contract(label, gauge, τ, τ0, group)`, used for τ = 0 and every measured τ:
  - 2pts `sh/ss/hh/lh/ls/ll` named `<flav>_<SM>_t<t>_<τ>` at every τ (bug 6; `hh` is always `PP`); `SS` 2pts at τ = 0 only;
  - dim-6 `WeakNonEye3pt` for every (`gammaIn`, `gammaOut`) in the probes, named `<q>h_<Op>_<G5|GTG5>_<G5|GTG5>_t<t>_dt<Δt>_<τ>`;
  - `HeavyMeson3pt` two-quark instance `<q>h_HQE_t<t>_dt<Δt>_<τ>` if `twoQuark` (and `atTau0` at τ = 0): heavy legs at τ, unflowed
    spectator `q_<t_i>`, `sinkNoise` = the wall source of $t_f$ (smeared if `smear == 1`), `tSnk` = $t_f$, triads
    `probe × vertices × probe`, gauge `trans_gauge` at τ = 0 and `<flow>_U` at τ > 0, `bc = flowPar.bc`, `flowTime = τ`.
- **Flow.** Base props `b, s, (l)` of every wall. `step`: one `FermionFlow` per step; `interval`: one per measured step covering the
  steps since the last. Each module takes the previous outputs and `_U`; `outProps = <prop>_t<label>`. Every created flow module is
  in its `FlowTime_t<label>` group. Labels from `tauLabel` (≥ 2 decimals, more if `step_size` needs them).
- **Measured steps** from `measuredSteps` (`SParameterGF` schedule). For the P0 regression note that round 1 used
  `t % meas_interval == 0` with `meas_interval *= 4` at `ti == Tcoarsen`; the measured steps agree for the example input, but in
  general compare only the τ measured by both.

Deviations from the plan text and things to check on the first build:
- The example input uses `probeGammas = Gamma5` (decision §8.4 still open) and `flowChaining = step`.
- The `<ldwfPar>` values in the example input are placeholders.
- First-build suspects: `strToVec<Gamma::Algebra>` for the probe names; `MGradientFlow::FermionFlow::Par`/`WilsonFlow::Par` field
  names against the installed Hadrons; `std::map` in `Wall`.

### P2 (2026-10-28): written, not compiled

- `vertices` now defaults to `all ops` (16 Γ + the HQE operators); `gammaPairs` defaults to `diag`. `dim7` is no longer rejected.
- New helpers in `Functions.hpp`: `deriv(legs, dirs, couple)`, `dim7Specs()`, `hm3ptFourQuark(...)`, `make_hm3pt_fourquark(...)`.
- `dim7Specs()` is the Lenz–Nierste list of §3.2, in this order (labels as written by the module):
  `none`, `qLeft@3`, `qBarRight@3`, `qBarLeft@3`, `qRight@3`, `qBarLeft qLeft@0123` (A), `qRight qLeft@0123` (B),
  `qBarRight qRight@0123` (C), `qBarLeft qBarRight@0123` (D). All with `couple = none`; the index-coupled specs stay a future option (§8.2).
- If `dim7` (and `atTau0` at τ = 0): one `HeavyMeson3pt` instance `<q>h_DeltaHs2_dim7_t<t>_dt<Δt>_<τ>` per (spectator, wall pair, τ),
  with legs `qBarLeft = h_i`, `qLeft = q_i`, `qBarRight = q_f`, `qRight = h_f` (flowed, as `WeakNonEye3pt` `DeltaHs2`), one entry
  `DeltaHs2_<G5|GTG5>_<G5|GTG5>` per probe pair, gauge `trans_gauge` at τ = 0 and `<flow>_U` at τ > 0. Its `hqe` block is left empty
  (skipped by the module).
- Example input: `dim7 = true`, empty `vertices`/`gammaPairs`.
- To check on the first run: every `deriv = none` entry equals the matching `<q>h_DeltaHs2_<gIn>_<gOut>` `WeakNonEye3pt` (§6.2).

### Eyes and condensates (2026-11-02): written, not compiled

See §7.2. Hadrons: `MSink::SourceSink` got a required `sourceOp` (`adj` | `none`). App: `runPar.eyes`, `saveWalls`, `wallStem`,
`wallIo32`, `<stochasticPar>`; new executable `FermionFlow_Eyes`; the flow loop moved into `make_flow_chain`.

### Module split (2026-10-29): written, not compiled

- `HeavyMeson3pt` lost its four-quark part, and its former `hqe` block is now top-level parameters (`qSrc`, `qSnk`, `spectator`,
  `sinkNoise`, `tSnk`, `gamma`). Its result is `{flowTime, tSnk, momentum, gamma, corr}`.
- The four-quark part is the new module `MContraction::WeakNonEye3ptDeriv`: `WeakNonEye3pt` plus `gammaPairs`, `derivatives`,
  `gauge`, `bc` and `momProjector`. Like `WeakNonEye3pt`, one module computes one (`gammaIn`, `gammaOut`) contraction.
- App changes (`Functions.hpp`, `FermionFlow_Lifetimes.cpp`):
  - `make_hm3pt_hqe` became `make_hm3pt` with the flat parameters;
  - `hm3ptFourQuark` and `make_hm3pt_fourquark` were replaced by `make_weaknoneye_deriv`, with `deriv`/`dim7Specs` now
    returning `MContraction::WeakDerivSpec`;
  - the dim-7 measurement is one `WeakNonEye3ptDeriv` per probe pair, named `<q>h_DeltaHs2_dim7_<G5|GTG5>_<G5|GTG5>_t<t>_dt<Δt>_<τ>`
    (the probe tags moved from the entry `name` into the module name).
- Comments in `Inputs.hpp` and the README updated. No change to the XML input.
