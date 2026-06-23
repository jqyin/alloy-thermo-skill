---
name: alloy-thermo-skill
description: >-
  Run Monte Carlo thermodynamics simulations of multicomponent (high-entropy)
  alloys on HPC and interpret the results. Use this when the task involves
  predicting the order-disorder transition temperature, specific heat,
  short-range order, or thermodynamic stability of a substitutional alloy from
  DFT-derived pair (cluster-expansion) interactions. The skill builds a C++/MPI
  replica-exchange (parallel-tempering) engine, prepares its input files from a
  composition and per-shell J_ij coupling matrices, submits SLURM jobs, and
  plots specific heat / susceptibility / Warren-Cowley order parameters to
  locate phase transitions. Works for any number of elements and BCC/FCC/SC
  lattices. Trigger on requests about alloy ordering, HEA thermodynamics, Tc /
  order-disorder transitions, lattice Monte Carlo, or cluster-expansion MC.
license: MIT
---

# Alloy thermodynamics by parallel-tempering Monte Carlo

This skill drives a lattice Monte Carlo engine that computes the finite-
temperature thermodynamics of a substitutional alloy whose energy is a
**pair-interaction (cluster-expansion) Hamiltonian** parameterized by
DFT-derived effective pair interactions `J_ij`. It is built for high-entropy
alloys (HEAs) but works for any number of species (>= 2) on BCC, FCC, or SC
lattices.

The engine uses **replica-exchange / parallel tempering**: one MPI rank holds
one temperature, and configurations swap between adjacent temperatures to
equilibrate efficiently across an order-disorder transition. Moves are
composition-conserving Kawasaki atom swaps. Outputs are energy, specific heat,
Warren-Cowley short-range-order (SRO) parameters, susceptibility, and Binder
cumulant versus temperature.

## When to use

- "Find the order-disorder transition temperature of <alloy>."
- "Compute the specific heat / short-range order of a high-entropy alloy from
  these pair interactions."
- "Run lattice/cluster-expansion Monte Carlo for <composition> on HPC and tell
  me if it orders."

## Repository layout

```
engine/        C++17 MPI engine (the simulator)
  include/ src/   modular sources; build with make or cmake
  tests/          serial correctness self-test
scripts/
  make_inputs.py  build composition.input + coupling.input + control file
  submit.py       render/submit a SLURM job
  analyze.py      plot thermodynamics and estimate Tc
assets/
  job.slurm.template
references/
  physics.md      model, units, observables, Tc detection, caveats
  file_formats.md  exact input/output file formats
  hpc.md          building and running on SLURM clusters
examples/MoNbTaW/  a worked 4-element BCC example (spec.json)
```

## End-to-end workflow

Follow these steps. Read `references/` files as needed — do not guess file
formats; they are documented exactly.

### 1. Build the engine

On a cluster, load an MPI compiler first (see `references/hpc.md` for examples),
then:

```bash
cd engine && make            # produces engine/alloy_mc
make test                    # optional: build the serial self-test
```

`make` defaults to the `mpic++` wrapper; override with `make CXX=mpicxx` or
`make CXX=CC` (Cray). No external libraries are required beyond MPI and a C++17
compiler.

### 2. Prepare inputs from DFT data

Write a JSON **spec** describing the chemistry and the pair interactions. Get a
template with:

```bash
python scripts/make_inputs.py --print-template > spec.json
```

Edit `spec.json` to set: `elements`, `composition` (fractions, sum ~ 1),
`lattice` (`bcc`/`fcc`/`sc`), `n_shells`, `energy_unit` (`ry`/`ev`/`mev`),
`intercept`, the per-shell symmetric `couplings` matrices `J_ij` (one
`n_elements x n_elements` matrix per shell, **ordered from nearest shell
outward**), and the `simulation` block (`N`, `T_init`, `T_final`, `n_runs`,
`n_drop`, `n_samples`, ...). See `references/physics.md` for the **sign
convention** of `J` and `references/file_formats.md` for the spec schema.

Then generate the run directory:

```bash
python scripts/make_inputs.py spec.json --outdir run/ --control-name control.input
```

This emits `run/composition.input`, `run/coupling.input`, `run/control.input`,
and `run/manifest.json` (provenance). The script prints the shell distances and
neighbor counts — sanity-check these against the expected lattice coordination
(BCC: 8,6,12,...; FCC: 12,6,24,...).

### 3. Submit to HPC

The **number of MPI ranks equals the number of temperature replicas** in the
parallel-tempering ladder (geometrically spaced between `T_init` and
`T_final`). Choose enough replicas to (a) bracket the expected transition and
(b) keep adjacent-temperature swap acceptance healthy (~20-40%). 32-128 replicas
is typical.

```bash
python scripts/submit.py --rundir run/ --binary engine/alloy_mc \
    --control control.input --account <PROJECT> --ntasks 64 \
    --tasks-per-node 56 --time 02:00:00 \
    --module "module load PrgEnv-gnu" --submit
```

Omit `--submit` to only render `run/job.slurm` for inspection. For local
testing without SLURM:

```bash
cd run && mpirun -n 8 ../engine/alloy_mc control.input 12345
```

### 4. Monitor

Use `squeue --me` / `sacct`. The engine checkpoints near the wall-time budget
(`walltime_hours`, `checkpoint_every`) and can resume by setting `restart 1` in
the control file and resubmitting (see `references/hpc.md`).

### 5. Analyze and interpret

When the job finishes, `run/` contains `thermo_run<i>.csv` (one per independent
run). Analyze:

```bash
python scripts/analyze.py --dir run/ --out run/analysis
```

This averages over runs, writes `thermo.png` (energy, specific heat,
susceptibility, Binder), `order.png` (SRO parameters), and `summary.{json,md}`
with the estimated transition temperature from the **specific-heat peak** and
the **susceptibility peak**. Read `summary.md` and report:

- the estimated Tc and whether the two estimators agree,
- the SRO trend (large |alpha| at low T -> ordering; alpha -> 0 at high T ->
  random solid solution),
- caveats (see below).

### 6. Iterate

If the analysis flags problems, adjust and rerun:

- **Peak not bracketed / off-grid:** widen or re-center `T_init`/`T_final`, or
  add replicas to densify the ladder near the peak.
- **Low swap acceptance:** add replicas (finer ladder) or narrow the range.
- **Noisy / unconverged:** increase `n_drop` (equilibration) and `n_samples`,
  or `n_runs` for error bars.
- **Pinning Tc precisely:** repeat at several `N` and compare specific-heat
  peak positions or Binder-cumulant crossings (finite-size scaling).

## Key guidance

- `N` is the linear lattice size; total sites = `N^3`. Start at `N=12-16` for
  exploration; go larger (24-32) to sharpen the transition. Cost scales with
  `N^3 * n_samples`.
- One rank per temperature: more ranks = finer ladder = better sampling, more
  cost. Cost is independent of the number of shells for the energy evaluation.
- Always verify composition is conserved and acceptance ratios are sane
  (`swap_accept`, `move_accept` columns in the CSV).

See `references/physics.md` for the full model and `references/file_formats.md`
for every file format. A complete worked example is in `examples/MoNbTaW/`.
