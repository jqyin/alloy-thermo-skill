# Worked example: MoNbTaW refractory high-entropy alloy

`spec.json` is a complete specification for the equimolar-ish BCC refractory HEA
**Mo-Nb-Ta-W** (composition 0.36 / 0.22 / 0.22 / 0.20), with **6 coordination
shells** of DFT-derived effective pair interactions (Rydberg units). The shell
geometry is BCC (coordination 8, 6, 12, 24, 8, 6).

## Run it end to end

```bash
# 1. Build the engine (load your MPI module first on a cluster)
cd ../../engine && make && cd -

# 2. Generate input files into a fresh run directory
python ../../scripts/make_inputs.py spec.json --outdir run/ --control-name control.input

# 3a. Quick local test (8 replicas, reduce N/samples in spec.json first)
cd run && mpirun -n 8 ../../../engine/alloy_mc control.input 2024 && cd ..

# 3b. Or submit to SLURM with one rank per temperature replica
python ../../scripts/submit.py --rundir run/ --binary ../../engine/alloy_mc \
    --control control.input --account <PROJECT> --ntasks 64 \
    --tasks-per-node 56 --time 02:00:00 \
    --module "module load PrgEnv-gnu cray-mpich" --submit

# 4. Analyze: plots + estimated transition temperature
python ../../scripts/analyze.py --dir run/ --out run/analysis
cat run/analysis/summary.md
```

## What to expect

On cooling, the unlike-pair interactions (negative `J` in the near shells)
drive **B2/B32-like chemical ordering**: the Warren-Cowley parameters grow in
magnitude and the specific heat develops a peak at the order-disorder
transition. The exact `Tc` depends on system size and sampling — see the
caveats in `references/physics.md` and use several `N` values for a rigorous
estimate.

## Notes

- `spec.json` requests `N=16`, `n_runs=4`, `n_samples=80000`. For a laptop smoke
  test, lower these (e.g. `N=8`, `n_drop=2000`, `n_samples=4000`, `n_runs=1`).
- The `couplings` here were extracted directly from the original validated
  `coupling.input`, so the generated geometry is identical to the legacy file.
