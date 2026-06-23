# Building and running on HPC (SLURM)

## Building

The engine needs only an MPI C++ compiler and C++17. Load your MPI stack, then
build:

```bash
cd engine
make                 # uses mpic++ by default -> engine/alloy_mc
# overrides:
make CXX=mpicxx      # OpenMPI/MPICH wrapper
make CXX=CC          # Cray programming environment wrapper
```

Module examples (adapt to your machine):

```bash
# Generic GNU + OpenMPI
module load gcc openmpi

# Cray (e.g. ORNL Frontier / Andes-like)
module load PrgEnv-gnu cray-mpich
make CXX=CC
```

CMake alternative:

```bash
cd engine
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j      # -> engine/build/alloy_mc
```

## Ranks, replicas, and nodes

The engine runs **one MPI rank per temperature replica**. The temperature ladder
is generated automatically as a geometric sequence from `T_init` to `T_final`
with `nranks` points. Therefore:

- More ranks => finer ladder => better replica-exchange mixing and more
  temperature resolution for locating Tc.
- Aim for adjacent-temperature swap acceptance ~20-40% (see `swap_accept` in the
  output CSV). If it is too low, add ranks or narrow `[T_init, T_final]`.
- 32-128 replicas is a typical sweet spot for an exploratory study.

`submit.py` computes nodes = ceil(ntasks / tasks_per_node). Set
`--tasks-per-node` to your machine's cores-per-node (or the value that fits
memory).

## Submitting

```bash
python scripts/submit.py \
    --rundir run/ --binary engine/alloy_mc --control control.input \
    --account <PROJECT> --ntasks 64 --tasks-per-node 56 --time 02:00:00 \
    --partition <queue> \
    --module "module load PrgEnv-gnu cray-mpich" \
    --submit
```

- Omit `--submit` to only write `run/job.slurm`; inspect, then `sbatch` it
  yourself.
- Repeat `--module` for multiple setup lines; they are written verbatim into the
  script's environment section.
- The script warns if `coupling.input`, `composition.input`, or the control file
  are missing from the run directory.

The rendered script `cd`s into the run directory and calls
`srun -n <ntasks> alloy_mc <control> <seed>`.

## Memory and cost

- Memory per rank is dominated by the neighbor-index table: roughly
  `N^3 * n_neighbors * 8 bytes`. Example: `N=32`, 64 neighbors -> ~16 MB/rank.
  This is small; you are compute-bound, not memory-bound.
- Wall-clock scales as `n_runs * (n_drop + n_samples * n_separation) * N^3`.
  The per-move energy evaluation is independent of `N` and of the number of
  shells, so the dominant factor is the number of sweeps times `N^3`.

## Checkpoint / restart

For long runs, set a realistic `walltime_hours` (a little under the SLURM
`-t`) and a nonzero `checkpoint_every`. When the soft budget is reached the
engine writes `checkpoint_run<i>_rank<r>.bin` and exits cleanly. To resume:

1. Set `restart 1` in the control file.
2. Resubmit the same job in the same run directory.

The engine reloads each rank's configuration and accumulators and continues the
measurement phase. (Restart resumes the current run's sampling; for multi-run
campaigns, prefer sizing each run to fit one allocation.)

## Local testing without SLURM

```bash
cd run
mpirun -n 8 ../engine/alloy_mc control.input 12345
# On a laptop with fewer cores than ranks, add --oversubscribe (OpenMPI).
```

Use a small `N` and short `n_drop`/`n_samples` for a smoke test before scaling
up.

## Troubleshooting

- **`mpi.h not found`:** the build used a non-MPI compiler. Ensure `mpic++`/`CC`
  resolves after `module load`, or pass `make CXX=...`.
- **All ranks abort with a parse error:** check the control file and that
  `composition.input` has exactly `ne` values matching `elements`.
- **Neighbor counts look wrong:** verify `lattice` and `n_shells`; the generator
  prints shell distances and coordination numbers — compare to the known
  coordination sequence for that lattice.
- **Specific-heat peak at the edge of the ladder:** extend or re-center
  `[T_init, T_final]` and rerun.
- **Results not reproducible:** reproducibility is per `(base_seed, rank)` and
  per rank count; changing the number of ranks changes the ladder and streams.
