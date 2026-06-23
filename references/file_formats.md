# File formats reference

All paths are relative to the **run directory** unless noted. The engine reads
`composition.input` and `coupling.input` by those fixed names, plus a control
file passed as `argv[1]`.

## spec.json (input to `make_inputs.py`)

A single JSON object. `make_inputs.py` turns it into the engine input files.

| key            | type                | meaning |
|----------------|---------------------|---------|
| `elements`     | list[str]           | species names; order defines indices 1..ne |
| `composition`  | list[float]         | fractions per element, sum ~ 1 |
| `lattice`      | `"bcc"\|"fcc"\|"sc"`| crystal lattice (sets neighbor geometry) |
| `n_shells`     | int                 | number of coordination shells to include |
| `energy_unit`  | `"ry"\|"ev"\|"mev"` | unit of the `couplings` values |
| `intercept`    | float               | constant energy offset / formula unit |
| `couplings`    | list of ne x ne     | one symmetric J matrix per shell, **nearest shell first** |
| `search`       | int (optional)      | offset search radius (default 4); raise for many shells |
| `simulation`   | object              | control parameters (see control file below) |

The `couplings[s][i][j]` value is `J` between element `i` and element `j` in
shell `s` (0-based shell, 0 = nearest). Matrices must be symmetric; the diagonal
is unused by the energy (set to 0). See `physics.md` for the sign convention.

Get a ready-to-edit template:

```bash
python scripts/make_inputs.py --print-template > spec.json
```

## composition.input

Whitespace-separated species fractions, one value per element, in element order:

```
0.36 0.22 0.22 0.2
```

## coupling.input

One row per neighbor offset. Rows are grouped by ascending shell distance; the
file ends with a row whose first field is negative (terminator). Each row:

```
r  J11 1 1  J12 1 2 ... J(ne)(ne) ne ne  dx dy dz
```

- `r` — neighbor distance in lattice-constant units (groups offsets into
  shells; **not** used in the energy).
- A `(value, ti, tj)` triple for every unordered pair `ti <= tj`, i.e.
  `ne*(ne+1)/2` triples, in the order `(1,1),(1,2),...,(1,ne),(2,2),...`.
  `ti`/`tj` are 1-based element indices.
- `dx dy dz` — integer neighbor offset in the primitive index basis.

The reader determines the number of shells at run time from the distinct `r`
values, and the number of neighbors per shell from how many rows share each `r`.
This is generated automatically; you normally never hand-write it.

## Control file (e.g. control.input)

Key/value lines; `#` comments and blank lines ignored. Unknown keys are an
error. Keys:

| key                | type   | meaning |
|--------------------|--------|---------|
| `elements`         | names  | ne species names (required) |
| `lattice`          | str    | informational label |
| `energy_unit`      | str    | `ry`/`ev`/`mev` (interprets `coupling.input`) |
| `intercept`        | float  | energy offset / formula unit |
| `N`                | int    | linear size; sites = N^3 |
| `T_init`,`T_final` | float  | ladder endpoints in K (geometric spacing) |
| `n_runs`           | int    | independent statistical runs |
| `n_drop`           | float  | equilibration sweeps before sampling |
| `n_samples`        | float  | number of measurement samples |
| `n_separation`     | float  | sweeps between samples |
| `init_order`       | 0/1    | 0 = disordered start, 1 = ordered/segregated start |
| `restart`          | 0/1    | 1 = resume from checkpoint files |
| `checkpoint_every` | int    | samples between wall-time checks (0 = never) |
| `walltime_hours`   | float  | soft budget; checkpoints near 0.9x and exits |

## Output: thermo_run<i>.csv

One file per independent run (`i = 0 .. n_runs-1`), CSV with a header row and
one data row per temperature replica. Columns:

```
T, E_per_site, C_per_site,
alpha_<el a>-<el b> (one per unordered pair),
alpha_mean, chi, binder, swap_accept, move_accept
```

See `physics.md` for definitions. `analyze.py` consumes these directly.

## Output: checkpoint_run<i>_rank<r>.bin

Binary per-rank checkpoint written when the wall-time budget is reached:
sample index, current energy, accumulators, and the species array. With
`restart 1` the engine reloads these and resumes the measurement phase. Do not
edit by hand.

## Output: manifest.json (from make_inputs.py)

Provenance: the resolved elements, composition, lattice, shell distances and
neighbor counts, energy unit, intercept, and the simulation block. Keep it with
the run for reproducibility.

## Command line

```
mpirun -n <nT> alloy_mc <control-file> [base_seed]
```

`<nT>` (MPI ranks) = number of temperature replicas. `base_seed` (default
12345) seeds a per-rank RNG stream so runs are reproducible from
`(base_seed, rank)`.
