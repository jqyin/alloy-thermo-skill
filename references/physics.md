# Physics and model reference

## What is being simulated

A **substitutional alloy** on a fixed lattice: every lattice site is occupied by
one of `ne` chemical species, and the only degrees of freedom are *which*
species sits on each site. The composition (fraction of each species) is fixed.
The simulation samples the equilibrium distribution of chemical arrangements as
a function of temperature to reveal **chemical ordering** and the
**order-disorder transition**.

This is the standard lattice model behind cluster-expansion + Monte Carlo
studies of high-entropy alloys.

## Lattice representation

Sites live on a simple-cubic index grid `N x N x N` with periodic boundary
conditions (total `N^3` sites). The physical crystal (BCC/FCC/SC) is encoded
entirely by the set of **neighbor offset vectors** in `coupling.input`. The
input generator (`scripts/make_inputs.py`) produces these offsets by enumerating
integer displacements in the crystal's *primitive* basis, so the cubic index
grid faithfully represents the chosen lattice:

- BCC nearest-neighbor shells have coordination 8, 6, 12, 24, 8, 6, ...
- FCC shells have coordination 12, 6, 24, 12, ...
- SC shells have coordination 6, 12, 8, ...

Neighbor "distances" in `coupling.input` are in units of the lattice constant
`a` and are used only to group offsets into coordination shells (not in the
energy).

## Hamiltonian

The total energy is a pair-interaction (two-body cluster expansion):

```
E = N^3 * ( sum_shell sum_{a<b} (W_ab^shell / NS_shell / N^3) * J_ab^shell
            + intercept )
```

where:

- `J_ab^shell` is the **effective pair interaction** between species `a` and `b`
  in a given coordination `shell` (symmetric: `J_ab = J_ba`), supplied from DFT.
- `W_ab^shell` is the number of `a-b` neighbor pairs in that shell (each bond
  counted from both endpoints, i.e. twice).
- `NS_shell` is the number of neighbors per site in that shell.
- `intercept` is a constant energy offset per formula unit (does not affect any
  transition temperature — it only shifts the absolute energy).

Only **unlike** pairs (`a < b`, off-diagonal `J`) contribute to the energy; the
diagonal `J_aa` are read for completeness but unused. Because the energy depends
only on the aggregate counts `W`, a full energy evaluation is `O(ne^2 *
n_shells)` and **independent of system size**, so each trial move recomputes the
exact total energy after an incremental `W` update. (The serial self-test in
`engine/tests/` verifies this incremental bookkeeping against a from-scratch
rebuild.)

### Sign convention of J (important)

The energy adds `+W_ab * J_ab`. With `W_ab >= 0`, a **negative** `J_ab` for an
unlike pair *lowers* the energy when there are more `a-b` neighbors, i.e. it
**favors ordering / mixing** (unlike neighbors preferred). A **positive** `J_ab`
favors **clustering / segregation** (like neighbors preferred). Feed DFT
effective pair interactions with this convention. If your ECIs use the opposite
sign convention, negate them.

## Units

- `J_ij` are read in `energy_unit` = `ry` (Rydberg), `ev`, or `mev`, and
  converted internally to meV (`1 Ry = 13605.69 meV`, `1 eV = 1000 meV`).
- Temperatures are in **Kelvin**; converted to energy via Boltzmann's constant
  `kB = 8.617333e-2 meV/K`.
- Acceptance factors are therefore the dimensionless `exp(-dE_meV / (kB*T))`.

## Monte Carlo dynamics

- **Single-site moves:** a Kawasaki swap picks a random site and a random
  neighbor; if they are unlike species it proposes exchanging them and accepts
  with the Metropolis probability `min(1, exp(-dE/(kB T)))`. This conserves the
  global composition exactly. One *sweep* attempts `N^3` such moves.
- **Replica exchange (parallel tempering):** every couple of sweeps, adjacent
  temperature replicas attempt to swap their entire configurations with
  probability `min(1, exp(-(beta_i - beta_j)(E_j - E_i)))`. This lets
  configurations diffuse in temperature and escape low-T metastable states,
  which is essential for sampling near an ordering transition.

## Observables (columns of `thermo_run<i>.csv`)

- `E_per_site` — mean energy per site (native energy unit).
- `C_per_site` — specific heat per site, in units of `kB`:
  `C = Var(E) * E_scale^2 / (kB T)^2 / N^3`, from energy fluctuations.
- `alpha_<A>-<B>` — **Warren-Cowley short-range-order** parameter for the pair
  A-B over the first coordination shell:
  `alpha_AB = 1 - P_AB / (c_A c_B)`, where `P_AB` is the average fraction of a
  site's first-shell neighbors of type B given the site is type A. `alpha = 0`
  means random (ideal solid solution); `alpha < 0` means A-B pairs are favored
  (ordering tendency); `alpha > 0` means A-B avoided (clustering).
- `alpha_mean` — mean of `|alpha_AB|` over all pairs; a scalar "amount of order".
- `chi` — susceptibility of `alpha_mean`: `Var(alpha_mean) * N^3 / (kB T)`.
- `binder` — Binder cumulant of `alpha_mean`: `1 - <m^4>/(3<m^2>^2)`.
- `swap_accept`, `move_accept` — replica-exchange and single-move acceptance
  ratios (diagnostics; not physics).

## Locating the transition (Tc)

The order-disorder transition shows up as:

1. A **peak in the specific heat** `C_per_site(T)` (latent-heat-like fluctuation
   peak). `analyze.py` reports its position with a local quadratic refinement.
2. A **peak in the susceptibility** `chi(T)` of the order parameter.
3. An inflection / drop in the **Binder cumulant**; for a true continuous
   transition, Binder curves for different `N` **cross** at Tc.

`analyze.py` reports Tc from estimators (1) and (2) and flags whether they
agree. They should roughly coincide for a well-resolved transition.

### Caveats

- **Finite size:** a single `N` cannot pin Tc exactly. The specific-heat peak
  shifts and sharpens as `N` grows. For a rigorous Tc, run several sizes and use
  Binder-cumulant crossings or finite-size scaling of the peak position.
- **Equilibration:** ensure `n_drop` is large enough and replica-exchange
  acceptance is healthy (~20-40%), otherwise low-T replicas may be stuck and the
  apparent ordering is a sampling artifact.
- **Ladder coverage:** the temperature ladder must bracket the peak; a peak at
  the edge of the range means you should extend or re-center `T_init`/`T_final`.
- **No transition:** a broad, low specific-heat hump with `alpha` smoothly
  decaying can indicate only short-range order with no sharp transition (common
  for many HEAs). Report this honestly rather than forcing a Tc.
