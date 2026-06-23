// alloy_model.hpp -- substitutional-alloy chemical-ordering model.
//
// State is the species assignment Atom[site] in {0, .., ne-1} on the lattice.
// The energy is a pair-interaction (cluster-expansion) form
//
//     E = N^3 * ( sum_shell sum_{a<b} (W_ab^shell / NS_shell / N^3) * J_ab^shell
//                 + intercept )
//
// where W_ab^shell counts a-b neighbor pairs in a coordination shell.  Moves are
// composition-conserving Kawasaki swaps of unlike neighbors.  Because the energy
// depends only on the aggregate pair counts W, a full energy evaluation is
// O(ne^2 * n_shells) and independent of system size, so each trial move
// recomputes the total energy exactly after an incremental update of W.
#ifndef ALLOY_MC_ALLOY_MODEL_HPP
#define ALLOY_MC_ALLOY_MODEL_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "lattice.hpp"
#include "rng.hpp"

namespace alloy {

class AlloyModel {
public:
    AlloyModel(const Lattice& lat, double intercept);

    // Populate the lattice from target composition fractions (size ne).  When
    // `ordered` is false the species are randomly shuffled; when true they are
    // laid down in contiguous blocks (a segregated initial state).
    void init_configuration(const std::vector<double>& composition, Rng& rng,
                             bool ordered);

    // Set the species array directly (used on replica exchange / restart) and
    // rebuild all derived quantities.
    void set_atoms(const std::vector<std::int16_t>& atoms);
    const std::vector<std::int16_t>& atoms() const { return atom_; }
    std::vector<std::int16_t>& atoms() { return atom_; }
    const std::vector<int>& counts() const { return nt_; }

    void rebuild_pair_counts();          // recompute W from scratch
    double energy() const;               // total energy in J's native unit

    // Outcome of a trial move, so callers can separate genuine attempts (an
    // unlike pair) from no-ops when reporting acceptance ratios.
    enum class Move { NoOp = -1, Rejected = 0, Accepted = 1 };

    // Attempt one Kawasaki swap.  `beta_eff` = E_scale / (kB * T) in the native
    // energy unit so that the Metropolis factor is exp(-dE * beta_eff).
    // Updates `current_energy` in place.
    Move attempt_swap(Rng& rng, double beta_eff, double& current_energy);

    // Warren-Cowley-type short-range-order parameters alpha_ab evaluated over
    // the first `order_shells` shells.  Returns a vector of length
    // ne*(ne-1)/2 + 1; the final element is the mean of |alpha| over all pairs.
    std::vector<double> order_parameters(int order_shells = 1) const;

    // Snapshot output: "<species> i j k" per site (cubic index coordinates).
    void write_xyz(const std::string& path, const std::vector<std::string>& names,
                   int frame, double T, double energy_per_site) const;

private:
    const Lattice& lat_;
    int ne_;
    double intercept_;

    std::vector<std::int16_t> atom_;     // [n_sites] species id
    std::vector<int> nt_;                // [ne] species counts
    std::vector<long long> W_;           // [ne*ne*n_shells] pair counts (a<=b)
    std::vector<long long> W_backup_;    // reusable rollback buffer

    long long& Wref(int a, int b, int s) {
        return W_[(static_cast<long>(a) * ne_ + b) * lat_.n_shells() + s];
    }
    long long Wat(int a, int b, int s) const {
        return W_[(static_cast<long>(a) * ne_ + b) * lat_.n_shells() + s];
    }

    // Incremental W update when site `pi` changes species ai->aj (its partner
    // pj is skipped because the pi-pj bond is handled by the paired call).
    void update_pair_counts_site(int ai, int aj, long pi, long pj);
};

}  // namespace alloy

#endif  // ALLOY_MC_ALLOY_MODEL_HPP
