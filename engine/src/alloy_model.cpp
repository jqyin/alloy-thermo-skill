#include "alloy_model.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <stdexcept>

namespace alloy {

AlloyModel::AlloyModel(const Lattice& lat, double intercept)
    : lat_(lat), ne_(lat.n_elements()), intercept_(intercept) {
    atom_.assign(lat_.n_sites(), 0);
    nt_.assign(ne_, 0);
    W_.assign(static_cast<long>(ne_) * ne_ * lat_.n_shells(), 0);
    W_backup_.resize(W_.size());
}

void AlloyModel::init_configuration(const std::vector<double>& composition,
                                    Rng& rng, bool ordered) {
    if (static_cast<int>(composition.size()) != ne_)
        throw std::runtime_error("composition size != number of elements");

    const long n = lat_.n_sites();
    std::vector<int> list;
    list.reserve(n);
    // Allocate floor(n * c_t) sites to each species, then fill the remainder
    // with random species (mirrors the legacy rounding behavior).
    for (int t = 0; t < ne_; ++t) {
        long ni = static_cast<long>(n * composition[t]);
        for (long c = 0; c < ni; ++c) list.push_back(t);
    }
    while (static_cast<long>(list.size()) < n)
        list.push_back(static_cast<int>(rng.below(ne_)));

    if (!ordered) {
        // Fisher-Yates shuffle using the shared RNG.
        for (long i = n - 1; i > 0; --i) {
            long j = static_cast<long>(rng.below(static_cast<std::size_t>(i + 1)));
            std::swap(list[i], list[j]);
        }
    }

    std::fill(nt_.begin(), nt_.end(), 0);
    for (long s = 0; s < n; ++s) {
        atom_[s] = static_cast<std::int16_t>(list[s]);
        nt_[list[s]]++;
    }
    rebuild_pair_counts();
}

void AlloyModel::set_atoms(const std::vector<std::int16_t>& atoms) {
    if (atoms.size() != atom_.size())
        throw std::runtime_error("set_atoms: size mismatch");
    atom_ = atoms;
    std::fill(nt_.begin(), nt_.end(), 0);
    for (auto a : atom_) nt_[a]++;
    rebuild_pair_counts();
}

void AlloyModel::rebuild_pair_counts() {
    std::fill(W_.begin(), W_.end(), 0);
    const int nn = lat_.n_neighbors();
    for (long site = 0; site < lat_.n_sites(); ++site) {
        int ai = atom_[site];
        int cnt = 0;
        for (int shell = 0; shell < lat_.n_shells(); ++shell) {
            for (int ii = 0; ii < lat_.shell_count(shell); ++ii) {
                int aj = atom_[lat_.neighbor(site, cnt++)];
                if (ai <= aj)
                    Wref(ai, aj, shell)++;
                else
                    Wref(aj, ai, shell)++;
            }
        }
        (void)nn;
    }
}

double AlloyModel::energy() const {
    const long n3 = lat_.n_sites();
    double e = 0.0;
    for (int shell = 0; shell < lat_.n_shells(); ++shell) {
        double inv = 1.0 / (static_cast<double>(lat_.shell_count(shell)) * n3);
        for (int a = 0; a < ne_ - 1; ++a)
            for (int b = a + 1; b < ne_; ++b)
                e += (Wat(a, b, shell) * inv) * lat_.J(a, b, shell);
    }
    e += intercept_;
    return e * n3;
}

void AlloyModel::update_pair_counts_site(int ai, int aj, long pi, long pj) {
    int cnt = 0;
    for (int shell = 0; shell < lat_.n_shells(); ++shell) {
        for (int j = 0; j < lat_.shell_count(shell); ++j) {
            long id = lat_.neighbor(pi, cnt++);
            if (id == pj) continue;       // the pi-pj bond is handled separately
            int a = atom_[id];
            // Each bond is counted from both endpoints, hence the factor of 2.
            if (ai <= a) Wref(ai, a, shell) -= 2; else Wref(a, ai, shell) -= 2;
            if (aj <= a) Wref(aj, a, shell) += 2; else Wref(a, aj, shell) += 2;
        }
    }
}

AlloyModel::Move AlloyModel::attempt_swap(Rng& rng, double beta_eff,
                                          double& current_energy) {
    long i = static_cast<long>(rng.below(static_cast<std::size_t>(lat_.n_sites())));
    int n = static_cast<int>(rng.below(static_cast<std::size_t>(lat_.n_neighbors())));
    long j = lat_.neighbor(i, n);

    int ai = atom_[i], aj = atom_[j];
    if (ai == aj) return Move::NoOp;      // swapping like atoms changes nothing

    double e_before = current_energy;
    W_backup_ = W_;                       // cheap rollback snapshot (reused buffer)

    update_pair_counts_site(ai, aj, i, j);
    update_pair_counts_site(aj, ai, j, i);
    atom_[i] = static_cast<std::int16_t>(aj);
    atom_[j] = static_cast<std::int16_t>(ai);

    double e_after = energy();
    double dE = e_after - e_before;

    bool accept = (dE <= 0.0) || (rng.uniform() < std::exp(-dE * beta_eff));
    if (accept) {
        current_energy = e_after;
        return Move::Accepted;
    }
    atom_[i] = static_cast<std::int16_t>(ai);
    atom_[j] = static_cast<std::int16_t>(aj);
    W_.swap(W_backup_);                   // restore previous counts
    return Move::Rejected;
}

std::vector<double> AlloyModel::order_parameters(int order_shells) const {
    const long n3 = lat_.n_sites();
    std::vector<double> occ(ne_);
    for (int t = 0; t < ne_; ++t) occ[t] = static_cast<double>(nt_[t]) / n3;

    int neighbors = lat_.neighbors_within(order_shells);
    if (neighbors == 0) neighbors = 1;
    std::vector<double> op(static_cast<long>(ne_) * ne_, 0.0);

    for (long site = 0; site < n3; ++site) {
        int ai = atom_[site];
        for (int k = 0; k < neighbors; ++k) {
            int aj = atom_[lat_.neighbor(site, k)];
            op[ai * ne_ + aj] += 1.0 / neighbors;
        }
    }
    for (int t = 0; t < ne_; ++t)
        for (int it = 0; it < ne_; ++it) {
            double v = op[t * ne_ + it] / n3;
            op[t * ne_ + it] =
                (occ[t] > 0 && occ[it] > 0) ? 1.0 - v / (occ[t] * occ[it]) : 0.0;
        }

    int npair = ne_ * (ne_ - 1) / 2;
    std::vector<double> M(npair + 1, 0.0);
    int c = 0;
    for (int a = 0; a < ne_ - 1; ++a)
        for (int b = a + 1; b < ne_; ++b) M[c++] = op[a * ne_ + b];
    double mean = 0.0;
    for (int p = 0; p < npair; ++p) mean += std::fabs(M[p]);
    M[npair] = npair > 0 ? mean / npair : 0.0;
    return M;
}

void AlloyModel::write_xyz(const std::string& path,
                           const std::vector<std::string>& names, int frame,
                           double T, double energy_per_site) const {
    std::ofstream f(path, std::ios::app);
    if (!f) return;
    f << lat_.n_sites() << "\n";
    f << "Energy_per_site=" << energy_per_site << " T=" << T
      << " frame=" << frame << "\n";
    const int N = lat_.N(), N2 = N * N;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) {
                long s = static_cast<long>(i) * N2 + j * N + k;
                f << names[atom_[s]] << ' ' << i << ' ' << j << ' ' << k << "\n";
            }
}

}  // namespace alloy
