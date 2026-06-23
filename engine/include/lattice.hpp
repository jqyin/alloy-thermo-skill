// lattice.hpp -- lattice geometry and pair couplings.
//
// The geometry is fully data-driven: the set of neighbor offset vectors (and
// hence the lattice type -- BCC, FCC, SC, ...) is read from coupling.input,
// together with the per-shell symmetric pair-interaction matrices J_ij.  Sites
// live on a simple-cubic index grid of N x N x N points with periodic boundary
// conditions; the physical lattice is encoded purely in the offset vectors.
#ifndef ALLOY_MC_LATTICE_HPP
#define ALLOY_MC_LATTICE_HPP

#include <array>
#include <string>
#include <vector>

namespace alloy {

class Lattice {
public:
    // Read coupling.input for `n_elements` species and build the neighbor-index
    // table for an N x N x N grid.  Throws std::runtime_error on malformed data.
    static Lattice load(const std::string& coupling_path, int n_elements, int N);

    int N() const { return N_; }
    long n_sites() const { return n_sites_; }
    int n_shells() const { return n_shells_; }
    int n_neighbors() const { return n_neighbors_; }
    int n_elements() const { return ne_; }

    int shell_count(int shell) const { return NS_[shell]; }   // neighbors in shell
    double shell_dist(int shell) const { return dist_[shell]; }

    // Symmetric pair coupling J between species a and b in coordination shell.
    double J(int a, int b, int shell) const {
        return J_[(static_cast<long>(a) * ne_ + b) * n_shells_ + shell];
    }

    // Global flat index of the k-th neighbor of `site`.
    long neighbor(long site, int k) const {
        return neighbor_index_[site * n_neighbors_ + k];
    }

    // Total neighbors contained in shells [0, n_shells_of) -- used for the
    // short-range-order parameter, which is evaluated over the first shell(s).
    int neighbors_within(int n_shells_of) const;

    void print_summary() const;  // rank-0 diagnostic of shells and J

private:
    int N_ = 0;
    long n_sites_ = 0;
    int ne_ = 0;
    int n_shells_ = 0;
    int n_neighbors_ = 0;

    std::vector<int> NS_;             // [n_shells] neighbors per shell
    std::vector<double> dist_;        // [n_shells] representative distance
    std::vector<std::array<int, 3>> offsets_;  // [n_neighbors]
    std::vector<double> J_;           // [ne*ne*n_shells], symmetric in (a,b)
    std::vector<long> neighbor_index_;  // [n_sites * n_neighbors]

    void build_neighbor_index();
};

}  // namespace alloy

#endif  // ALLOY_MC_LATTICE_HPP
