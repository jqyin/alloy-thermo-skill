#include "lattice.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace alloy {

namespace {
[[noreturn]] void fail(const std::string& msg) {
    throw std::runtime_error("coupling.input: " + msg);
}
}  // namespace

Lattice Lattice::load(const std::string& coupling_path, int n_elements, int N) {
    Lattice lat;
    lat.ne_ = n_elements;
    lat.N_ = N;
    lat.n_sites_ = static_cast<long>(N) * N * N;

    std::ifstream f(coupling_path);
    if (!f) fail("cannot open '" + coupling_path + "'");

    // Each row lists one neighbor offset:
    //   r  [J ti tj]*npairs  dx dy dz
    // with npairs = ne*(ne+1)/2 covering every unordered (ti,tj) including
    // diagonals.  Rows are grouped by ascending shell distance r; a row with
    // r < 0 (or EOF) terminates the list.  The number of distinct distances
    // determines the number of coordination shells at run time.
    const int npairs = n_elements * (n_elements + 1) / 2;

    // Temporary per-shell J storage keyed by shell index, grown on demand.
    std::vector<double> Jflat;  // (shell, a, b) -> value, laid out shell-major
    auto ensure_shells = [&](int nsh) {
        long want = static_cast<long>(nsh) * n_elements * n_elements;
        if (static_cast<long>(Jflat.size()) < want) Jflat.resize(want, 0.0);
    };

    int shell = -1;
    double cur_dist = 0.0;
    constexpr double dist_tol = 1e-6;

    std::string line;
    while (std::getline(f, line)) {
        std::istringstream is(line);
        double r;
        if (!(is >> r)) continue;       // blank line
        if (r < 0) break;               // explicit terminator

        // Decide whether this row opens a new shell.
        if (shell < 0 || std::fabs(r - cur_dist) > dist_tol) {
            ++shell;
            cur_dist = r;
            lat.NS_.push_back(0);
            lat.dist_.push_back(r);
            ensure_shells(shell + 1);
        }

        // Read the npairs (value, ti, tj) triples.
        for (int p = 0; p < npairs; ++p) {
            double v;
            int ti, tj;
            if (!(is >> v >> ti >> tj))
                fail("row has fewer than " + std::to_string(npairs) +
                     " pair entries (expected for " +
                     std::to_string(n_elements) + " elements)");
            if (ti < 1 || ti > n_elements || tj < 1 || tj > n_elements)
                fail("species index out of range 1.." +
                     std::to_string(n_elements));
            int a = ti - 1, b = tj - 1;
            long base = static_cast<long>(shell) * n_elements * n_elements;
            Jflat[base + a * n_elements + b] = v;
            Jflat[base + b * n_elements + a] = v;  // enforce symmetry
        }

        int dx, dy, dz;
        if (!(is >> dx >> dy >> dz))
            fail("row missing offset vector (dx dy dz)");

        lat.offsets_.push_back({dx, dy, dz});
        lat.NS_[shell] += 1;
    }

    if (shell < 0) fail("no coupling rows found");

    lat.n_shells_ = shell + 1;
    lat.n_neighbors_ = static_cast<int>(lat.offsets_.size());

    // Repack J into the public (a,b,shell) layout.
    lat.J_.assign(static_cast<long>(n_elements) * n_elements * lat.n_shells_, 0.0);
    for (int s = 0; s < lat.n_shells_; ++s)
        for (int a = 0; a < n_elements; ++a)
            for (int b = 0; b < n_elements; ++b) {
                long src = static_cast<long>(s) * n_elements * n_elements +
                           a * n_elements + b;
                lat.J_[(static_cast<long>(a) * n_elements + b) * lat.n_shells_ + s] =
                    Jflat[src];
            }

    lat.build_neighbor_index();
    return lat;
}

void Lattice::build_neighbor_index() {
    auto wrap = [](int x, int n) { return ((x % n) + n) % n; };
    neighbor_index_.assign(n_sites_ * n_neighbors_, 0);
    const int N = N_, N2 = N_ * N_;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k) {
                long site = static_cast<long>(i) * N2 + j * N + k;
                for (int n = 0; n < n_neighbors_; ++n) {
                    int si = wrap(i + offsets_[n][0], N);
                    int sj = wrap(j + offsets_[n][1], N);
                    int sk = wrap(k + offsets_[n][2], N);
                    neighbor_index_[site * n_neighbors_ + n] =
                        static_cast<long>(si) * N2 + sj * N + sk;
                }
            }
}

int Lattice::neighbors_within(int n_shells_of) const {
    int total = 0;
    for (int s = 0; s < n_shells_of && s < n_shells_; ++s) total += NS_[s];
    return total;
}

void Lattice::print_summary() const {
    std::fprintf(stderr, "# lattice: %d shells, %d neighbors, %d^3 = %ld sites\n",
                 n_shells_, n_neighbors_, N_, n_sites_);
    std::fprintf(stderr, "# shell  count  distance   J[0][1]\n");
    for (int s = 0; s < n_shells_; ++s)
        std::fprintf(stderr, "#  %3d  %5d  %8.4f  %12.6g\n", s, NS_[s], dist_[s],
                     ne_ > 1 ? J(0, 1, s) : 0.0);
}

}  // namespace alloy
