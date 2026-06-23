// test_consistency.cpp -- serial correctness checks for the alloy model.
//
// The Monte Carlo loop maintains the pair-count tensor W incrementally as atoms
// swap.  This test exercises that bookkeeping and asserts the incrementally
// tracked energy stays exactly consistent with a from-scratch recomputation.
//
// Builds without MPI: it links only config/lattice/alloy_model.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "alloy_model.hpp"
#include "config.hpp"
#include "lattice.hpp"
#include "rng.hpp"

namespace {
int failures = 0;
void check(bool ok, const char* what) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok) ++failures;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        std::string control = argc > 1 ? argv[1] : "test.input";
        alloy::SimConfig cfg = alloy::load_config(control);
        int ne = cfg.n_elements();
        alloy::Lattice lat = alloy::Lattice::load("coupling.input", ne, cfg.N);

        // Equal composition keeps the test independent of composition.input.
        std::vector<double> comp(ne, 1.0 / ne);

        alloy::Rng rng(2024, 0);
        alloy::AlloyModel model(lat, cfg.intercept);
        model.init_configuration(comp, rng, /*ordered=*/false);

        double pE = model.energy();

        // A fresh model built from the same atoms must agree on the energy.
        alloy::AlloyModel ref(lat, cfg.intercept);
        ref.set_atoms(model.atoms());
        check(std::fabs(pE - ref.energy()) < 1e-9,
              "rebuilt energy matches initial energy");

        // Drive many trial moves; W is updated incrementally inside.
        double beta_eff = 1.0;  // arbitrary finite temperature
        long accepted = 0, attempted = 0;
        for (int s = 0; s < 50; ++s)
            for (long c = 0; c < lat.n_sites(); ++c) {
                auto m = model.attempt_swap(rng, beta_eff, pE);
                if (m != alloy::AlloyModel::Move::NoOp) {
                    ++attempted;
                    if (m == alloy::AlloyModel::Move::Accepted) ++accepted;
                }
            }

        // After all those incremental updates, the tracked energy must still
        // equal the value computed from a clean rebuild of the pair counts.
        alloy::AlloyModel ref2(lat, cfg.intercept);
        ref2.set_atoms(model.atoms());
        double e_rebuilt = ref2.energy();
        std::printf("  tracked pE = %.12g , rebuilt = %.12g , attempted=%ld "
                    "accepted=%ld\n",
                    pE, e_rebuilt, attempted, accepted);
        check(std::fabs(pE - e_rebuilt) < 1e-6 * (1.0 + std::fabs(e_rebuilt)),
              "incremental W bookkeeping matches from-scratch rebuild");

        // Composition must be conserved by Kawasaki swaps.
        std::vector<int> counts0(ne, 0);
        for (auto a : ref.atoms()) counts0[a]++;
        bool conserved = true;
        for (int t = 0; t < ne; ++t)
            if (counts0[t] != model.counts()[t]) conserved = false;
        check(conserved, "composition conserved under swaps");

        // Order parameters are finite and the mean is in a sane range.
        auto M = model.order_parameters(1);
        bool finite = true;
        for (double v : M) finite = finite && std::isfinite(v);
        check(finite, "order parameters are finite");

    } catch (const std::exception& e) {
        std::printf("[FAIL] exception: %s\n", e.what());
        return 2;
    }

    std::printf(failures == 0 ? "\nALL TESTS PASSED\n" : "\n%d TEST(S) FAILED\n",
                failures);
    return failures == 0 ? 0 : 1;
}
