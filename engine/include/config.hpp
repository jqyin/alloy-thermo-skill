// config.hpp -- simulation control parameters.
//
// All run-time knobs live in a single key/value text file (see
// references/file_formats.md).  Unlike the legacy build, nothing about the
// chemistry (number of elements, their names, energy units) is fixed at compile
// time -- it is all parsed here.
#ifndef ALLOY_MC_CONFIG_HPP
#define ALLOY_MC_CONFIG_HPP

#include <string>
#include <vector>

namespace alloy {

struct SimConfig {
    // --- Chemistry -------------------------------------------------------
    std::vector<std::string> elements;   // species names; size == n_elements
    std::string energy_unit = "ry";      // unit of J_ij in coupling.input
    double intercept = 0.0;              // constant energy offset / formula unit

    // --- Lattice ---------------------------------------------------------
    int N = 8;                           // linear size; total sites = N^3
    std::string lattice = "bcc";         // informational label only

    // --- Temperature ladder (Kelvin) ------------------------------------
    double T_init  = 100.0;              // first replica temperature
    double T_final = 2000.0;             // last replica temperature

    // --- Monte Carlo schedule -------------------------------------------
    int    n_runs       = 1;             // independent statistical runs
    double n_drop       = 1.0e4;         // equilibration sweeps
    double n_samples    = 5.0e4;         // measurement samples
    double n_separation = 1.0;           // sweeps between measurements

    // --- Misc ------------------------------------------------------------
    int    init_order      = 0;          // 0 = disordered start, 1 = ordered
    int    restart         = 0;          // 0 = fresh, 1 = resume checkpoint
    int    checkpoint_every = 1000;      // sweeps between wall-time checks
    double walltime_hours   = 48.0;      // soft wall-time budget

    int n_elements() const { return static_cast<int>(elements.size()); }
};

// Parse a simulation control file.  Throws std::runtime_error on malformed
// input or out-of-range values.
SimConfig load_config(const std::string& path);

// Echo the resolved configuration (rank 0 only, typically) for the log.
void print_config(const SimConfig& cfg);

}  // namespace alloy

#endif  // ALLOY_MC_CONFIG_HPP
