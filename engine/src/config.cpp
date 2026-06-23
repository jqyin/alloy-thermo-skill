#include "config.hpp"

#include <cctype>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace alloy {

namespace {

// Trim leading/trailing whitespace.
std::string trim(const std::string& s) {
    std::size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    std::size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

[[noreturn]] void fail(const std::string& msg) {
    throw std::runtime_error("config: " + msg);
}

}  // namespace

SimConfig load_config(const std::string& path) {
    std::ifstream f(path);
    if (!f) fail("cannot open control file '" + path + "'");

    SimConfig cfg;
    bool have_elements = false;
    std::string line;
    while (std::getline(f, line)) {
        std::string t = trim(line);
        if (t.empty() || t[0] == '#') continue;

        std::istringstream is(t);
        std::string key;
        is >> key;

        if (key == "elements") {
            cfg.elements.clear();
            std::string name;
            while (is >> name) cfg.elements.push_back(name);
            if (cfg.elements.empty()) fail("'elements' needs >= 1 name");
            have_elements = true;
            continue;  // consuming names to EOL sets failbit; skip the check
        } else if (key == "energy_unit") {
            is >> cfg.energy_unit;
        } else if (key == "intercept") {
            is >> cfg.intercept;
        } else if (key == "lattice") {
            is >> cfg.lattice;
        } else if (key == "N") {
            is >> cfg.N;
        } else if (key == "T_init") {
            is >> cfg.T_init;
        } else if (key == "T_final") {
            is >> cfg.T_final;
        } else if (key == "n_runs") {
            is >> cfg.n_runs;
        } else if (key == "n_drop") {
            is >> cfg.n_drop;
        } else if (key == "n_samples") {
            is >> cfg.n_samples;
        } else if (key == "n_separation") {
            is >> cfg.n_separation;
        } else if (key == "init_order") {
            is >> cfg.init_order;
        } else if (key == "restart") {
            is >> cfg.restart;
        } else if (key == "checkpoint_every") {
            is >> cfg.checkpoint_every;
        } else if (key == "walltime_hours") {
            is >> cfg.walltime_hours;
        } else {
            fail("unknown parameter '" + key + "'");
        }
        if (is.fail()) fail("malformed value for '" + key + "'");
    }

    // --- Validation ------------------------------------------------------
    if (!have_elements) fail("missing required 'elements' line");
    if (cfg.N < 1) fail("N must be >= 1");
    if (cfg.T_init <= 0.0 || cfg.T_final <= 0.0)
        fail("T_init and T_final must be > 0");
    if (cfg.n_runs < 1) fail("n_runs must be >= 1");
    if (cfg.n_drop < 0 || cfg.n_samples < 0 || cfg.n_separation < 0)
        fail("MC step counts must be >= 0");
    if (cfg.init_order != 0 && cfg.init_order != 1)
        fail("init_order must be 0 or 1");
    if (cfg.restart != 0 && cfg.restart != 1)
        fail("restart must be 0 or 1");
    return cfg;
}

void print_config(const SimConfig& cfg) {
    std::printf("# ===== alloy-mc configuration =====\n");
    std::printf("# elements      :");
    for (const auto& e : cfg.elements) std::printf(" %s", e.c_str());
    std::printf("  (n_elements = %d)\n", cfg.n_elements());
    std::printf("# lattice       : %s\n", cfg.lattice.c_str());
    std::printf("# energy_unit   : %s\n", cfg.energy_unit.c_str());
    std::printf("# intercept     : %.12g\n", cfg.intercept);
    std::printf("# N             : %d  (sites = %d)\n", cfg.N,
                cfg.N * cfg.N * cfg.N);
    std::printf("# T_init/T_final: %g / %g K\n", cfg.T_init, cfg.T_final);
    std::printf("# n_runs        : %d\n", cfg.n_runs);
    std::printf("# n_drop        : %g\n", cfg.n_drop);
    std::printf("# n_samples     : %g\n", cfg.n_samples);
    std::printf("# n_separation  : %g\n", cfg.n_separation);
    std::printf("# init_order    : %d\n", cfg.init_order);
    std::printf("# restart       : %d\n", cfg.restart);
    std::printf("# ==================================\n");
}

}  // namespace alloy
