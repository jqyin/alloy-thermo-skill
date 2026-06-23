// parallel_tempering.hpp -- replica-exchange Monte Carlo driver.
//
// One MPI rank owns one replica at a fixed temperature drawn from a geometric
// ladder spanning [T_init, T_final].  Within a replica we run composition-
// conserving Kawasaki sweeps; between replicas we periodically attempt
// configuration exchanges between adjacent temperatures (parallel tempering),
// which dramatically improves sampling across ordering transitions.
#ifndef ALLOY_MC_PARALLEL_TEMPERING_HPP
#define ALLOY_MC_PARALLEL_TEMPERING_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "alloy_model.hpp"
#include "config.hpp"
#include "lattice.hpp"
#include "rng.hpp"

namespace alloy {

class ParallelTempering {
public:
    ParallelTempering(const SimConfig& cfg, const Lattice& lat,
                       std::vector<double> composition, int rank, int nprocs,
                       std::uint64_t base_seed);

    // Run all independent runs; writes thermo_run<i>.csv per run on rank 0.
    void run();

private:
    void build_temperature_ladder();
    void sweep();                       // one MC sweep = n_sites move attempts
    void replica_exchange(bool even_round);
    void run_once(int irun);
    void write_results(int irun, const std::vector<double>& gathered) const;

    bool write_checkpoint(int irun, long sample_index, double avgE, double avgE2,
                          const std::vector<double>& avgM,
                          const std::vector<double>& avgM2,
                          const std::vector<double>& avgM4) const;
    bool read_checkpoint(int irun, long& sample_index, double& avgE,
                         double& avgE2, std::vector<double>& avgM,
                         std::vector<double>& avgM2, std::vector<double>& avgM4);

    const SimConfig& cfg_;
    const Lattice& lat_;
    std::vector<double> composition_;
    int rank_;
    int nprocs_;

    Rng rng_;
    AlloyModel model_;
    std::vector<double> T_;             // full ladder (known to every rank)
    double pT_ = 0.0;                   // this replica's temperature (K)
    double pE_ = 0.0;                   // this replica's current energy (native)
    double e_scale_ = 1.0;             // native energy unit -> meV

    long long att_swap_ = 0, acc_swap_ = 0;   // replica-exchange stats
    long long att_move_ = 0, acc_move_ = 0;   // single-site move stats
    double t_start_ = 0.0;             // MPI_Wtime at run start
};

}  // namespace alloy

#endif  // ALLOY_MC_PARALLEL_TEMPERING_HPP
