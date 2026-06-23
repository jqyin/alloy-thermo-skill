// main.cpp -- entry point for the alloy chemical-ordering Monte Carlo engine.
//
// Usage:
//   mpirun -n <nT> alloy_mc <control.input> [base_seed]
//
// The number of MPI ranks <nT> is the number of temperature replicas in the
// parallel-tempering ladder.  Input files expected in the working directory:
//   <control.input>   simulation parameters (see config.hpp / file_formats.md)
//   composition.input  target species fractions (one per element)
//   coupling.input     per-neighbor distances, J_ij matrices and offset vectors
#include <mpi.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "config.hpp"
#include "lattice.hpp"
#include "parallel_tempering.hpp"

namespace {

std::vector<double> read_composition(const std::string& path, int ne) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open '" + path + "'");
    std::vector<double> c;
    double x;
    while (f >> x) c.push_back(x);
    if (static_cast<int>(c.size()) != ne)
        throw std::runtime_error("composition.input has " +
                                 std::to_string(c.size()) + " values but " +
                                 std::to_string(ne) + " elements were declared");
    double sum = 0.0;
    for (double v : c) sum += v;
    if (sum < 0.99 || sum > 1.01)
        std::fprintf(stderr,
                     "# warning: composition sums to %g (expected ~1.0)\n", sum);
    return c;
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    try {
        if (argc < 2)
            throw std::runtime_error(
                "usage: mpirun -n <nT> alloy_mc <control.input> [base_seed]");

        std::string control = argv[1];
        std::uint64_t base_seed = (argc >= 3)
                                      ? std::strtoull(argv[2], nullptr, 10)
                                      : 12345ULL;

        alloy::SimConfig cfg = alloy::load_config(control);
        std::vector<double> comp =
            read_composition("composition.input", cfg.n_elements());
        alloy::Lattice lat =
            alloy::Lattice::load("coupling.input", cfg.n_elements(), cfg.N);

        if (rank == 0) {
            alloy::print_config(cfg);
            lat.print_summary();
            std::fprintf(stderr, "# replicas (MPI ranks) = %d\n", nprocs);
        }

        alloy::ParallelTempering pt(cfg, lat, comp, rank, nprocs, base_seed);
        pt.run();

        if (rank == 0) std::fprintf(stderr, "# done\n");
    } catch (const std::exception& e) {
        std::fprintf(stderr, "[rank %d] error: %s\n", rank, e.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    MPI_Finalize();
    return 0;
}
