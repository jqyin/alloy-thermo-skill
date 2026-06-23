// rng.hpp -- single, seedable pseudo-random number generator.
//
// The legacy code mixed three generators (a KISS generator, a Mersenne Twister
// from random.h, and the C library rand()).  This consolidates everything onto
// one std::mt19937_64 instance per MPI rank so that a run is fully reproducible
// from (base_seed, rank).
#ifndef ALLOY_MC_RNG_HPP
#define ALLOY_MC_RNG_HPP

#include <cstdint>
#include <random>

namespace alloy {

class Rng {
public:
    // Each rank derives a distinct stream from the shared base seed so that
    // replicas explore configuration space independently but reproducibly.
    explicit Rng(std::uint64_t base_seed, int rank = 0) {
        std::seed_seq seq{base_seed, static_cast<std::uint64_t>(rank),
                          0x9e3779b97f4a7c15ULL};
        gen_.seed(seq);
    }

    // Uniform double in [0, 1).
    double uniform() { return unit_(gen_); }

    // Uniform integer in [0, n).
    std::size_t below(std::size_t n) {
        std::uniform_int_distribution<std::size_t> d(0, n - 1);
        return d(gen_);
    }

    std::mt19937_64& engine() { return gen_; }

private:
    std::mt19937_64 gen_;
    std::uniform_real_distribution<double> unit_{0.0, 1.0};
};

}  // namespace alloy

#endif  // ALLOY_MC_RNG_HPP
