#include "parallel_tempering.hpp"

#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>

#include "units.hpp"

namespace alloy {

namespace {
constexpr int kTag = 123;
}

ParallelTempering::ParallelTempering(const SimConfig& cfg, const Lattice& lat,
                                     std::vector<double> composition, int rank,
                                     int nprocs, std::uint64_t base_seed)
    : cfg_(cfg),
      lat_(lat),
      composition_(std::move(composition)),
      rank_(rank),
      nprocs_(nprocs),
      rng_(base_seed, rank),
      model_(lat, cfg.intercept) {
    e_scale_ = energy_unit_to_meV(cfg_.energy_unit);
    build_temperature_ladder();
    pT_ = T_[rank_];
}

void ParallelTempering::build_temperature_ladder() {
    T_.assign(nprocs_, cfg_.T_init);
    if (nprocs_ == 1) {
        T_[0] = cfg_.T_init;
        return;
    }
    // Geometric (log-uniform) spacing keeps swap acceptance roughly constant.
    double l0 = std::log(cfg_.T_init);
    double l1 = std::log(cfg_.T_final);
    for (int i = 0; i < nprocs_; ++i)
        T_[i] = std::exp(l0 + (l1 - l0) * i / (nprocs_ - 1));
}

void ParallelTempering::sweep() {
    double beta_eff = e_scale_ / (kB_meV_per_K * pT_);
    for (long c = 0; c < lat_.n_sites(); ++c) {
        AlloyModel::Move m = model_.attempt_swap(rng_, beta_eff, pE_);
        if (m == AlloyModel::Move::NoOp) continue;
        ++att_move_;
        if (m == AlloyModel::Move::Accepted) ++acc_move_;
    }
}

void ParallelTempering::replica_exchange(bool even_round) {
    const int low_parity = even_round ? 0 : 1;
    const bool i_am_low = (rank_ % 2 == low_parity) && (rank_ + 1 < nprocs_);
    const bool i_am_high =
        (rank_ % 2 != low_parity) && (rank_ - 1 >= 0) &&
        ((rank_ - 1) % 2 == low_parity);

    MPI_Status st;
    if (i_am_low) {
        int hi = rank_ + 1;
        double Ei = pE_, Ej = 0.0;
        MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, hi, kTag, &Ej, 1, MPI_DOUBLE, hi, kTag,
                     MPI_COMM_WORLD, &st);
        // Metropolis criterion for exchanging configurations between the
        // (low T = rank_) and (high T = hi) replicas.
        double delta = (Ei - Ej) * e_scale_ *
                       (1.0 / (kB_meV_per_K * T_[hi]) -
                        1.0 / (kB_meV_per_K * T_[rank_]));
        int flag = (delta <= 0.0 || rng_.uniform() < std::exp(-delta)) ? 1 : 0;
        MPI_Send(&flag, 1, MPI_INT, hi, kTag, MPI_COMM_WORLD);
        ++att_swap_;
        if (flag) {
            std::vector<std::int16_t> recv(lat_.n_sites());
            MPI_Sendrecv(model_.atoms().data(), static_cast<int>(lat_.n_sites()),
                         MPI_INT16_T, hi, kTag, recv.data(),
                         static_cast<int>(lat_.n_sites()), MPI_INT16_T, hi, kTag,
                         MPI_COMM_WORLD, &st);
            model_.set_atoms(recv);
            pE_ = Ej;
            ++acc_swap_;
        }
    } else if (i_am_high) {
        int lo = rank_ - 1;
        double Ei = pE_, Ej = 0.0;
        MPI_Sendrecv(&Ei, 1, MPI_DOUBLE, lo, kTag, &Ej, 1, MPI_DOUBLE, lo, kTag,
                     MPI_COMM_WORLD, &st);
        int flag = 0;
        MPI_Recv(&flag, 1, MPI_INT, lo, kTag, MPI_COMM_WORLD, &st);
        if (flag) {
            std::vector<std::int16_t> recv(lat_.n_sites());
            MPI_Sendrecv(model_.atoms().data(), static_cast<int>(lat_.n_sites()),
                         MPI_INT16_T, lo, kTag, recv.data(),
                         static_cast<int>(lat_.n_sites()), MPI_INT16_T, lo, kTag,
                         MPI_COMM_WORLD, &st);
            model_.set_atoms(recv);
            pE_ = Ej;  // Ej is the low replica's energy
        }
    }
}

void ParallelTempering::run() {
    t_start_ = MPI_Wtime();
    for (int irun = 0; irun < cfg_.n_runs; ++irun) run_once(irun);
}

void ParallelTempering::run_once(int irun) {
    const int npair = lat_.n_elements() * (lat_.n_elements() - 1) / 2;
    const int nM = npair + 1;  // pair alphas + mean|alpha|

    double avgE = 0.0, avgE2 = 0.0;
    std::vector<double> avgM(nM, 0.0), avgM2(nM, 0.0), avgM4(nM, 0.0);

    long start_sample = 0;
    bool resumed = false;
    if (cfg_.restart == 1)
        resumed = read_checkpoint(irun, start_sample, avgE, avgE2, avgM, avgM2,
                                  avgM4);

    bool even_round = true;
    if (!resumed) {
        // Fresh start: build a configuration and equilibrate.
        model_.init_configuration(composition_, rng_, cfg_.init_order == 1);
        pE_ = model_.energy();
        for (long mcs = 0; mcs < static_cast<long>(cfg_.n_drop); ++mcs) {
            sweep();
            if (mcs % 2 == 0) {
                replica_exchange(even_round);
                even_round = !even_round;
            }
        }
    }

    // --- Measurement phase ---------------------------------------------
    const long n_samples = static_cast<long>(cfg_.n_samples);
    const long sep = static_cast<long>(cfg_.n_separation);
    for (long mcs = start_sample; mcs < n_samples; ++mcs) {
        for (long i = 0; i < sep; ++i) {
            sweep();
            if ((mcs * sep + i) % 2 == 0) {
                replica_exchange(even_round);
                even_round = !even_round;
            }
        }

        avgE += pE_;
        avgE2 += pE_ * pE_;
        std::vector<double> M = model_.order_parameters(/*order_shells=*/1);
        for (int p = 0; p < nM; ++p) {
            avgM[p] += M[p];
            avgM2[p] += M[p] * M[p];
            avgM4[p] += M[p] * M[p] * M[p] * M[p];
        }

        // Periodic wall-time check / checkpoint.
        if (cfg_.checkpoint_every > 0 && mcs % cfg_.checkpoint_every == 0) {
            int stop = 0;
            if (rank_ == 0) {
                double elapsed_h = (MPI_Wtime() - t_start_) / 3600.0;
                if (elapsed_h > 0.9 * cfg_.walltime_hours) stop = 1;
            }
            MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (stop) {
                write_checkpoint(irun, mcs, avgE, avgE2, avgM, avgM2, avgM4);
                if (rank_ == 0)
                    std::fprintf(stderr,
                                 "# wall-time budget reached; checkpoint written "
                                 "at sample %ld of run %d\n",
                                 mcs, irun);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Finalize();
                std::exit(0);
            }
        }
    }

    // --- Reduce to per-replica observables ------------------------------
    double inv = 1.0 / static_cast<double>(n_samples);
    avgE *= inv;
    avgE2 *= inv;
    for (int p = 0; p < nM; ++p) {
        avgM[p] *= inv;
        avgM2[p] *= inv;
        avgM4[p] *= inv;
    }

    const long n3 = lat_.n_sites();
    double kT = kB_meV_per_K * pT_;
    // Specific heat per site (in units of kB): Var(E)/(kB T)^2 / N^3.
    double C = (avgE2 - avgE * avgE) * e_scale_ * e_scale_ / (kT * kT);
    double C_per_site = C / n3;
    double E_per_site = avgE / n3;
    // Susceptibility and Binder cumulant of the mean |alpha| order parameter.
    double chi = (avgM2[npair] - avgM[npair] * avgM[npair]) * n3 / kT;
    double binder =
        avgM2[npair] > 0
            ? 1.0 - avgM4[npair] / (avgM2[npair] * avgM2[npair]) / 3.0
            : 0.0;
    double swap_acc = att_swap_ > 0 ? static_cast<double>(acc_swap_) / att_swap_ : 0.0;
    double move_acc = att_move_ > 0 ? static_cast<double>(acc_move_) / att_move_ : 0.0;

    // Row layout gathered on rank 0 (T is added from the ladder there):
    //   E_per_site, C_per_site, alpha_0..alpha_{npair-1}, alpha_mean, chi,
    //   binder, swap_acc, move_acc
    const int row = 2 + nM + 4;
    std::vector<double> local(row);
    int c = 0;
    local[c++] = E_per_site;
    local[c++] = C_per_site;
    for (int p = 0; p < nM; ++p) local[c++] = avgM[p];
    local[c++] = chi;
    local[c++] = binder;
    local[c++] = swap_acc;
    local[c++] = move_acc;

    std::vector<double> gathered;
    if (rank_ == 0) gathered.resize(static_cast<long>(row) * nprocs_);
    MPI_Gather(local.data(), row, MPI_DOUBLE,
               rank_ == 0 ? gathered.data() : nullptr, row, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);

    if (rank_ == 0) write_results(irun, gathered);
}

void ParallelTempering::write_results(int irun,
                                      const std::vector<double>& gathered) const {
    const int ne = lat_.n_elements();
    const int npair = ne * (ne - 1) / 2;
    const int nM = npair + 1;
    const int row = 2 + nM + 4;

    char fname[256];
    std::snprintf(fname, sizeof(fname), "thermo_run%d.csv", irun);
    std::ofstream f(fname);
    if (!f) throw std::runtime_error("cannot write results file");

    // Header with descriptive, machine-parseable column names.
    f << "T,E_per_site,C_per_site";
    for (int a = 0; a < ne - 1; ++a)
        for (int b = a + 1; b < ne; ++b)
            f << ",alpha_" << cfg_.elements[a] << "-" << cfg_.elements[b];
    f << ",alpha_mean,chi,binder,swap_accept,move_accept\n";

    f.setf(std::ios::scientific);
    f.precision(8);
    for (int n = 0; n < nprocs_; ++n) {
        const double* r = &gathered[static_cast<long>(n) * row];
        f << T_[n];
        for (int k = 0; k < row; ++k) f << "," << r[k];
        f << "\n";
    }
}

// --- Checkpoint I/O (binary, one file per rank) ------------------------
bool ParallelTempering::write_checkpoint(
    int irun, long sample_index, double avgE, double avgE2,
    const std::vector<double>& avgM, const std::vector<double>& avgM2,
    const std::vector<double>& avgM4) const {
    char fname[256];
    std::snprintf(fname, sizeof(fname), "checkpoint_run%d_rank%d.bin", irun,
                  rank_);
    std::ofstream f(fname, std::ios::binary);
    if (!f) return false;
    long n = lat_.n_sites();
    int nM = static_cast<int>(avgM.size());
    f.write(reinterpret_cast<const char*>(&sample_index), sizeof(long));
    f.write(reinterpret_cast<const char*>(&pE_), sizeof(double));
    f.write(reinterpret_cast<const char*>(&avgE), sizeof(double));
    f.write(reinterpret_cast<const char*>(&avgE2), sizeof(double));
    f.write(reinterpret_cast<const char*>(&nM), sizeof(int));
    f.write(reinterpret_cast<const char*>(avgM.data()), sizeof(double) * nM);
    f.write(reinterpret_cast<const char*>(avgM2.data()), sizeof(double) * nM);
    f.write(reinterpret_cast<const char*>(avgM4.data()), sizeof(double) * nM);
    f.write(reinterpret_cast<const char*>(model_.atoms().data()),
            sizeof(std::int16_t) * n);
    return true;
}

bool ParallelTempering::read_checkpoint(int irun, long& sample_index,
                                        double& avgE, double& avgE2,
                                        std::vector<double>& avgM,
                                        std::vector<double>& avgM2,
                                        std::vector<double>& avgM4) {
    char fname[256];
    std::snprintf(fname, sizeof(fname), "checkpoint_run%d_rank%d.bin", irun,
                  rank_);
    std::ifstream f(fname, std::ios::binary);
    if (!f) return false;
    long n = lat_.n_sites();
    int nM = 0;
    f.read(reinterpret_cast<char*>(&sample_index), sizeof(long));
    f.read(reinterpret_cast<char*>(&pE_), sizeof(double));
    f.read(reinterpret_cast<char*>(&avgE), sizeof(double));
    f.read(reinterpret_cast<char*>(&avgE2), sizeof(double));
    f.read(reinterpret_cast<char*>(&nM), sizeof(int));
    avgM.resize(nM);
    avgM2.resize(nM);
    avgM4.resize(nM);
    f.read(reinterpret_cast<char*>(avgM.data()), sizeof(double) * nM);
    f.read(reinterpret_cast<char*>(avgM2.data()), sizeof(double) * nM);
    f.read(reinterpret_cast<char*>(avgM4.data()), sizeof(double) * nM);
    std::vector<std::int16_t> atoms(n);
    f.read(reinterpret_cast<char*>(atoms.data()), sizeof(std::int16_t) * n);
    if (!f) return false;
    model_.set_atoms(atoms);
    pE_ = model_.energy();
    return true;
}

}  // namespace alloy
