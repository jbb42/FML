#ifndef LAGRANGIANBIAS_HEADER
#define LAGRANGIANBIAS_HEADER

#include <FML/ComputePowerSpectra/ComputePowerSpectrum.h>
#include <FML/FFTWGrid/FFTWGrid.h>
#include <FML/Global/Global.h>
#include <FML/Interpolation/ParticleGridInterpolation.h>
#include <FML/MPIParticles/MPIParticles.h>

#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace FML {

    /// Basis fields of the hybrid Lagrangian bias expansion. Each particle carries the weights
    ///   w = {1, delta_L, delta_L^2, s^2, nabla^2 delta_L}
    /// evaluated at its Lagrangian position q (all but the first minus their particle mean). Painting the
    /// particles at their Eulerian positions with these weights gives the 5 advected fields, whose 15 auto
    /// and cross power spectra are written to file.
    ///
    /// The particle type T must have double pos[3], double q[3], double bias_weights[n_bias_fields],
    /// static int active_bias_index, and get_mass() returning 1 for index 0 and 1 + bias_weights[index] otherwise.
    namespace LAGRANGIANBIAS {

        template <int N>
        using FFTWGrid = FML::GRID::FFTWGrid<N>;
        template <class T>
        using MPIParticles = FML::PARTICLE::MPIParticles<T>;

        constexpr int n_bias_fields = 5;

        /// Ordered (name, value) lines describing the simulation, written to the info file
        using InfoList = std::vector<std::pair<std::string, std::string>>;

        template <class V>
        std::string info_string(const V & value) {
            std::ostringstream ss;
            ss << std::setprecision(10) << value;
            return ss.str();
        }

        /// Set the weights to {1, 0, 0, 0, 0}: matter only (used when there is no linear field)
        template <class T>
        void reset_bias_weights(MPIParticles<T> & part) {
            for (auto & particle : part)
                for (int i = 0; i < n_bias_fields; i++)
                    particle.bias_weights[i] = (i == 0) ? 1.0 : 0.0;
        }

        //=================================================================================
        /// Compute the weights at the particles' q from delta_L(k) = growth_ratio_of_k(k) * delta_ini(k), where
        /// growth_ratio_of_k = D(k, a_out) / D(k, a_ini) with k in h/Mpc, which is exact also for scale-dependent
        /// growth. Does nothing if delta_ini_fourier is empty (particles read from file).
        //=================================================================================
        template <int NDIM, class T>
        void compute_bias_weights(MPIParticles<T> & part,
                                  const FFTWGrid<NDIM> & delta_ini_fourier,
                                  double boxsize,
                                  const std::function<double(double)> & growth_ratio_of_k) {
            static_assert(NDIM == 3, "The Lagrangian bias weights are only implemented for NDIM = 3");
            const int nmesh = delta_ini_fourier.get_nmesh();
            if (nmesh == 0)
                return;

            // Swapping pos and q (and moving the particles to the MPI task owning their new pos) before and
            // after makes all interpolation below happen at q
            auto swap_pos_and_q = [&]() {
                for (auto & particle : part)
                    for (int idim = 0; idim < NDIM; idim++)
                        std::swap(particle.pos[idim], particle.q[idim]);
                part.communicate_particles();
            };
            swap_pos_and_q();
            reset_bias_weights(part);

            // The growth ratio depends only on |k| = kfund * sqrt(n2) with integer n2 = nx^2 + ny^2 + nz^2
            const double kfund = 2.0 * M_PI / boxsize;
            std::vector<double> growth_of_n2(3 * (nmesh / 2) * (nmesh / 2) + 1);
            for (size_t n2 = 0; n2 < growth_of_n2.size(); n2++)
                growth_of_n2[n2] = growth_ratio_of_k(kfund * std::sqrt(double(n2)));

            // Interpolate kernel(k) * delta_L(k) to every particle and call update(bias_weights, value)
            FFTWGrid<NDIM> grid(nmesh, delta_ini_fourier.get_n_extra_slices_left(), delta_ini_fourier.get_n_extra_slices_right());
            std::vector<double> values;
            auto add_field = [&](auto && kernel, auto && update) {
                const int local_nx = grid.get_local_nx();
                const int local_x_start = grid.get_local_x_start();
                grid.set_grid_status_real(false);
                #ifdef USE_OMP
                #pragma omp parallel for collapse(2)
                #endif
                for (int ix = 0; ix < local_nx; ix++) {
                    for (int iy = 0; iy < nmesh; iy++) {
                        const int nx = ix + local_x_start > nmesh / 2 ? ix + local_x_start - nmesh : ix + local_x_start;
                        const int ny = iy > nmesh / 2 ? iy - nmesh : iy;
                        for (int nz = 0; nz <= nmesh / 2; nz++) {
                            const int n2 = nx * nx + ny * ny + nz * nz;
                            const double k[3] = {nx * kfund, ny * kfund, nz * kfund};
                            const size_t c = ((size_t)ix * nmesh + iy) * (nmesh / 2 + 1) + nz;
                            grid.set_fourier_from_index(c, delta_ini_fourier.get_fourier_from_index(c) *
                                FML::GRID::FloatType(growth_of_n2[n2] * kernel(k, n2 * kfund * kfund)));
                        }
                    }
                }
                grid.fftw_c2r();
                FML::INTERPOLATION::interpolate_grid_to_particle_positions(grid, part.get_particles_ptr(), part.get_npart(), values, "CIC");
                for (size_t p = 0; p < part.get_npart(); p++)
                    update(part[p].bias_weights, values[p]);
            };

            // delta_L and delta_L^2, nabla^2 delta_L, and s^2 = s_ij s_ij with s_ij = (k_i k_j / k^2 - delta_ij / 3) delta_L
            add_field([](const double *, double) { return 1.0; }, [](double * w, double v) { w[1] = v; w[2] = v * v; });
            add_field([](const double *, double k2) { return -k2; }, [](double * w, double v) { w[4] = v; });
            const std::pair<int, int> tidal_components[] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}};
            for (auto [i, j] : tidal_components) {
                const double mult = (i == j) ? 1.0 : 2.0; // s_ij = s_ji
                add_field([i = i, j = j](const double * k, double k2) {
                    return k2 > 0.0 ? k[i] * k[j] / k2 - (i == j ? 1.0 / 3.0 : 0.0) : 0.0;
                }, [mult](double * w, double v) { w[3] += mult * v * v; });
            }

            // Subtract the particle mean of the bias weights (for delta_L^2 this is sigma^2), so the mean of
            // get_mass() = 1 + weight is exactly 1, as the density assignment assumes
            double mean[n_bias_fields] = {0.0};
            for (auto & particle : part)
                for (int i = 1; i < n_bias_fields; i++)
                    mean[i] += particle.bias_weights[i];
            FML::SumArrayOverTasks(mean, n_bias_fields);
            for (auto & particle : part)
                for (int i = 1; i < n_bias_fields; i++)
                    particle.bias_weights[i] -= mean[i] / double(part.get_npart_total());

            swap_pos_and_q();
        }

        //=================================================================================
        /// Paint the 5 weighted fields at the particles' Eulerian positions and write output_folder/pofk_ij.txt
        /// for all 15 pairs ij, plus output_folder/pofk_bias_info.txt with the lines in info and the Poisson shot
        /// noise of each pair. Call compute_bias_weights first.
        //=================================================================================
        template <int NDIM, class T>
        void compute_bias_power_spectra(MPIParticles<T> & part,
                                        double boxsize,
                                        int nmesh,
                                        std::string density_assignment_method,
                                        bool interlacing,
                                        std::string output_folder,
                                        const InfoList & info) {
            constexpr int n = n_bias_fields;
            const auto [nleft, nright] =
                FML::INTERPOLATION::get_extra_slices_needed_for_density_assignment(density_assignment_method);

            // Paint field i with mass 1 + weight i, and subtract the unit-mass field 0 to get the weighted field.
            // All grids are kept so the 15 spectra need no further FFTs
            std::vector<FFTWGrid<NDIM>> grids;
            grids.reserve(n);
            for (int i = 0; i < n; i++) {
                grids.emplace_back(nmesh, nleft, nright);
                T::active_bias_index = i;
                FML::INTERPOLATION::particles_to_fourier_grid(part.get_particles_ptr(), part.get_npart(),
                    part.get_npart_total(), grids[i], density_assignment_method, interlacing);
                FML::INTERPOLATION::deconvolve_window_function_fourier<NDIM>(grids[i], density_assignment_method);
                if (i > 0)
                    for (auto && c : grids[i].get_fourier_range())
                        grids[i].set_fourier_from_index(c, grids[i].get_fourier_from_index(c) - grids[0].get_fourier_from_index(c));
            }
            T::active_bias_index = 0;

            // Poisson shot noise P_ij^SN = V/N <w_i w_j>
            std::vector<double> shotnoise(n * n, 0.0);
            for (auto & particle : part)
                for (int i = 0; i < n; i++)
                    for (int j = i; j < n; j++)
                        shotnoise[i * n + j] += particle.bias_weights[i] * particle.bias_weights[j];
            FML::SumArrayOverTasks(shotnoise.data(), int(shotnoise.size()));
            const double npart_total = double(part.get_npart_total());
            for (auto & sn : shotnoise)
                sn *= std::pow(boxsize, NDIM) / npart_total / npart_total;

            for (int i = 0; i < n; i++) {
                for (int j = i; j < n; j++) {
                    FML::CORRELATIONFUNCTIONS::PowerSpectrumBinning<NDIM> pofk(nmesh / 2);
                    pofk.subtract_shotnoise = false;
                    if (i == j)
                        FML::CORRELATIONFUNCTIONS::bin_up_power_spectrum(grids[i], pofk);
                    else
                        FML::CORRELATIONFUNCTIONS::bin_up_cross_power_spectrum(grids[i], grids[j], pofk);
                    pofk.scale(boxsize);

                    if (FML::ThisTask == 0) {
                        std::ofstream fp(output_folder + "/pofk_" + std::to_string(i) + std::to_string(j) + ".txt");
                        for (int k = 0; k < pofk.n; k++)
                            fp << pofk.kbin[k] << " " << pofk.pofk[k] << "\n";
                    }
                }
            }

            if (FML::ThisTask == 0) {
                std::ofstream fp(output_folder + "/pofk_bias_info.txt");
                fp << std::setprecision(10) << std::left;
                fp << "# Reference for the files pofk_ij.txt in this folder\n";
                fp << "# Columns: k (h/Mpc)   P_ij(k) (Mpc/h)^3\n";
                fp << "# Fields: 0 = matter, 1 = delta_L, 2 = delta_L^2 - <delta_L^2>, 3 = s^2 - <s^2>,\n";
                fp << "#         4 = nabla^2 delta_L (in (Mpc/h)^2). Weights are linear fields at q (mean subtracted)\n";
                fp << "#         advected with the particles, built from the linear field at this redshift (D(k,z) per mode)\n";
                fp << "# Window function deconvolved, shot noise NOT subtracted (Poisson estimate below)\n";
                for (const auto & [name, value] : info)
                    fp << std::setw(31) << name << value << "\n";
                fp << "# Poisson shot noise P_ij^SN = V/N <u_i u_j> in (Mpc/h)^3\n";
                for (int i = 0; i < n; i++)
                    for (int j = i; j < n; j++)
                        fp << std::setw(31) << "shotnoise_" + std::to_string(i) + std::to_string(j) << shotnoise[i * n + j] << "\n";
            }
        }

    } // namespace LAGRANGIANBIAS
} // namespace FML

#endif
