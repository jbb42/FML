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

    /// This namespace deals with the basis fields of the hybrid Lagrangian bias expansion: the weights
    /// {delta_L, delta_L^2 - <delta_L^2>, s^2 - <s^2>, nabla^2 delta_L} evaluated at the particles'
    /// Lagrangian positions q, and the auto and cross power spectra of these fields advected with the particles.
    ///
    /// The particle type T must have:
    ///   double pos[NDIM], double q[NDIM]       (Eulerian and Lagrangian position)
    ///   double bias_weights[n_bias_fields]     (bias_weights[0] is unused, the unit weight is implicit)
    ///   static int active_bias_index           (which weight get_mass() uses)
    ///   double get_mass() const                returning 1 for active_bias_index == 0, else 1 + bias_weights[index]
    namespace LAGRANGIANBIAS {

        template <int N>
        using FFTWGrid = FML::GRID::FFTWGrid<N>;
        template <class T>
        using MPIParticles = FML::PARTICLE::MPIParticles<T>;

        /// Matter (unit weight) + the 4 bias fields
        constexpr int n_bias_fields = 5;

        /// Ordered (name, value) lines added to the info file next to the spectra
        using InfoList = std::vector<std::pair<std::string, std::string>>;

        /// Format a value for an InfoList with the same precision as the rest of the info file
        template <class V>
        std::string info_string(const V & value) {
            std::ostringstream ss;
            ss << std::setprecision(10) << value;
            return ss.str();
        }

        //=================================================================================
        /// Set all bias weights to zero (get_mass() then returns 1 for every field). Used when there is no
        /// linear field to build the weights from, e.g. when particles are read from file.
        //=================================================================================
        template <class T>
        void zero_bias_weights(MPIParticles<T> & part) {
            auto * particle_array = part.get_particles_ptr();
            for (size_t p = 0; p < part.get_npart(); p++)
                for (auto & w : particle_array[p].bias_weights)
                    w = 0.0;
        }

        //=================================================================================
        /// Compute the Lagrangian bias weights {delta_L, delta_L^2, s^2, nabla^2 delta_L} at the particles' q.
        ///
        /// The weights are built from the linear field at the output time,
        /// delta_L(k, a_out) = D(k, a_out) / D(k, a_ini) * delta_ini(k), so they are exact also for
        /// scale-dependent growth (f(R), neutrinos). Each weight has its particle mean subtracted.
        ///
        /// @param[in] part The particles. They are temporarily moved to the MPI domain owning their q
        /// @param[in] delta_ini_fourier The initial density field delta(k, a_ini) the particles were generated from
        /// @param[in] boxsize The boxsize in Mpc/h
        /// @param[in] growth_ratio_of_k D(k, a_out) / D(k, a_ini) as function of k in h/Mpc
        /// @param[in] growth_is_scaledependent If false growth_ratio_of_k is only evaluated at k = 0
        //=================================================================================
        template <int NDIM, class T>
        void compute_bias_weights(MPIParticles<T> & part,
                                  const FFTWGrid<NDIM> & delta_ini_fourier,
                                  double boxsize,
                                  const std::function<double(double)> & growth_ratio_of_k,
                                  bool growth_is_scaledependent) {
            static_assert(NDIM == 3, "The Lagrangian bias weights are only implemented for NDIM = 3");
            if (delta_ini_fourier.get_nmesh() == 0)
                return; // No initial field (particles read from file): the weights stay zero

            // 1. Route particles to the MPI domain that owns their Lagrangian 'q' coordinate
            auto * particle_array = part.get_particles_ptr();
            for (size_t p = 0; p < part.get_npart(); p++) {
                for (int idim = 0; idim < NDIM; idim++) std::swap(particle_array[p].pos[idim], particle_array[p].q[idim]);
            }
            part.communicate_particles();

            // Re-fetch pointer/count as particles have migrated across ranks!
            particle_array = part.get_particles_ptr();
            const size_t local_npart = part.get_npart();

            // 2. Linear growth from the IC redshift to a_out, tabulated in n^2 = |k / kfund|^2 if scale-dependent
            const int nmesh = delta_ini_fourier.get_nmesh();
            const int nmesh_half = nmesh / 2 + 1;
            const int local_nx = delta_ini_fourier.get_local_nx();
            const int local_x_start = delta_ini_fourier.get_local_x_start();
            const double kfund = 2.0 * M_PI / boxsize;

            const double growth_ini_to_out = growth_ratio_of_k(0.0);
            std::vector<double> growth_of_n2;
            if (growth_is_scaledependent) {
                growth_of_n2.resize(3 * (nmesh / 2) * (nmesh / 2) + 1);
                #ifdef USE_OMP
                #pragma omp parallel for
                #endif
                for (size_t n2 = 0; n2 < growth_of_n2.size(); n2++) {
                    growth_of_n2[n2] = growth_ratio_of_k(kfund * std::sqrt(double(n2)));
                }
            }

            // 3. Fill the scratch grid with kernel(k) * D(k) * delta_ini(k), transform and interpolate at q
            FFTWGrid<NDIM> scratch(nmesh, delta_ini_fourier.get_n_extra_slices_left(), delta_ini_fourier.get_n_extra_slices_right());
            auto interpolate_filtered_field = [&](auto && kernel, std::vector<double> & values) {
                scratch.set_grid_status_real(false);
                #ifdef USE_OMP
                #pragma omp parallel for collapse(2)
                #endif
                for (int ix = 0; ix < local_nx; ix++) {
                    for (int iy = 0; iy < nmesh; iy++) {
                        const int ix_global = ix + local_x_start;
                        const int nx = ix_global > nmesh / 2 ? ix_global - nmesh : ix_global;
                        const int ny = iy > nmesh / 2 ? iy - nmesh : iy;
                        for (int iz = 0; iz < nmesh_half; iz++) {
                            const int n2 = nx * nx + ny * ny + iz * iz;
                            const double k[3] = {nx * kfund, ny * kfund, iz * kfund};
                            const double k2 = n2 * kfund * kfund;
                            const double growth = growth_of_n2.empty() ? growth_ini_to_out : growth_of_n2[n2];
                            const size_t c = ((size_t)ix * nmesh + iy) * nmesh_half + iz;
                            scratch.set_fourier_from_index(c, delta_ini_fourier.get_fourier_from_index(c) *
                                                                  FML::GRID::FloatType(growth * kernel(k, k2)));
                        }
                    }
                }
                scratch.fftw_c2r();
                values.resize(local_npart);
                FML::INTERPOLATION::interpolate_grid_to_particle_positions(
                    scratch,
                    particle_array,
                    local_npart,
                    values,
                    "CIC"
                );
            };

            // 4. delta_L, its Laplacian and the tidal tensor invariant s^2 = s_ij s_ij
            std::vector<double> delta_vals, nabla2_vals, comp_vals;
            std::vector<double> s2_vals(local_npart, 0.0);
            interpolate_filtered_field([](const double *, double) { return 1.0; }, delta_vals);
            interpolate_filtered_field([](const double *, double k2) { return -k2; }, nabla2_vals);

            const std::vector<std::pair<int, int>> tidal_comps = { {0,0}, {1,1}, {2,2}, {0,1}, {0,2}, {1,2} };
            for (const auto & comp : tidal_comps) {
                const int i = comp.first, j = comp.second;
                interpolate_filtered_field([i, j](const double * k, double k2) {
                    return k2 > 0.0 ? k[i] * k[j] / k2 - (i == j ? 1.0 / 3.0 : 0.0) : 0.0;
                }, comp_vals);

                const double mult = (i == j) ? 1.0 : 2.0;
                #ifdef USE_OMP
                #pragma omp parallel for
                #endif
                for (size_t p = 0; p < local_npart; p++) {
                    s2_vals[p] += mult * comp_vals[p] * comp_vals[p];
                }
            }

            // 5. Compute the particle means <delta>, sigma^2, <s^2> and <nabla^2 delta>. All weights are made
            // mean-free so that the mean mass used to normalize the density assignment is exactly 1
            double sums[4] = {0.0, 0.0, 0.0, 0.0};
            #ifdef USE_OMP
            #pragma omp parallel for reduction(+:sums[:4])
            #endif
            for (size_t p = 0; p < local_npart; p++) {
                sums[0] += delta_vals[p];
                sums[1] += delta_vals[p] * delta_vals[p];
                sums[2] += s2_vals[p];
                sums[3] += nabla2_vals[p];
            }
            FML::SumArrayOverTasks(sums, 4);

            const double total_parts = static_cast<double>(part.get_npart_total());
            const double delta_mean  = sums[0] / total_parts;
            const double sigma_sq    = sums[1] / total_parts;
            const double s2_mean     = sums[2] / total_parts;
            const double nabla2_mean = sums[3] / total_parts;

            // 6. Assign the bias weights
            #ifdef USE_OMP
            #pragma omp parallel for
            #endif
            for (size_t p = 0; p < local_npart; p++) {
                const double delta_val = delta_vals[p];
                particle_array[p].bias_weights[0] = 0.0;
                particle_array[p].bias_weights[1] = delta_val - delta_mean;
                particle_array[p].bias_weights[2] = delta_val * delta_val - sigma_sq;
                particle_array[p].bias_weights[3] = s2_vals[p] - s2_mean;
                particle_array[p].bias_weights[4] = nabla2_vals[p] - nabla2_mean;
            }

            // 7. Swap back to Eulerian coordinates and restore original MPI domains
            for (size_t p = 0; p < local_npart; p++) {
                for (int idim = 0; idim < NDIM; idim++) std::swap(particle_array[p].pos[idim], particle_array[p].q[idim]);
            }
            part.communicate_particles();
        }

        //=================================================================================
        /// Deposit the matter field and the 4 weighted fields at the particles' Eulerian positions and write
        /// all 15 auto and cross power spectra to output_folder/pofk_ij.txt, plus a reference file
        /// output_folder/pofk_bias_info.txt with the extra lines in info and the Poisson shot noise of each pair.
        /// The weights must have been computed first (compute_bias_weights).
        ///
        /// @param[in] part The particles
        /// @param[in] boxsize The boxsize in Mpc/h
        /// @param[in] nmesh Grid size for the density assignment
        /// @param[in] density_assignment_method NGP, CIC, TSC, PCS or PQS
        /// @param[in] interlacing Use interlaced grids for alias reduction
        /// @param[in] output_folder Folder to write the files to
        /// @param[in] info (name, value) lines describing the simulation, written to the info file
        //=================================================================================
        template <int NDIM, class T>
        void compute_bias_power_spectra(MPIParticles<T> & part,
                                        double boxsize,
                                        int nmesh,
                                        std::string density_assignment_method,
                                        bool interlacing,
                                        std::string output_folder,
                                        const InfoList & info) {
            const int active_fields = n_bias_fields;
            const auto [nleft, nright] =
                FML::INTERPOLATION::get_extra_slices_needed_for_density_assignment(density_assignment_method);

            // Subtract two Fourier space grids a - b and store the result in a
            auto subtract_grid = [](FFTWGrid<NDIM> & a, const FFTWGrid<NDIM> & b) {
                const size_t ncell = (size_t)a.get_local_nx() * a.get_nmesh() * (a.get_nmesh() / 2 + 1);
                #ifdef USE_OMP
                #pragma omp parallel for
                #endif
                for (size_t c = 0; c < ncell; c++)
                    a.set_fourier_from_index(c, a.get_fourier_from_index(c) - b.get_fourier_from_index(c));
            };

            // 1. Deposit each field. Keep all grids so every auto/cross spectrum needs no extra FFTs
            std::vector<FFTWGrid<NDIM>> grids;
            grids.reserve(active_fields);
            for (int i = 0; i < active_fields; i++) {
                grids.emplace_back(nmesh, nleft, nright);
                T::active_bias_index = i;

                FML::INTERPOLATION::particles_to_fourier_grid(part.get_particles_ptr(),
                                                              part.get_npart(),
                                                              part.get_npart_total(),
                                                              grids[i],
                                                              density_assignment_method,
                                                              interlacing);
                FML::INTERPOLATION::deconvolve_window_function_fourier<NDIM>(grids[i], density_assignment_method);

                // get_mass() returns 1 + weight for i > 0. Subtract the unit-mass field to isolate the weighted field
                if (i > 0)
                    subtract_grid(grids[i], grids[0]);
            }
            T::active_bias_index = 0;

            // 2. Poisson shot noise of each pair, P_ij^SN = V/N <u_i u_j> with u_0 = 1 and u_i = weight i
            std::vector<double> shotnoise(active_fields * active_fields, 0.0);
            {
                const auto * p = part.get_particles_ptr();
                for (size_t ip = 0; ip < part.get_npart(); ip++) {
                    double u[active_fields];
                    u[0] = 1.0;
                    for (int i = 1; i < active_fields; i++)
                        u[i] = p[ip].bias_weights[i];
                    for (int i = 0; i < active_fields; i++)
                        for (int j = i; j < active_fields; j++)
                            shotnoise[i * active_fields + j] += u[i] * u[j];
                }
                FML::SumArrayOverTasks(shotnoise.data(), int(shotnoise.size()));
                const double volume_per_particle = std::pow(boxsize, NDIM) / double(part.get_npart_total());
                for (auto & sn : shotnoise)
                    sn *= volume_per_particle / double(part.get_npart_total());
            }

            // 3. Compute auto and cross spectra from the deposited grids
            for (int i = 0; i < active_fields; i++) {
                for (int j = i; j < active_fields; j++) {
                    FML::CORRELATIONFUNCTIONS::PowerSpectrumBinning<NDIM> p_ij(nmesh / 2);
                    p_ij.subtract_shotnoise = false;

                    if (i == j)
                        FML::CORRELATIONFUNCTIONS::bin_up_power_spectrum(grids[i], p_ij);
                    else
                        FML::CORRELATIONFUNCTIONS::bin_up_cross_power_spectrum(grids[i], grids[j], p_ij);
                    p_ij.scale(boxsize);

                    if (FML::ThisTask == 0) {
                        std::string filename = output_folder + "/pofk_" + std::to_string(i) + std::to_string(j) + ".txt";
                        std::ofstream fp(filename);
                        for (int k = 0; k < p_ij.n; k++) {
                            fp << p_ij.kbin[k] << " " << p_ij.pofk[k] << "\n";
                        }
                    }
                }
            }

            // 4. Write a reference file describing the pofk_ij.txt files
            if (FML::ThisTask == 0) {
                std::ofstream fp(output_folder + "/pofk_bias_info.txt");
                fp << std::setprecision(10);
                fp << "# Reference for the files pofk_ij.txt in this folder\n";
                fp << "# Columns: k (h/Mpc)   P_ij(k) (Mpc/h)^3\n";
                fp << "# Fields: 0 = matter, 1 = delta_L, 2 = delta_L^2 - <delta_L^2>, 3 = s^2 - <s^2>,\n";
                fp << "#         4 = nabla^2 delta_L (in (Mpc/h)^2). Weights are linear fields at q (mean subtracted)\n";
                fp << "#         advected with the particles, built from the linear field at this redshift (D(k,z) per mode)\n";
                fp << "# Window function deconvolved, shot noise NOT subtracted (Poisson estimate below)\n";
                for (const auto & [name, value] : info)
                    fp << std::left << std::setw(31) << name << value << "\n";
                fp << "# Poisson shot noise P_ij^SN = V/N <u_i u_j> in (Mpc/h)^3\n";
                for (int i = 0; i < active_fields; i++)
                    for (int j = i; j < active_fields; j++)
                        fp << std::left << std::setw(31) << ("shotnoise_" + std::to_string(i) + std::to_string(j))
                           << shotnoise[i * active_fields + j] << "\n";
            }
        }

    } // namespace LAGRANGIANBIAS
} // namespace FML

#endif
