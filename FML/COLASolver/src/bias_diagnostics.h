#pragma once
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace bias_diag {

inline void msg(bool ok, const std::string &s) {
    if (FML::ThisTask == 0) std::cout << (ok ? "  [ ok ] " : "  [WARN] ") << s << "\n";
}
inline double gsum(double x) { double r; MPI_Allreduce(&x, &r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return r; }
inline double gmin(double x) { double r; MPI_Allreduce(&x, &r, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); return r; }
inline double gmax(double x) { double r; MPI_Allreduce(&x, &r, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); return r; }

// ---------------------------------------------------------------------------
// 1) Weight integrity
// ---------------------------------------------------------------------------
template <class T>
void check_weights(FML::PARTICLE::MPIParticles<T> &part,
                   FML::GRID::FFTWGrid<NDIM_SPACE> *delta_real,
                   double sigma_sq_grid) {
    T *p = part.get_particles_ptr();
    const size_t np = part.get_npart();

    int nmesh = 0, lx0 = 0, lnx = 0, nz_real = 0;
    double sigma_phys = -1, sigma_raw = -1, grid_mean = 0;
    if (delta_real) {
        nmesh = delta_real->get_nmesh();
        lx0 = delta_real->get_local_x_start();
        lnx = delta_real->get_local_nx();
        nz_real = 2 * (nmesh / 2 + 1);
        const long long ntot = (long long)nmesh * nmesh * nmesh;

        double ss = 0, sm = 0;
        for (int ix = 0; ix < lnx; ix++)
            for (int iy = 0; iy < nmesh; iy++)
                for (int iz = 0; iz < nmesh; iz++) {
                    long long idx = ((long long)ix * nmesh + iy) * nz_real + iz;
                    double v = delta_real->get_real_from_index(idx);
                    ss += v * v; sm += v;
                }
        sigma_phys = std::sqrt(gsum(ss) / ntot);
        grid_mean = gsum(sm) / ntot;

        const long long n_raw = (long long)lnx * nmesh * nz_real;
        const double *rd = delta_real->get_real_grid();
        double sr = 0;
        for (long long i = 0; i < n_raw; i++) sr += rd[i] * rd[i];
        sigma_raw = std::sqrt(gsum(sr) / ntot);
    }

    double n_bad = 0, n_zero1 = 0, n_qout = 0, n_nonlocal = 0, n_mis1 = 0;
    double s1 = 0, s2 = 0, s11 = 0, s22 = 0, s12 = 0, s1111 = 0;
    double mn1 = 1e300, mx1 = -1e300, mn2 = 1e300, mx2 = -1e300;
    const double sg = sigma_sq_grid > 0 ? std::sqrt(sigma_sq_grid) : 1.0;

    for (size_t i = 0; i < np; i++) {
        const double w1 = p[i].bias_weights[1], w2 = p[i].bias_weights[2];
        if (!std::isfinite(w1) || !std::isfinite(w2)) { n_bad++; continue; }
        if (w1 == 0.0) n_zero1++;
        s1 += w1; s2 += w2; s11 += w1 * w1; s22 += w2 * w2; s12 += w1 * w2; s1111 += w1 * w1 * w1 * w1;
        mn1 = std::min(mn1, w1); mx1 = std::max(mx1, w1);
        mn2 = std::min(mn2, w2); mx2 = std::max(mx2, w2);

        for (int d = 0; d < NDIM_SPACE; d++)
            if (p[i].q[d] < 0.0 || p[i].q[d] >= 1.0) { n_qout++; break; }

        if (delta_real) {
            auto cell = [&](double q) { return (int)((((long long)std::floor(q * nmesh + 0.5)) % nmesh + nmesh) % nmesh); };
            int ix = cell(p[i].q[0]), iy = cell(p[i].q[1]), iz = cell(p[i].q[2]);
            if (ix < lx0 || ix >= lx0 + lnx) { n_nonlocal++; continue; }
            long long idx = ((long long)(ix - lx0) * nmesh + iy) * nz_real + iz;
            double ref = delta_real->get_real_from_index(idx);
            if (std::fabs(ref - w1) > 1e-4 * std::fabs(ref) + 1e-6 * sg) n_mis1++;
        }
    }

    const double N = gsum((double)np) - gsum(n_bad);
    n_bad = gsum(n_bad); n_zero1 = gsum(n_zero1); n_qout = gsum(n_qout);
    n_nonlocal = gsum(n_nonlocal); n_mis1 = gsum(n_mis1);
    s1 = gsum(s1); s2 = gsum(s2); s11 = gsum(s11); s22 = gsum(s22); s12 = gsum(s12); s1111 = gsum(s1111);
    mn1 = gmin(mn1); mx1 = gmax(mx1); mn2 = gmin(mn2); mx2 = gmax(mx2);

    const double m1 = s1 / N, m2 = s2 / N;
    const double v1 = s11 / N - m1 * m1, v2 = s22 / N - m2 * m2;
    const double sd1 = std::sqrt(std::max(v1, 0.0)), sd2 = std::sqrt(std::max(v2, 0.0));
    const double corr = (s12 / N - m1 * m2) / (sd1 * sd2 + 1e-300);
    const double kurt = (s1111 / N) / (v1 * v1 + 1e-300);

    if (FML::ThisTask == 0) {
        std::cout << std::scientific << std::setprecision(4)
                  << "\n===== bias weight diagnostics =====\n"
                  << "  N particles (finite)     : " << N << "\n"
                  << "  w1: mean " << m1 << "  rms " << sd1 << "  min " << mn1 << "  max " << mx1 << "\n"
                  << "  w2: mean " << m2 << "  rms " << sd2 << "  min " << mn2 << "  max " << mx2 << "\n"
                  << "  corr(w1,w2) " << corr << " (Gaussian: ~0)   kurtosis(w1) " << kurt << " (Gaussian: ~3)\n";
        if (delta_real)
            std::cout << "  sigma: physical cells " << sigma_phys << " | raw incl. padding " << sigma_raw
                      << " | particle-sampled " << sd1 << " | grid mean " << grid_mean << "\n";
    }
    msg(n_bad == 0, "NaN/Inf weights: " + std::to_string((long long)n_bad));
    msg(n_zero1 < 0.01 * N, "particles with w1 == 0 exactly: " + std::to_string((long long)n_zero1));
    msg(n_qout == 0, "particles with q outside [0,1): " + std::to_string((long long)n_qout));
    msg(std::fabs(m2) < 0.05 * sd2 + 1e-12, "mean(w2) should be ~0 (sampling noise only)");
    msg(std::fabs(kurt - 3.0) < 0.3, "w1 kurtosis close to Gaussian value 3");
    msg(sd2 > 1e-2, "rms(w2) > 1e-2. If not, mass = 1 + w2 is indistinguishable from mass = 1");
    msg(mn1 + 1.0 > 0 && mn2 + 1.0 > 0, "1 + w stays positive (negative masses if this fails)");
    if (delta_real) {
        msg(n_nonlocal == 0, "particles whose q lies outside this task's x-slab: " + std::to_string((long long)n_nonlocal));
        msg(n_mis1 == 0, "stored w1 != delta_L re-sampled at q: " + std::to_string((long long)n_mis1));
        msg(std::fabs(sigma_raw / sigma_phys - 1) < 1e-3, "padding does not change sigma (raw vs physical)");
        msg(std::fabs(sd1 / sigma_phys - 1) < 0.02, "particle-sampled rms(w1) matches grid sigma");
        if (sigma_sq_grid > 0)
            msg(std::fabs(m2 - (s11 / N - sigma_sq_grid)) < 1e-4 * std::max(sd2, sigma_sq_grid),
                "w2 == w1^2 - sigma_sq_used for all particles (consistent sigma_sq)");
    }
}

// ---------------------------------------------------------------------------
// 2) What does get_mass() actually return?
// ---------------------------------------------------------------------------
template <class T>
void check_mass(FML::PARTICLE::MPIParticles<T> &part) {
    T *p = part.get_particles_ptr();
    const size_t np = part.get_npart();
    const double Ntot = part.get_npart_total();
    const int saved = T::active_bias_index;

    if (FML::ThisTask == 0) std::cout << "\n===== get_mass() per field =====\n";
    for (int f = 0; f < 3; f++) {
        T::active_bias_index = f;
        double sm = 0, smm = 0, mn = 1e300, mx = -1e300, nneg = 0;
        for (size_t i = 0; i < np; i++) {
            double m = p[i].get_mass();
            sm += m; smm += m * m; mn = std::min(mn, m); mx = std::max(mx, m); if (m <= 0) nneg++;
        }
        sm = gsum(sm); smm = gsum(smm); mn = gmin(mn); mx = gmax(mx); nneg = gsum(nneg);
        const double mean = sm / Ntot, sd = std::sqrt(std::max(smm / Ntot - mean * mean, 0.0));
        const double neff_frac = (sm * sm / smm) / Ntot;
        if (FML::ThisTask == 0)
            std::cout << std::scientific << std::setprecision(4) << "  field " << f << ": mean " << mean
                      << "  sd " << sd << "  min " << mn << "  max " << mx << "  Neff/N " << neff_frac << "\n";
        if (f > 0) {
            msg(nneg == 0, "field " + std::to_string(f) + ": non-positive masses: " + std::to_string((long long)nneg));
            msg(sd / std::fabs(mean) > 1e-3, "field " + std::to_string(f) + ": relative mass spread > 1e-3");
        }
    }
    T::active_bias_index = saved;
}

// ---------------------------------------------------------------------------
// 3) Quick deposit spectra test
// ---------------------------------------------------------------------------
template <class T, class Method>
void check_spectra(FML::PARTICLE::MPIParticles<T> &part, int nmesh, const Method &method, bool interlacing) {
    const auto [nleft, nright] = FML::INTERPOLATION::get_extra_slices_needed_for_density_assignment(method);
    const int saved = T::active_bias_index;

    std::unique_ptr<FML::GRID::FFTWGrid<NDIM_SPACE>> g[3];
    for (int i = 0; i < 3; i++) {
        T::active_bias_index = i;
        g[i] = std::make_unique<FML::GRID::FFTWGrid<NDIM_SPACE>>(nmesh, nleft, nright);
        FML::INTERPOLATION::particles_to_fourier_grid(part.get_particles_ptr(), part.get_npart(),
                                                      part.get_npart_total(), *g[i], method, interlacing);
        FML::INTERPOLATION::deconvolve_window_function_fourier<NDIM_SPACE>(*g[i], method);
    }
    T::active_bias_index = saved;

    std::vector<double> k, P[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = i; j < 3; j++) {
            FML::CORRELATIONFUNCTIONS::PowerSpectrumBinning<NDIM_SPACE> b(nmesh / 2);
            b.subtract_shotnoise = false;
            if (i == j) FML::CORRELATIONFUNCTIONS::bin_up_power_spectrum(*g[i], b);
            else        FML::CORRELATIONFUNCTIONS::bin_up_cross_power_spectrum(*g[i], *g[j], b);
            P[i][j].resize(b.n);
            for (int kk = 0; kk < b.n; kk++) P[i][j][kk] = b.pofk[kk];
            if (i == 0 && j == 0) { k.resize(b.n); for (int kk = 0; kk < b.n; kk++) k[kk] = b.kbin[kk]; }
        }

    const int nb = (int)k.size();
    double d20 = 0, d21 = 0, d10 = 0, rmin12 = 1, rmin02 = 1;
    for (int kk = 0; kk < nb; kk++) {
        if (P[0][0][kk] <= 0 || P[1][1][kk] <= 0 || P[2][2][kk] <= 0) continue;
        d20 = std::max(d20, std::fabs(P[2][2][kk] / P[0][0][kk] - 1));
        d21 = std::max(d21, std::fabs(P[2][2][kk] / P[1][1][kk] - 1));
        d10 = std::max(d10, std::fabs(P[1][1][kk] / P[0][0][kk] - 1));
        rmin12 = std::min(rmin12, P[1][2][kk] / std::sqrt(P[1][1][kk] * P[2][2][kk]));
        rmin02 = std::min(rmin02, P[0][2][kk] / std::sqrt(P[0][0][kk] * P[2][2][kk]));
    }
    if (FML::ThisTask == 0) {
        std::cout << std::scientific << std::setprecision(3)
                  << "\n===== spectra (no shot-noise subtraction) =====\n"
                  << "       k        P00        P11        P22     P11/P00    P22/P00    P22/P11    r12        r02\n";
        const int step = std::max(1, nb / 8);
        for (int kk = 0; kk < nb; kk += step) {
            if (P[0][0][kk] <= 0) continue;
            std::cout << std::setw(10) << k[kk] << " " << P[0][0][kk] << " " << P[1][1][kk] << " " << P[2][2][kk]
                      << " " << P[1][1][kk] / P[0][0][kk] << " " << P[2][2][kk] / P[0][0][kk]
                      << " " << P[2][2][kk] / P[1][1][kk]
                      << " " << P[1][2][kk] / std::sqrt(P[1][1][kk] * P[2][2][kk])
                      << " " << P[0][2][kk] / std::sqrt(P[0][0][kk] * P[2][2][kk]) << "\n";
        }
    }
    msg(d10 > 1e-3, "field 1 differs from field 0");
    msg(d20 > 1e-3, "field 2 differs from field 0");
    msg(d21 > 1e-3, "field 2 differs from field 1 (if not: index 2 is not being applied)");
}

} // namespace bias_diag