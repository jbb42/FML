"""Reference solutions for the exercises in implement_control_variates.ipynb. Try the exercises first."""
import numpy as np


# Exercise 1
def cv_estimate(x, y, mu_y, beta):
    """The control-variate estimate x - beta (y - mu_y), elementwise."""
    return x - beta * (y - mu_y)


# Exercise 2
def optimal_beta(x, y):
    """beta* = Cov(x, y) / Var(y), from samples (along the first axis)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    dx, dy = x - x.mean(0), y - y.mean(0)
    return (dx * dy).sum(0) / (dy * dy).sum(0)


# Exercise 3
def beta_tanh(k, k0=0.618, dk=0.167):
    """The damping of beta(k) of Kokron et al. 2022 (Eq. 4.7), as in FML: 1 on large scales, 0 on small."""
    return 0.5 * (1.0 - np.tanh((np.asarray(k) - k0) / dk))


# Exercise 4
def binned_mean(P_of_k, K, which, weight, nbins):
    """The mean of a function of |k| over the modes of each bin, weighted like the measurement."""
    vals = P_of_k(np.where(K > 0, K, 1.0)).ravel() * weight
    ok = (which >= 0) & (which < nbins)
    sums = np.bincount(which[ok], weights=vals[ok], minlength=nbins)
    counts = np.bincount(which[ok], weights=weight[ok], minlength=nbins)
    return sums / np.maximum(counts, 1)


# Exercise 5
def variance_reduction(P_sims, P_cvs):
    """Var(P) / Var(P_cv) per k bin, over the realisations (first axis)."""
    return np.var(P_sims, axis=0, ddof=1) / np.var(P_cvs, axis=0, ddof=1)


# Exercise 7
PAIRS = [(i, j) for i in range(5) for j in range(i, 5)]


def load_hybridbias(path):
    """k, number of modes and the 15 spectra [15, nk] of a FML hybridbias*_z*.txt file."""
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1], data[:, 2:2 + len(PAIRS)].T


def control_variate(P, P_field, P_pred, k, k0=0.618, dk=0.167):
    """P_cv = P - beta(k) (P_field - P_pred) for all 15 spectra at once ([15, nk] arrays)."""
    return P - beta_tanh(k, k0, dk) * (P_field - P_pred)


def boost_cv(P_mg, Pf_mg, Ppred_mg, P_gr, Pf_gr, Ppred_gr, k, k0=0.618, dk=0.167):
    """The MG/GR boost with the control variate applied to each model separately (each with its own prediction)."""
    return control_variate(P_mg, Pf_mg, Ppred_mg, k, k0, dk) / control_variate(P_gr, Pf_gr, Ppred_gr, k, k0, dk)
