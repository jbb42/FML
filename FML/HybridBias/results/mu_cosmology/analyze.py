"""Which cosmological parameters change the mu - GR change of the spectra, and by how much compared with the noise?

D_ij = (P_mu - P_GR) / sqrt(P_ii^GR P_jj^GR) in each cosmology (GR and mu from the same seed, z = 0 or 1, no control
variate: the same phase in all runs, so the noise largely cancels between cosmologies). For each model and cosmology:
the change of D relative to the fiducial cosmology, |D_c - D_fid| / max_k |D_fid| (in % of the size of the effect), in
k bins, for the main spectra; and the realisation noise from the seeds 1001, 1002, 1003 of the fiducial cosmology.
    python3 analyze.py [z]
"""
import glob
import os
import sys

import numpy as np

import run_all

HERE = os.path.dirname(os.path.abspath(__file__))
PAIRS = [(i, j) for i in range(5) for j in range(i, 5)]
N = ["1", "d", "d2", "s2", "n"]
SHOW = [(0, 0), (1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3)]  # nabla^2 left out (see stitch.py)
BINS = [(0.01, 0.1), (0.1, 0.3), (0.3, 0.6), (0.6, 1.0)]
Z = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0


def load(label, model):
    f = glob.glob(os.path.join(HERE, label, model, f"snapshot_*_z{Z:.3f}", f"hybridbias_z{Z:.3f}.txt"))
    if not f:
        return None, None
    d = np.loadtxt(f[0])
    return d[:, 0], d[:, 2:17].T


def D(label, model):
    k, g = load(label, "GR")
    _, m = load(label, model)
    if g is None or m is None:
        return None, None
    diag = [np.abs(g[PAIRS.index((i, i))]) for i in range(5)]
    return k, (m - g) / np.array([np.sqrt(diag[i] * diag[j]) for i, j in PAIRS])


def change(Dc, Dref, k):
    """Per shown pair and k bin: mean |Dc - Dref| in the bin, in % of max |Dref| over 0.01 < k < 1."""
    out = np.zeros((len(SHOW), len(BINS)))
    for a, p in enumerate(SHOW):
        i = PAIRS.index(p)
        scale = np.max(np.abs(Dref[i][(k > 0.01) & (k < 1.0)]))
        for b, (lo, hi) in enumerate(BINS):
            m = (k > lo) & (k <= hi)
            out[a, b] = 100 * np.mean(np.abs(Dc[i][m] - Dref[i][m])) / scale
    return out


if __name__ == "__main__":
    labels = list(run_all.COSMOLOGIES) + ["Ob_lo", "Ob_hi"]
    print(f"z = {Z}. Entries: |D_cosmology - D_fiducial| as % of max|D_fiducial|, worst of the spectra "
          f"{', '.join(N[i] + N[j] for i, j in SHOW)}, in k bins {BINS}")
    summary = {}
    for model in [m for m in run_all.MODELS if m != "GR"]:
        k, Df = D("fiducial", model)
        if Df is None:
            continue
        print(f"\n=== {model}")
        noise = [change(D(s, model)[1], Df, k) for s in ("fiducial_s1002", "fiducial_s1003") if D(s, model)[1] is not None]
        if noise:
            n = np.max(np.mean(noise, axis=0), axis=0) / np.sqrt(2)  # One seed vs another: sqrt(2) x the noise of one
            print(f"  {'noise (seeds)':14s}" + "".join(f"{v:8.2f}" for v in n))
        for label in labels[1:]:
            _, Dc = D(label, model)
            if Dc is None:
                continue
            c = change(Dc, Df, k)
            worst = np.max(c, axis=0)
            summary.setdefault(label, []).append(worst)
            which = [N[SHOW[a][0]] + N[SHOW[a][1]] for a in np.argmax(c, axis=0)]
            print(f"  {label:14s}" + "".join(f"{v:8.2f}" for v in worst) + "   (" + ",".join(which) + ")")
    print("\nWorst over all models (%):")
    for label, rows in summary.items():
        print(f"  {label:14s}" + "".join(f"{v:8.2f}" for v in np.max(rows, axis=0)) + f"   ({len(rows)} models)")
