"""Analysis of the euclid25 tests (run_all.py). Usage: python3 analyze.py A|B|C|D"""
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MUC = os.path.join(os.path.dirname(HERE), "mu_convergence")
PAIRS = [(i, j) for i in range(5) for j in range(i, 5)]
N = ["1", "d", "d2", "s2", "n"]
SHOW = [(0, 0), (1, 1), (2, 2), (3, 3), (0, 2), (1, 3), (0, 4), (4, 4)]


def ld(folder, kind="", z=0.0):
    f = glob.glob(os.path.join(folder, f"snapshot_*_z{z:.3f}", f"hybridbias{kind}_z{z:.3f}.txt"))
    if not f:
        return None, None
    d = np.loadtxt(f[0])
    return d[:, 0], d[:, 2:17].T


def table(title, k, f, bins, fmt="{:+8.2f}"):
    print(title)
    print("  k (h/Mpc)  " + "".join(f"{N[i] + N[j]:>8s}" for i, j in SHOW))
    for lo, hi in bins:
        m = (k > lo) & (k <= hi)
        print(f"  {lo:.2f}-{hi:.2f}  " + "".join(fmt.format(f(PAIRS.index(p), m)) for p in SHOW))


def log_info(folder):
    log = open(os.path.join(folder, "log.txt")).read()
    t = lambda w: float((re.findall(r"Total elapsed time for \[\s*" + re.escape(w) + r"\]:\s+([\d.]+)", log) or [0])[-1])
    errs = sorted(set(l.strip()[:110] for l in log.splitlines() if "ERROR" in l or "Abort" in l))
    peak = (re.findall(r"Peak RSS:\s+([\d.]+) MB \(max", log) or ["?"])[-1]
    mu = re.findall(r"mu\(a = [^\n]*", log)[:2]
    return t, errs, peak, mu


def A():
    bins = [(0.02, 0.1), (0.1, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 1.6)]
    for tag in ["A0.3_n3_planck", "A1.0_n1_w0wa"]:
        k, pred = ld(os.path.join(HERE, "A", f"{tag}_pred", "r"), "_zeldovich")
        runs = sorted(glob.glob(os.path.join(HERE, "A", f"{tag}_s*", "r")))
        F = np.array([ld(r, "_zeldovichfield")[1] for r in runs])
        mean, err = F.mean(0), F.std(0, ddof=1) / np.sqrt(len(F))
        beta = 0.5 * (1 - np.tanh((k - 0.618) / 0.167))
        print(f"\n=== {tag}: {len(F)} Zel'dovich fields against the prediction")
        table(" bias the control variate would add, beta (mean - prediction) / prediction (%)", k,
              lambda p, m: 100 * np.mean(beta[m] * (mean[p][m] - pred[p][m]) / pred[p][m]), bins)
        table(" (mean - prediction) in units of the error of the mean (sigma)", k,
              lambda p, m: np.mean((mean[p][m] - pred[p][m]) / err[p][m]), bins)


def B():
    for folder in sorted(glob.glob(os.path.join(HERE, "B", "*", "r"))):
        name = folder.split(os.sep)[-2]
        done = bool(glob.glob(os.path.join(folder, "snapshot_*_z0.000", "hybridbias_z0.000.txt")))
        t, errs, peak, mu = log_info(folder)
        print(f"{name:28s} finished: {done} | {t('The whole simulation'):.0f} s | peak/task {peak} MB | "
              f"{'; '.join(m.strip() for m in mu)}" + (f" | ERRORS: {errs[:2]}" if errs else ""))


def C():
    for nm in (256, 384, 512):
        folder = os.path.join(HERE, "C", f"nmesh{nm}", "r")
        snap = glob.glob(os.path.join(folder, "snapshot_*_z0.000"))
        files = {os.path.basename(f): os.path.getmtime(f) for f in glob.glob(os.path.join(snap[0], "*"))} if snap else {}
        t, errs, peak, _ = log_info(folder)
        zf, zp = files.get("hybridbias_zeldovichfield_z0.000.txt"), files.get("hybridbias_zeldovich_z0.000.txt")
        print(f"hybridbias_nmesh {nm} (k_Ny {np.pi * nm / 1024:.2f} h/Mpc): prediction "
              f"{(zp - zf) / 60 if zp and zf else float('nan'):.1f} min after the field | whole run "
              f"{t('The whole simulation'):.0f} s | peak/task {peak} MB" + (f" | ERRORS {errs[:2]}" if errs else ""))


def D():
    bins = [(0.01, 0.1), (0.1, 0.3), (0.3, 0.6), (0.6, 1.0)]
    fol = lambda m, c: os.path.join(MUC, f"{'gr' if m == 'GR' else 'mu_' + m}_{c}_L1024_f2048_n1024_s30_normal")

    def P(m, c):
        k, raw = ld(fol(m, c))
        _, zf = ld(fol(m, c), "_zeldovichfield")
        _, zp = ld(fol(m, c), "_zeldovich")
        return k, raw - 0.5 * (1 - np.tanh((k - 0.618) / 0.167)) * (zf - zp)

    k, gf = P("GR", "fiducial")
    _, mf = P("A0.3_n3", "fiducial")
    for c in sorted(os.path.basename(f).split("_")[1] for f in glob.glob(os.path.join(MUC, "gr_ns*_normal"))):
        _, g = P("GR", c)
        _, m = P("A0.3_n3", c)
        table(f"\n{c}: boost / fiducial boost - 1 (%), A0.3_n3, z = 0, with control variate", k,
              lambda p, mm: 100 * np.mean((m[p][mm] / g[p][mm]) / (mf[p][mm] / gf[p][mm]) - 1), bins)
        table(f"{c}: GR spectra / fiducial - 1 (%)", k, lambda p, mm: 100 * np.mean(g[p][mm] / gf[p][mm] - 1), bins)


if __name__ == "__main__":
    for part in sys.argv[1:]:
        globals()[part]()
