#!/usr/bin/env python3
"""Omega_m = Omega_cdm + Omega_b = 0.2, 0.3 and 0.4 (Omega_b = 0.049, h and A_s as in parameterfile.lua), two seeds.
Run as python3 scripts/scan_omega_m.py [--dry-run] [--plot-only]; copy this file to make a new scan.
"""
from scan import scan

scan("omega_m", vary={"cosmology_OmegaCDM": [0.151, 0.251, 0.351]}, seeds=[1001, 1002])
