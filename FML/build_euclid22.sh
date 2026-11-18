#!/bin/bash
# Native (-march=native = znver3) builds of SimulationCode on euclid22, in double and single precision.
# Run on euclid22 from anywhere: bash /mn/stornext/u3/jonasbbe/pc/Dokumenter/FML/FML/build_euclid22.sh
# Executables: FML/build_euclid22_double/nbody and FML/build_euclid22_float/nbody (ignored by git)
cd /mn/stornext/u3/jonasbbe/pc/Dokumenter/FML/FML
H=/mn/stornext/u3/hansw/local
COMMON="-DCMAKE_BUILD_TYPE=Release -DFML_BUILD_TESTS=OFF -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_PREFIX_PATH=$H;$HOME/local"
# Double precision: FFTW (with threads) from ~/local as the Makefile build, GSL and Lua from hansw
cmake -S SimulationCode -B build_euclid22_double $COMMON -DFFTW_ROOT=$HOME/local > build_euclid22_double.cmake.log 2>&1 \
    && cmake --build build_euclid22_double -j 32 > build_euclid22_double.build.log 2>&1
echo "double: exit $?"
# Single precision: the float FFTW of hansw (no threads library there); since FMLDevelop of 2026-11-18 the perturbation
# theory (FFTLog) also needs the double FFTW, linked statically from hansw
cmake -S SimulationCode -B build_euclid22_float $COMMON -DFFTW_ROOT=$H -DFML_FFTW_PRECISION=float \
    -DFML_USE_FFTW_THREADS=OFF -DCMAKE_CXX_STANDARD_LIBRARIES=$H/lib/libfftw3.a > build_euclid22_float.cmake.log 2>&1 \
    && cmake --build build_euclid22_float -j 32 > build_euclid22_float.build.log 2>&1
echo "float: exit $?"
