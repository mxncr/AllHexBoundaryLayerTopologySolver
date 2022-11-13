# Paper

    Robust topological construction of all-hexahedral boundary layer meshes
    Reberol, Maxence; Verhetsel, Kilian; Henrotte, François; Bommes, David; Remacle, Jean-François
    ACM Transactions on Mathematical Software, 2022

- [Preprint](https://raw.githubusercontent.com/mxncr/mxncr.github.io/master/pdf/hexbl_2021_v2.pdf)
- [Supplemental](https://www.hextreme.eu/data/hexbl_supplemental.pdf)
- [Slides](https://mxncr.github.io/pdf/slides_hexbl_frames2021.pdf)

# Overview

This library contains the integer solver of the paper "Robust topological
construction of all-hexahedral boundary layer meshes".  Specifically, it solves
the global integer problem (Eq. 3) with the Gecode library.  The setup of the
branch-and-bound search (including the conditions of Proposition 1) is
implemented in `src/solver.cpp`. To call the integer solver from your code,
call `solveAllHexLayerTopology(...)` which is declared in `src/solver.h`.

This code is distributed to supplement the manuscript. This is not a standalone
meshing tool. Domain decomposition of the surface mesh and hexahedral mesh
construction, which are necessary in practice, are not part of this repository.

The complete implementation (used to generate the results in the paper) is
available in the official [Gmsh](https://gmsh.info/) repository. See the hexbl
branch: https://gitlab.onelab.info/gmsh/gmsh/-/tree/hexbl

# Compilation

This library depends on the `Gecode` library (https://www.gecode.org/, MIT
license).  CMake will look for a system-wide installation, which usually can be
installed via the package manager of your distribution (or homebrew on MacOS).

On Linux / MacOS:

    mkdir build 
    cd build/ 
    cmake ..
    make 

Tests are implemented with Catch2 (https://github.com/catchorg/Catch2, BSL.1-0 License),
whose single header is included in `tests/`.

# Executables

Basic tests are available in `tests/basic_tests.cpp` and can be run with the
executable `basic_tests`.

A list of disk triangulations (generated with `plantri`) is available in the
file `asset/disk_trgls_3.data`.

# Experimental proof of Proposition 2

The proof is implemented in `tests/eproof.cpp` and can be run with the
executable `experimental_proof`. 

If one uses large valence ranges, the executable may say that a disk
triangulation boundary has not been found (e.g. n=[5,5,5,5,5]). This is because
the disk triangulation is not in `assert/disk_trgls.data`, as the range is too
large, and there would be too many disk triangulations (literally billions,
trillions).


