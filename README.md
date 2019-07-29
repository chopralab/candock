# CANDOCK

This branch contains a version of CANDOCK which has been thoroughly benchmarked.
Please see [this link](https://www.biorxiv.org/content/10.1101/442897v1) for details.

## Dependancies

This version of CANDOCK requires the following software packages to compile:

* A C++-11 compiler (GCC 5.0 or above, Clang, MSVC 2017 or above, etc)
* CMAKE
* The Boost libraries

All other dependancies (GSL and OpenMM) will be downloaded and build by CANDOCK
during compilation.

## How to compile:

cmake . -Bbuild
cd build
make candock

For release compile (with optimization):

cmake3 . -Bbuild -DCMAKE_BUILD_TYPE=Release

