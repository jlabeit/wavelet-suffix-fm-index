Description
==========
This repository contains the experimental setup from
> Julian Labeit, Julian Shun, and Guy E. Blelloch. Parallel Lightweight Wavelet Tree, Suffix Array and FM-Index Construction. DCC 2015.


Installation
==========
Currently the code is only tested on Ubuntu machine with GCC 5.2.0 and cilkplus for parallelization.
- To enable parallelization set the environment variable GCILK=1 before compilation.
- To enable the usage of 64 bits integers (for example when building large wavelet trees) set the environment variable LONG=1 before compilation.
