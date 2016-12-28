# Divsufsort
mkdir libdivsufsort/build
cd libdivsufsort/build
cmake -DCMAKE_INSTALL_PREFIX=..  -DBUILD_DIVSUFSORT64=ON ..
make
make install
cd ../../

# Parallel-divsufsort
mkdir parallel-divsufsort/build
cd parallel-divsufsort/build
cmake -DCMAKE_INSTALL_PREFIX=.. ..
make
make install
cd ../../

# All suffix arrays.
cd suffixArray
sh install.sh
cd ../
