# Divsufsort
mkdir libdivsufsort/build
cd libdivsufsort/build
cmake -DCMAKE_INSTALL_PREFIX=.. ..
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
