for algo in divsufsort parallelDivsufsort parallelKS parallelRange pScan
do
    cd $algo
    make clean
    make
    cd ..
done
