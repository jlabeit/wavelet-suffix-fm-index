for algo in ddWT levelWT recWT sdslWT serialWT
do
    cd $algo
    make clean
    make
    cd ..
done
