for algo in rankSerial rank selectSerial select
do
    cd $algo
    rm runTests.py
    make clean
    make
    cd ..
done
