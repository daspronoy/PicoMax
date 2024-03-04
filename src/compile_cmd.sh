
export EIGEN="../include/eigen-3.4.0"


g++ -O3 -fopenmp -DMKL_ILP64 -m64 -I "${MKLROOT}/include" -I ${EIGEN} pmx.cpp -o picomax

#g++ -O3 -fopenmp -I ../include/eigen-3.4.0 pmx.cpp -o picomax