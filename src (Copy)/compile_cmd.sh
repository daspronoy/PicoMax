#compile_cmd.sh
#!/bin/bash
BASE=/home/sauron/projects/picomax


# -ffast-math
# -funsafe-math-optimizations
# -fno-signed-zeros
# -ffinite-math-only
# -fno-math-errno -fno-trapping-math

#-L/usr/lib/x86_64-linux-gnu/hdf5/serial \
# -I /usr/include/hdf5/serial \
# -I ${BASE}/include/hdf5-1.14.3/src -I ${BASE}/include/hdf5 
g++ -O3 -fopenmp -ffast-math \
    -I ${BASE}/include/eigen-3.4.0\
    -o picomax pmx.cpp pmx_io.cpp pmx_chi.cpp pmx_basis.cpp pmx_epm.cpp