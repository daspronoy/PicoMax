#!/bin/bash

g++ -O3 -fopenmp -I ../include/eigen-3.4.0 -o picomax pmx.cpp chi.cpp basis.cpp epm.cpp io.cpp