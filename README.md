# threaded-FFT3D-bench
3D FFT benchmark
This benchmark is designed to measure FFT3D performance using various FFT libraries: FFTW, MKL, ArmPL. Modify the `Makefile` or use the following commands to compile.

1. use ArmPL, gcc compiler on arm64
```
mpif90 -cpp -O3 -Wall -I${ARMPL_INCLUDES} -c fftw_test_3d_multithread.f90
mpif90 -o fftw_test_3d_multithread fftw_test_3d_multithread.o -O3 -L${ARMPL_LIBRARIES} -larmpl_mp -fopenmp
```

2. use FFTW, gcc compiler on arm64
```
mpif90 -cpp -O3 -Wall -I${FFTW_ROOT}/include -c fftw_test_3d_multithread.f90
mpif90 -o fftw_test_3d_multithread fftw_test_3d_multithread.o -O3 -L${FFTW_ROOT}/lib -lfftw3_threads -lfftw3 
```

3. use MKL, ifort compiler on x86
```
mpiifort -cpp -O3 -Wall -I${FFTW_ROOT}/include -c fftw_test_3d_multithread.f90
mpiifort -o fftw_test_3d_multithread fftw_test_3d_multithread.o -O3 -L${FFTW_ROOT}/lib -lfftw3_threads -lfftw3 
```

4. use FFTW, ifort compiler on x86
```
mpiifort -cpp -O3 -Wall -I${MKL_ROOT}/include -c fftw_test_3d_multithread.f90
mpiifort -o fftw_test_3d_multithread fftw_test_3d_multithread.o -O3 -qmkl
```

commands to run the benchmark (3D complex-to-complex FFT with the dimension of 70 x 70 x 210) with 8 ranks per node and 8 threads per rank
Open MPI, hpc7g
```
/shared/tools/threaded_fft_bench$ mpirun -np 8 --map-by node:PE=8 --bind-to core -x OMP_NUM_THREADS=8 ./fftw_test_3d_multithread-fftw 70 70 210 8 100
```

Intel MPI, hpc6id
```
 mpirun -np 8 -genv OMP_NUM_THREADS=8 -genv I_MPI_PIN_DOMAIN=omp ./fftw_test_3d_multithread 70 70 210 8 100
```
