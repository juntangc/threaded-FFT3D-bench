!!!! -*- Mode: F90 -*- !!!!


program fftw_test_3d_multithread

  include 'fftw3.f'
  include 'mpif.h'

  integer*8          :: plan
  integer*4          :: iRet
  character(LEN=100) :: option

  integer*4          :: N_THREADS
  integer*4          :: VAR
  integer*4          :: N_ITER
  integer*4          :: i
  integer*4          :: dim_x, dim_y, dim_z
  integer*4          :: isign

  complex*16   , dimension(:,:,:), allocatable :: in
!$    INTEGER THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS
  integer*4 rank, size, ierror
  real*8                :: ttt0, dt, sumdt
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

narg = COMMAND_ARGUMENT_COUNT()
  if (narg < 4) then
     write(*,*) ""
     write(*,*) "Test/bench for multithread FFTW speed"
     write(*,*) "Note: you must have compiled FFTW with './configure --enable-threads'"
     write(*,*) ""
     write(*,*) "fftw_test_3d_mt_double <dim_x> <dim_y> <dim_z> <nthreads> [<niter>]"
     write(*,*) " l : number of FFT points"
     write(*,*) " n : number of SMP threads"
     write(*,*) " i : how many FFT calls to make and average their exec times"
     stop
  end if 

  ! retrieve first 3 input parameters : dim_x, dim_y, dim_z
  call GET_COMMAND_ARGUMENT(1,option)
  read(option,*)  VAR
  dim_x    = VAR
  call GET_COMMAND_ARGUMENT(2,option)
  read(option,*)  VAR
  dim_y    = VAR
  call GET_COMMAND_ARGUMENT(3,option)
  read(option,*)  VAR
  dim_z    = VAR

  ! 4th input parameter : number of threads
  call GET_COMMAND_ARGUMENT(4,option)
  read(option,*) N_THREADS

  ! retrieve 5rd input parameter : number of iteration
  N_ITER = 1
  if (narg>=5) then
     call GET_COMMAND_ARGUMENT(5,option)
     read(option,*) N_ITER
  end if

  if (dim_x < 8 .or. dim_x > 1024 .or. dim_y < 8 .or. dim_y > 1024 .or. dim_z < 8 .or. dim_z > 1024) then
     write(*,*) "FFT dimension not between 8 and 1024"
     stop
  end if
  
  if (N_THREADS < 1) N_THREADS=1
  if (N_THREADS > 16) N_THREADS=16

if(rank.eq.0) then
  write(*,*) "Using ", dim_x, "*", dim_y, "*", dim_z, " as fft length"
  write(*,*) "Using ", N_THREADS, " threads"
  write(*,*) "Using ", N_ITER, " iteration"
endif

  ! allocate memory and initialize
  allocate (in  (     dim_x     , dim_y , dim_z) )

  ! 
  DO k = 1,dim_z
     DO j = 1,dim_y
        DO i = 1,dim_x
           in(i,j,k) = (1.d0, 1.d0)
        END DO
     END DO
  END DO

  ! initialize FFTW
  call dfftw_init_threads(iRet)
if(rank.eq.0) then
  write(*,*) 'fftw init thread status : ',iRet
endif

  call dfftw_plan_with_nthreads(N_THREADS)

  !   initialize the plan: this can take a long time !
  
if(rank.eq.0) then
  write(*,*) 'Creating complex-to-complex 3D DFT plan'
endif
  isgin = 0

!$!=======================================================================
!$!  initialise openMP FFT, has to be done here and not in main.F in
!$!  in general, since the FFTs are called outside and inside openMP
!$!  parallel regions.
!$!=======================================================================
!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

   !  now run one or more FFT interations 
   ttt0 = mpi_wtime()

   isign = 0
   DO i = 1,N_ITER
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS) 
      if (isign.eq.0) then
         call dfftw_plan_dft_3d(plan, dim_x, dim_y, dim_z, in, in, FFTW_FORWARD, FFTW_ESTIMATE)
         isign = 1
      else
         call dfftw_plan_dft_3d(plan, dim_x, dim_y, dim_z, in, in, FFTW_BACKWARD, FFTW_ESTIMATE)
         isign = 0
      endif
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      call dfftw_execute_dft(plan, in, in)
#if !defined(fftw_cache_plans) && !defined(__NEC_TUNE__)
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      call dfftw_destroy_plan(plan)
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
#endif // fftw_cache_plans or __NEC_TUNE__
   end DO

   dt = (mpi_wtime()-ttt0)/N_ITER
   call mpi_reduce(dt, sumdt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
  
if(rank.eq.0) then
   write(*,*) 'Plan exec time        : ',sumdt/size,' sec (average over ', N_ITER,' iteration(s))'
endif

   deallocate(in)
   
   call dfftw_cleanup_threads();
   
call MPI_FINALIZE(ierror)
end program fftw_test_3d_multithread

