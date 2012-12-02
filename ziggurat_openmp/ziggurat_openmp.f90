program main

!*****************************************************************************80
!
!! MAIN is the main program for ZIGGURAT_OPENMP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2009
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer ( kind = 4 ) proc_num
  integer ( kind = 4 ) thread_num

  call timestamp ( )

  proc_num = omp_get_num_procs ( )
  thread_num = omp_get_max_threads ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZIGGURAT_OPENMP:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of processors available is: ', proc_num
  write ( *, '(a,i8)' ) '  The number of threads available is:    ', thread_num

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZIGGURAT_OPENMP:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests SHR3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 December 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) jsr_value
  real ( kind = 8 ) mega_rate_par
  real ( kind = 8 ) mega_rate_seq
  integer ( kind = 4 ) r
  integer ( kind = 4 ), parameter :: r_num = 1000
  integer ( kind = 4 ), allocatable, dimension ( : ) :: result_par
  integer ( kind = 4 ), allocatable, dimension ( : ) :: result_seq
  integer ( kind = 4 ) s
  integer ( kind = 4 ), parameter :: s_num = 10000
  integer ( kind = 4 ), allocatable, dimension ( : ) :: seed
  integer ( kind = 4 ) shr3
  integer ( kind = 4 ) thread
  integer ( kind = 4 ) thread_num
  real ( kind = 8 ) wtime_par
  real ( kind = 8 ) wtime_seq

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  SHR3 computes random integers.'
  write ( *, '(a)' ) '  Since the output is completely determined'
  write ( *, '(a)' ) '  by the input value of SEED, we can run in'
  write ( *, '(a)' ) '  parallel as long as we make an array of seeds.'
!
!  Set up the SEED array, which will be used for both sequential and
!  parallel computations.
!
!$omp parallel

!$omp master 

  thread_num = omp_get_num_threads ( )

  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  The number of threads is ', thread_num

!$omp end master

!$omp end parallel

  allocate ( seed(0:thread_num-1) )
  allocate ( result_seq(0:thread_num-1) )
  allocate ( result_par(0:thread_num-1) )
!
!  Sequential execution.
!  The sequential execution will only match the parallel execution if we can
!  guarantee that the parallel threads are scheduled to execute the R loop
!  consecutively.
!
  jsr = 123456789

  do thread = 0, thread_num - 1
    seed(thread) = shr3 ( jsr )
  end do

  wtime_seq = omp_get_wtime ( )

  do r = 1, r_num

    thread = mod ( r - 1, thread_num )

    jsr = seed(thread)

    do s = 1, s_num
      jsr_value = shr3 ( jsr )
    end do

    result_seq(thread) = jsr_value
    seed(thread) = jsr

  end do

  wtime_seq = omp_get_wtime ( ) - wtime_seq

  mega_rate_seq = real ( r_num, kind = 8 ) * real ( s_num, kind = 8 ) &
    / wtime_seq / 1000000.0D+00
!
!  Parallel.
!
  jsr = 123456789

  do thread = 0, thread_num - 1
    seed(thread) = shr3 ( jsr )
  end do

  wtime_par = omp_get_wtime ( )

!$omp parallel &
!$omp   shared ( result_par, seed ) &
!$omp   private ( jsr, jsr_value, r, s, thread )

!$omp do schedule ( static, 1 )

  do r = 1, r_num

    thread = omp_get_thread_num ( )

    jsr = seed(thread)

    do s = 1, s_num
      jsr_value = shr3 ( jsr )
    end do

    result_par(thread) = jsr_value
    seed(thread) = jsr

  end do

!$omp end do

!$omp end parallel

  wtime_par = omp_get_wtime ( ) - wtime_par

  mega_rate_par = real ( r_num, kind = 8 ) * real ( s_num, kind = 8 ) &
    / wtime_par / 1000000.0D+00
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Correctness check:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing values sequentially should reach the'
  write ( *, '(a)' ) '  same result as doing it in parallel:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    THREAD    Sequential      Parallel    Difference'
  write ( *, '(a)' ) ' '

  do thread = 0, thread_num - 1
    write ( *, '(2x,i8,2x,i12,2x,i12,2x,i12)' ) &
      thread, result_seq(thread), result_par(thread), &
      result_seq(thread) - result_par(thread)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Efficiency check:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing values in parallel should be faster:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '              Sequential      Parallel'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,2x,g14.6)' ) '      TIME:', wtime_seq,     wtime_par
  write ( *, '(a,g14.6,2x,g14.6)' ) '      RATE:', mega_rate_seq, mega_rate_par

  deallocate ( result_par )
  deallocate ( result_seq )
  deallocate ( seed )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests R4_UNI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 December 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) jsr_value
  real ( kind = 8 ) mega_rate_par
  real ( kind = 8 ) mega_rate_seq
  integer ( kind = 4 ) r
  integer ( kind = 4 ), parameter :: r_num = 1000
  real ( kind = 4 ) r4_uni
  real ( kind = 4 ) r4_value
  real ( kind = 4 ), allocatable, dimension ( : ) :: result_par
  real ( kind = 4 ), allocatable, dimension ( : ) :: result_seq
  integer ( kind = 4 ) s
  integer ( kind = 4 ), parameter :: s_num = 10000
  integer ( kind = 4 ), allocatable, dimension ( : ) :: seed
  integer ( kind = 4 ) shr3
  integer ( kind = 4 ) thread
  integer ( kind = 4 ) thread_num
  real ( kind = 8 ) wtime_par
  real ( kind = 8 ) wtime_seq

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  R4_UNI computes uniformly random single precision real values.'
  write ( *, '(a)' ) '  Since the output is completely determined'
  write ( *, '(a)' ) '  by the input value of SEED, we can run in'
  write ( *, '(a)' ) '  parallel as long as we make an array of seeds.'
!
!  Set up the SEED array, which will be used for both sequential and
!  parallel computations.
!
!$omp parallel

!$omp master 

  thread_num = omp_get_num_threads ( )

  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  The number of threads is ', thread_num

!$omp end master

!$omp end parallel

  allocate ( seed(0:thread_num-1) )
  allocate ( result_seq(0:thread_num-1) )
  allocate ( result_par(0:thread_num-1) )
!
!  Sequential execution.
!  The sequential execution will only match the parallel execution if we can
!  guarantee that the parallel threads are scheduled to execute the R loop
!  consecutively.
!
  jsr = 123456789

  do thread = 0, thread_num - 1
    seed(thread) = shr3 ( jsr )
  end do

  wtime_seq = omp_get_wtime ( )

  do r = 1, r_num

    thread = mod ( r - 1, thread_num )

    jsr = seed(thread)

    do s = 1, s_num
      r4_value = r4_uni ( jsr )
    end do

    result_seq(thread) = r4_value
    seed(thread) = jsr

  end do

  wtime_seq = omp_get_wtime ( ) - wtime_seq

  mega_rate_seq = real ( r_num, kind = 8 ) * real ( s_num, kind = 8 ) &
    / wtime_seq / 1000000.0D+00
!
!  Parallel.
!
  jsr = 123456789

  do thread = 0, thread_num - 1
    seed(thread) = shr3 ( jsr )
  end do

  wtime_par = omp_get_wtime ( )

!$omp parallel &
!$omp   shared ( result_par, seed ) &
!$omp   private ( jsr, r, r4_value, s, thread )

!$omp do schedule ( static, 1 )

  do r = 1, r_num

    thread = omp_get_thread_num ( )

    jsr = seed(thread)

    do s = 1, s_num
      r4_value = r4_uni ( jsr )
    end do

    result_par(thread) = r4_value
    seed(thread) = jsr

  end do

!$omp end do

!$omp end parallel

  wtime_par = omp_get_wtime ( ) - wtime_par

  mega_rate_par = real ( r_num, kind = 8 ) * real ( s_num, kind = 8 ) &
    / wtime_par / 1000000.0D+00
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Correctness check:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing values sequentially should reach the'
  write ( *, '(a)' ) '  same result as doing it in parallel:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    THREAD    Sequential        Parallel      Difference'
  write ( *, '(a)' ) ' '

  do thread = 0, thread_num - 1
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      thread, result_seq(thread), result_par(thread), &
      result_seq(thread) - result_par(thread)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Efficiency check:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing values in parallel should be faster:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '              Sequential      Parallel'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,2x,g14.6)' ) '      TIME:', wtime_seq,     wtime_par
  write ( *, '(a,g14.6,2x,g14.6)' ) '      RATE:', mega_rate_seq, mega_rate_par

  deallocate ( result_par )
  deallocate ( result_seq )
  deallocate ( seed )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests R4_NOR.
!
!  Discussion:
!
!    The arrays FN, KN and WN, once set up by R4_NOR_SETUP, are "read only" when
!    accessed by R4_NOR.  So we only need to have one copy of these arrays,
!    and they can be shared.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 December 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  real ( kind = 4 ) fn(128)
  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) jsr_value
  integer ( kind = 4 ) kn(128)
  real ( kind = 8 ) mega_rate_par
  real ( kind = 8 ) mega_rate_seq
  integer ( kind = 4 ) r
  integer ( kind = 4 ), parameter :: r_num = 1000
  real ( kind = 4 ) r4_nor
  real ( kind = 4 ) r4_value
  real ( kind = 4 ), allocatable, dimension ( : ) :: result_par
  real ( kind = 4 ), allocatable, dimension ( : ) :: result_seq
  integer ( kind = 4 ) s
  integer ( kind = 4 ), parameter :: s_num = 10000
  integer ( kind = 4 ), allocatable, dimension ( : ) :: seed
  integer ( kind = 4 ) shr3
  integer ( kind = 4 ) thread
  integer ( kind = 4 ) thread_num
  real ( kind = 4 ) wn(128)
  real ( kind = 8 ) wtime_par
  real ( kind = 8 ) wtime_seq

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  R4_NOR computes normal random single precision real values.'
  write ( *, '(a)' ) '  Since the output is completely determined'
  write ( *, '(a)' ) '  by the input value of SEED and the tables, we can run in'
  write ( *, '(a)' ) '  parallel as long as we make an array of seeds and share the tables.'
!
!  Set up the SEED array and the tables, which will be used for both sequential and
!  parallel computations.
!
!$omp parallel

!$omp master 

  thread_num = omp_get_num_threads ( )

  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  The number of threads is ', thread_num

!$omp end master

!$omp end parallel

  allocate ( seed(0:thread_num-1) )
  allocate ( result_seq(0:thread_num-1) )
  allocate ( result_par(0:thread_num-1) )

  call r4_nor_setup ( kn, fn, wn )
!
!  Sequential execution.
!  The sequential execution will only match the parallel execution if we can
!  guarantee that the parallel threads are scheduled to execute the R loop
!  consecutively.
!
  jsr = 123456789

  do thread = 0, thread_num - 1
    seed(thread) = shr3 ( jsr )
  end do

  wtime_seq = omp_get_wtime ( )

  do r = 1, r_num

    thread = mod ( r - 1, thread_num )

    jsr = seed(thread)

    do s = 1, s_num
      r4_value = r4_nor ( jsr, kn, fn, wn )
    end do

    result_seq(thread) = r4_value
    seed(thread) = jsr

  end do

  wtime_seq = omp_get_wtime ( ) - wtime_seq

  mega_rate_seq = real ( r_num, kind = 8 ) * real ( s_num, kind = 8 ) &
    / wtime_seq / 1000000.0D+00
!
!  Parallel.
!
  jsr = 123456789

  do thread = 0, thread_num - 1
    seed(thread) = shr3 ( jsr )
  end do

  wtime_par = omp_get_wtime ( )

!$omp parallel &
!$omp   shared ( fn, kn, result_par, seed, wn ) &
!$omp   private ( jsr, r, r4_value, s, thread )

!$omp do schedule ( static, 1 )

  do r = 1, r_num

    thread = omp_get_thread_num ( )

    jsr = seed(thread)

    do s = 1, s_num
      r4_value = r4_nor ( jsr, kn, fn, wn )
    end do

    result_par(thread) = r4_value
    seed(thread) = jsr

  end do

!$omp end do

!$omp end parallel

  wtime_par = omp_get_wtime ( ) - wtime_par

  mega_rate_par = real ( r_num, kind = 8 ) * real ( s_num, kind = 8 ) &
    / wtime_par / 1000000.0D+00
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Correctness check:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing values sequentially should reach the'
  write ( *, '(a)' ) '  same result as doing it in parallel:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    THREAD    Sequential        Parallel      Difference'
  write ( *, '(a)' ) ' '

  do thread = 0, thread_num - 1
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      thread, result_seq(thread), result_par(thread), &
      result_seq(thread) - result_par(thread)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Efficiency check:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing values in parallel should be faster:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '              Sequential      Parallel'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,2x,g14.6)' ) '      TIME:', wtime_seq,     wtime_par
  write ( *, '(a,g14.6,2x,g14.6)' ) '      RATE:', mega_rate_seq, mega_rate_par

  deallocate ( result_par )
  deallocate ( result_seq )
  deallocate ( seed )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests R4_EXP.
!
!  Discussion:
!
!    The arrays FE, KE and WE, once set up by R4_EXP_SETUP, are "read only" when
!    accessed by R4_EXP.  So we only need to have one copy of these arrays,
!    and they can be shared.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 December 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  real ( kind = 4 ) fe(256)
  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) jsr_value
  integer ( kind = 4 ) ke(256)
  real ( kind = 8 ) mega_rate_par
  real ( kind = 8 ) mega_rate_seq
  integer ( kind = 4 ) r
  integer ( kind = 4 ), parameter :: r_num = 1000
  real ( kind = 4 ) r4_exp
  real ( kind = 4 ) r4_value
  real ( kind = 4 ), allocatable, dimension ( : ) :: result_par
  real ( kind = 4 ), allocatable, dimension ( : ) :: result_seq
  integer ( kind = 4 ) s
  integer ( kind = 4 ), parameter :: s_num = 10000
  integer ( kind = 4 ), allocatable, dimension ( : ) :: seed
  integer ( kind = 4 ) shr3
  integer ( kind = 4 ) thread
  integer ( kind = 4 ) thread_num
  real ( kind = 4 ) we(256)
  real ( kind = 8 ) wtime_par
  real ( kind = 8 ) wtime_seq

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  R4_EXP computes exponential random single precision real values.'
  write ( *, '(a)' ) '  Since the output is completely determined'
  write ( *, '(a)' ) '  by the input value of SEED and the tables, we can run in'
  write ( *, '(a)' ) '  parallel as long as we make an array of seeds and share the tables.'
!
!  Set up the SEED array and the tables, which will be used for both sequential and
!  parallel computations.
!
!$omp parallel

!$omp master 

  thread_num = omp_get_num_threads ( )

  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  The number of threads is ', thread_num

!$omp end master

!$omp end parallel

  allocate ( seed(0:thread_num-1) )
  allocate ( result_seq(0:thread_num-1) )
  allocate ( result_par(0:thread_num-1) )

  call r4_exp_setup ( ke, fe, we )
!
!  Sequential execution.
!  The sequential execution will only match the parallel execution if we can
!  guarantee that the parallel threads are scheduled to execute the R loop
!  consecutively.
!
  jsr = 123456789

  do thread = 0, thread_num - 1
    seed(thread) = shr3 ( jsr )
  end do

  wtime_seq = omp_get_wtime ( )

  do r = 1, r_num

    thread = mod ( r - 1, thread_num )

    jsr = seed(thread)

    do s = 1, s_num
      r4_value = r4_exp ( jsr, ke, fe, we )
    end do

    result_seq(thread) = r4_value
    seed(thread) = jsr

  end do

  wtime_seq = omp_get_wtime ( ) - wtime_seq

  mega_rate_seq = real ( r_num, kind = 8 ) * real ( s_num, kind = 8 ) &
    / wtime_seq / 1000000.0D+00
!
!  Parallel.
!
  jsr = 123456789

  do thread = 0, thread_num - 1
    seed(thread) = shr3 ( jsr )
  end do

  wtime_par = omp_get_wtime ( )

!$omp parallel &
!$omp   shared ( fe, ke, result_par, seed, we ) &
!$omp   private ( jsr, r, r4_value, s, thread )

!$omp do schedule ( static, 1 )

  do r = 1, r_num

    thread = omp_get_thread_num ( )

    jsr = seed(thread)

    do s = 1, s_num
      r4_value = r4_exp ( jsr, ke, fe, we )
    end do

    result_par(thread) = r4_value
    seed(thread) = jsr

  end do

!$omp end do

!$omp end parallel

  wtime_par = omp_get_wtime ( ) - wtime_par

  mega_rate_par = real ( r_num, kind = 8 ) * real ( s_num, kind = 8 ) &
    / wtime_par / 1000000.0D+00
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Correctness check:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing values sequentially should reach the'
  write ( *, '(a)' ) '  same result as doing it in parallel:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    THREAD    Sequential        Parallel      Difference'
  write ( *, '(a)' ) ' '

  do thread = 0, thread_num - 1
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      thread, result_seq(thread), result_par(thread), &
      result_seq(thread) - result_par(thread)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Efficiency check:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing values in parallel should be faster:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '              Sequential      Parallel'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,2x,g14.6)' ) '      TIME:', wtime_seq,     wtime_par
  write ( *, '(a,g14.6,2x,g14.6)' ) '      RATE:', mega_rate_seq, mega_rate_par

  deallocate ( result_par )
  deallocate ( result_seq )
  deallocate ( seed )

  return
end
function r4_exp ( jsr, ke, fe, we )

!*****************************************************************************80
!
!! R4_EXP returns an exponentially distributed single precision real value.
!
!  Discussion:
!
!    The underlying algorithm is the ziggurat method.
!
!    Before the first call to this function, the user must call R4_EXP_SETUP
!    to determine the values of KE, FE and WE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) JSR, the seed.
!
!    Input, integer ( kind = 4 ) KE(256), data computed by R4_EXP_SETUP.
!
!    Input, real ( kind = 4 ) FE(256), WE(256), data computed by R4_EXP_SETUP.
!
!    Output, real ( kind = 4 ) R4_EXP, an exponentially distributed 
!    random value.
!
  implicit none

  real ( kind = 4 ) fe(256)
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) jz
  integer ( kind = 4 ) ke(256)
  real ( kind = 4 ) r4_exp
  real ( kind = 4 ) r4_uni
  integer ( kind = 4 ) shr3
  real ( kind = 4 ) value
  real ( kind = 4 ) we(256)
  real ( kind = 4 ) x

  jz = shr3 ( jsr )
  iz = iand ( jz, 255 )

  if ( abs ( jz  ) < ke(iz+1) ) then

    value = real ( abs ( jz ), kind = 4 ) * we(iz+1)

  else

    do

      if ( iz == 0 ) then
        value = 7.69711E+00 - log ( r4_uni ( jsr ) )
        exit
      end if

      x = real ( abs ( jz ), kind = 4 ) * we(iz+1)

      if ( fe(iz+1) + r4_uni ( jsr ) * ( fe(iz) - fe(iz+1) ) < exp ( - x ) )  then
        value = x
        exit
      end if

      jz = shr3 ( jsr )
      iz = iand ( jz, 255 )

      if ( abs ( jz ) < ke(iz+1) ) then
        value = real ( abs ( jz ), kind = 4 ) * we(iz+1)
        exit
      end if

    end do
    
  end if

  r4_exp = value

  return
end
subroutine r4_exp_setup ( ke, fe, we )

!*****************************************************************************80
!
!! R4_EXP_SETUP sets data needed by R4_EXP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KE(256), data needed by R4_EXP.
!
!    Output, real ( kind = 4 ) FE(256), WE(256), data needed by R4_EXP.
!
  implicit none

  real ( kind = 8 ) de
  real ( kind = 4 ) fe(256)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ke(256)
! real ( kind = 8 ), parameter :: m2 = 4294967296.0D+00
  real ( kind = 8 ), parameter :: m2 = 2147483648.0D+00
  real ( kind = 8 ) q
  real ( kind = 8 ) te
  real ( kind = 8 ), parameter :: ve = 3.949659822581572D-03
  real ( kind = 4 ) we(256)

  de = 7.697117470131487D+00
  te = 7.697117470131487D+00

  q = ve / exp ( - de )

  ke(1) = int ( ( de / q ) * m2 )
  ke(2) = 0

  we(1) = real ( q / m2, kind = 4 )
  we(256) = real ( de / m2, kind = 4 )

  fe(1) = 1.0E+00
  fe(256) = real ( exp ( - de ), kind = 4 )

  do i = 255, 2, -1
    de = - log ( ve / de + exp ( - de ) )
    ke(i+1) = int ( ( de / te ) * m2 )
    te = de
    fe(i) = real ( exp ( - de ), kind = 4 )
    we(i) = real ( de / m2, kind = 4 )
  end do

  return
end
function r4_nor ( jsr, kn, fn, wn )

!*****************************************************************************80
!
!! R4_NOR returns a normally distributed single precision real value.
!
!  Discussion:
!
!    The value returned is generated from a distribution with mean 0 and 
!    variance 1.
!
!    The underlying algorithm is the ziggurat method.
!
!    Before the first call to this function, the user must call R4_NOR_SETUP
!    to determine the values of KN, FN and WN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) JSR, the seed.
!
!    Input, integer ( kind = 4 ) KN(128), data computed by R4_NOR_SETUP.
!
!    Input, real ( kind = 4 ) FN(128), WN(128), data computed by R4_NOR_SETUP.
!
!    Output, real ( kind = 4 ) R4_NOR, a normally distributed random value.
!
  implicit none

  real ( kind = 4 ) fn(128)
  integer ( kind = 4 ) hz
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) kn(128)
  real ( kind = 4 ), parameter :: r = 3.442620E+00
  real ( kind = 4 ) r4_nor
  real ( kind = 4 ) r4_uni
  integer ( kind = 4 ) shr3
  real ( kind = 4 ) value
  real ( kind = 4 ) wn(128)
  real ( kind = 4 ) x
  real ( kind = 4 ) y

  hz = shr3 ( jsr )
  iz = iand ( hz, 127 )

  if ( abs ( hz ) < kn(iz+1) ) then

    value = real ( hz, kind = 4 ) * wn(iz+1)

  else

    do

      if ( iz == 0 ) then

        do
          x = - 0.2904764E+00 * log ( r4_uni ( jsr ) )
          y = - log ( r4_uni ( jsr ) )
          if ( x * x <= y + y ) then
            exit
          end if
        end do

        if ( hz <= 0 ) then
          value = - r - x
        else
          value = + r + x
        end if

        exit

      end if

      x = real ( hz, kind = 4 ) * wn(iz+1)

      if ( fn(iz+1) + r4_uni ( jsr ) * ( fn(iz) - fn(iz+1) ) &
         < exp ( - 0.5E+00 * x * x ) ) then
        value = x
        exit
      end if

      hz = shr3 ( jsr )
      iz = iand ( hz, 127 )

      if ( abs ( hz ) < kn(iz+1) ) then
        value = real ( hz, kind = 4 ) * wn(iz+1)
        exit
      end if

    end do

  end if

  r4_nor = value

  return
end
subroutine r4_nor_setup ( kn, fn, wn )

!*****************************************************************************80
!
!! R4_NOR_SETUP sets data needed by R4_NOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KN(128), data needed by R4_NOR.
!
!    Output, real ( kind = 4 ) FN(128), WN(128), data needed by R4_NOR.
!
  implicit none

  real ( kind = 8 ) dn
  real ( kind = 4 ) fn(128)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kn(128)
  real ( kind = 8 ), parameter :: m1 = 2147483648.0D+00
  real ( kind = 8 ) q
  real ( kind = 8 ) tn
  real ( kind = 8 ), parameter :: vn = 9.91256303526217D-03
  real ( kind = 4 ) wn(128)

  dn = 3.442619855899D+00
  tn = 3.442619855899D+00

  q = vn / exp ( - 0.5D+00 * dn * dn )

  kn(1) = int ( ( dn / q ) * m1 )
  kn(2) = 0

  wn(1) = real ( q / m1, kind = 4 )
  wn(128) = real ( dn / m1, kind = 4 )

  fn(1) = 1.0E+00
  fn(128) = real ( exp ( - 0.5D+00 * dn * dn ), kind = 4 )

  do i = 127, 2, -1
    dn = sqrt ( - 2.0D+00 * log ( vn / dn + exp ( - 0.5D+00 * dn * dn ) ) )
    kn(i+1) = int ( ( dn / tn ) * m1 )
    tn = dn
    fn(i) = real ( exp ( - 0.5D+00 * dn * dn ), kind = 4 )
    wn(i) = real ( dn / m1, kind = 4 )
  end do

  return
end
function r4_uni ( jsr )

!*****************************************************************************80
!
!! R4_UNI returns a uniformly distributed real value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) JSR, the seed.
!
!    Output, real ( kind = 4 ) R4_UNI, a uniformly distributed random value in
!    the range [0,1].
!
  implicit none

  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) jsr_input
  real ( kind = 4 ) r4_uni

  jsr_input = jsr

  jsr = ieor ( jsr, ishft ( jsr,   13 ) )
  jsr = ieor ( jsr, ishft ( jsr, - 17 ) )
  jsr = ieor ( jsr, ishft ( jsr,    5 ) )

! r4_uni = 0.5E+00 + 0.2328306E-09 * real ( jsr_input + jsr, kind = 4 )

  r4_uni = 0.5E+00 + real ( jsr_input + jsr, kind = 4 ) &
    / real ( 65536, kind = 4 ) / real ( 65536, kind = 4 )

  return
end
function shr3 ( jsr )

!*****************************************************************************80
!
!! SHR3 evaluates the SHR3 generator for integers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2008
!
!  Author:
!
!    Original C version by George Marsaglia, Wai Wan Tsang
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    George Marsaglia, Wai Wan Tsang,
!    The Ziggurat Method for Generating Random Variables,
!    Journal of Statistical Software,
!    Volume 5, Number 8, October 2000, seven pages.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) JSR, the seed.
!
!    Output, integer ( kind = 4 ) SHR3, the value of the SHR3 generator.
!
  implicit none

  integer ( kind = 4 ) jsr
  integer ( kind = 4 ) jsr_input
  integer ( kind = 4 ) shr3

  jsr_input = jsr

  jsr = ieor ( jsr, ishft ( jsr,   13 ) )
  jsr = ieor ( jsr, ishft ( jsr, - 17 ) )
  jsr = ieor ( jsr, ishft ( jsr,    5 ) )

  shr3 = jsr_input + jsr

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
