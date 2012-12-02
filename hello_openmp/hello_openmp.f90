program main

!*****************************************************************************80
!
!! MAIN is the main program for HELLO.
!
!  Discussion:
!
!    HELLO is a "Hello, World" program for OpenMP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 June 2010
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer id
  double precision wtime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HELLO_OPENMP'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'

  wtime = omp_get_wtime ( )

  write ( *, '(a,i8)' ) &
    '  The number of processors available = ', omp_get_num_procs ( )
  write ( *, '(a,i8)' ) &
    '  The number of threads available    = ', omp_get_max_threads ( )
!
!  OUTSIDE THE PARALLEL REGION, have each thread say hello (there's only 1).
!
  id = omp_get_thread_num ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  OUTSIDE the parallel region.'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8,a,i8)' ) '  HELLO from process ', id

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Going INSIDE the parallel region:'
  write ( *, '(a)' ) ' '
!
!  INSIDE THE PARALLEL REGION, have each thread say hello.
!
!$omp parallel &
!$omp private ( id )
  id = omp_get_thread_num ( )

  write ( *, '(a,i8,a,i8)' ) '  HELLO from process ', id

!$omp end parallel
!
!  Finish up by measuring the elapsed time.
!
  wtime = omp_get_wtime ( ) - wtime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Back OUTSIDE the parallel region.'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HELLO_OPENMP'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Elapsed wall clock time = ', wtime

  stop
end
