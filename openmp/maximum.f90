program main

!*****************************************************************************80
!
!! MAIN is the main program for MAXIMUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 June 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer, parameter :: n = 50

  integer ( kind = 4 ) i
  integer ( kind = 4 ) proc_num
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) thread_num
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) x_max1
  real    ( kind = 8 ) x_max2
  real    ( kind = 8 ) x_max3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MAXIMUM:'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  OpenMP allows FORTRAN programs to compute a maximum'
  write ( *, '(a)' ) '  as a reduction variable.  How does this work exactly?'

  write ( *, '(a,i8)' ) '  The number of processors available = ', omp_get_num_procs ( )
  write ( *, '(a,i8)' ) '  The number of threads available    = ', omp_get_max_threads ( )
!
!  Set X to some "random" but not too wild values.
!
  seed = 1325

  do i = 1, n
    seed = mod ( ( 3125 * seed ), 65536 )
    x(i) = ( seed - 32768.0 ) / 16384.0
  end do
!
!  Determine X_MAX1 the old fashioned way.
!
  x_max1 = - huge ( x_max1 )
  do i = 1, n
    if ( x_max1 < x(i) ) then
      x_max1 = x(i)
    end if
  end do
!
!  Determine X_MAX2 in a parallel loop explicitly invoking the MAX function.  
!  Note also that we do NOT initialize X_MAX2.
!
!$omp parallel &
!$omp   shared ( x ) &
!$omp   private ( i ) 

!$omp do reduction ( max: x_max2 )
  do i = 1, n
    x_max2 = max ( x_max2, x(i) )
  end do
!$omp end do

!$omp end parallel
!
!  Determine X_MAX3 in a parallel loop but don't actually call MAX, do it yourself. 
!
  x_max3 = - huge ( x_max3 )

!$omp parallel &
!$omp   shared ( x ) &
!$omp   private ( i )

!$omp do reduction ( max: x_max3 )
  do i = 1, n
    if ( x_max3 < x(i) ) then
      x_max3 = x(i)
    end if
  end do
!$omp end do

!$omp end parallel

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X_MAX should be   ', x_max1
  write ( *, '(a,g14.6)' ) '  Computed X_MAX2 = ', x_max2
  write ( *, '(a,g14.6)' ) '  Computed X_MAX3 = ', x_max3
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MAXIMUM:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
