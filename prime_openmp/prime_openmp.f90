program main

!*****************************************************************************80
!
!! MAIN is the main program for PRIME_NUMBER_OPENMP.
!
!  Discussion:
!
!    This program calls a version of PRIME_NUMBER that includes
!    OpenMP directives for parallel processing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2009
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer ( kind = 4 ) n_factor
  integer ( kind = 4 ) n_hi
  integer ( kind = 4 ) n_lo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRIME_NUMBER_OPENMP'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  Number of processors available = ', omp_get_num_procs ( )
  write ( * ,'(a,i8)' ) &
    '  Number of threads =              ', omp_get_max_threads ( )

  n_lo = 1
  n_hi = 131072
  n_factor = 2

  call prime_number_sweep_openmp ( n_lo, n_hi, n_factor )

  n_lo = 5
  n_hi = 500000
  n_factor = 10

  call prime_number_sweep_openmp ( n_lo, n_hi, n_factor )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRIME_NUMBER_OPENMP'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine prime_number_sweep_openmp ( n_lo, n_hi, n_factor )

!*****************************************************************************80
!
!! PRIME_NUMBER_SWEEP_OPENMP does repeated calls to PRIME_NUMBER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LO, the first value of N.
!
!    Input, integer ( kind = 4 ) N_HI, the last value of N.
!
!    Input, integer ( kind = 4 ) N_FACTOR, the factor by which to increase N after
!    each iteration.
!
  use omp_lib

  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_factor
  integer ( kind = 4 ) n_hi
  integer ( kind = 4 ) n_lo
  integer ( kind = 4 ) primes
  real    ( kind = 8 ) wtime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRIME_NUMBER_SWEEP_OPENMP'
  write ( *, '(a)' ) '  Call PRIME_NUMBER to count the primes from 1 to N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        Pi          Time'
  write ( *, '(a)' ) ' '

  n = n_lo

  do while ( n <= n_hi )

    wtime = omp_get_wtime ( )

    call prime_number ( n, primes )

    wtime = omp_get_wtime ( ) - wtime

    write ( *, '(2x,i8,2x,i8,g14.6)' ) n, primes, wtime

    n = n * n_factor

  end do
 
  return
end
subroutine prime_number ( n, total )

!*****************************************************************************80
!
!! PRIME_NUMBER returns the number of primes between 1 and N.
!
!  Discussion:
!
!    A naive algorithm is used.
!
!    Mathematica can return the number of primes less than or equal to N
!    by the command PrimePi[N].
!
!                N  PRIME_NUMBER
!
!                1           0
!               10           4
!              100          25
!            1,000         168
!           10,000       1,229
!          100,000       9,592
!        1,000,000      78,498
!       10,000,000     664,579
!      100,000,000   5,761,455
!    1,000,000,000  50,847,534
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum number to check.
!
!    Output, integer ( kind = 4 ) TOTAL, the number of prime numbers up to N.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) total

  total = 0

!$omp parallel &
!$omp shared ( n ) &
!$omp private ( i, j, prime )

!$omp do reduction ( + : total )

  do i = 2, n

    prime = 1

    do j = 2, i - 1
      if ( mod ( i, j ) == 0 ) then
        prime = 0
        exit
      end if
    end do

    total = total + prime

  end do

!$omp end do

!$omp end parallel

  return
end
