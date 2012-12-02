program main

!*****************************************************************************80
!
!! MAIN is the main program for PRIME_SERIAL_PRB.
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
  implicit none

  integer ( kind = 4 ) n_factor
  integer ( kind = 4 ) n_hi
  integer ( kind = 4 ) n_lo

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRIME_SERIAL_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PRIME_SERIAL library.'

  n_lo = 1
  n_hi = 131072
  n_factor = 2

  call prime_number_sweep ( n_lo, n_hi, n_factor )

  n_lo = 5
  n_hi = 500000
  n_factor = 10

  call prime_number_sweep ( n_lo, n_hi, n_factor )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRIME_SERIAL_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end



