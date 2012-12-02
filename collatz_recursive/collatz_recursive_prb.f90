program main

!*****************************************************************************80
!
!! MAIN is the main program for COLLATZ_RECURSIVE_PRB.
!
!  Discussion:
!
!    COLLATZ_RECURSIVE_PRB tests the COLLATZ_RECURSIVE library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLLATZ_RECURSIVE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the COLLATZ_RECURSIVE library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLLATZ_RECURSIVE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests COLLATZ_PATH;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) n
  integer ( kind = 4 ) :: n_test(test_num) = (/ 7, 8, 9, 10, 600 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  COLLATZ_PATH prints the members of a Collatz path.'

  do test = 1, test_num
    n = n_test(test)
    write ( *, '(a)' ) ' '
    write ( *, '(2x,i8,a)' ) n, ' is the starting point.'
    call collatz_path ( n )
  end do

  return
end
