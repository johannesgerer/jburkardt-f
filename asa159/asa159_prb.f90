program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA159_PRB.
!
!  Discussion:
!
!    ASA159_PRB tests the routines in ASA159.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA159_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA159 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA159_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests RCONT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ), save, dimension ( n ) :: c = (/ 2, 2, 2, 2, 1 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  logical              key
  integer ( kind = 4 ), parameter :: ntest = 10
  integer ( kind = 4 ), save, dimension ( m ) :: r = (/ 3, 2, 2, 1, 1 /)
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  RCONT2 constructs a random matrix with'
  write ( *, '(a)' ) '  given row and column sums.'

  call i4vec_print ( m, r, '  The rowsum vector:' )
  call i4vec_print ( n, c, '  The columnsum vector:' )

  key = .false.

  do i = 1, ntest

    call rcont2 ( m, n, r, c, key, seed, a, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  RCONT2 returned error flag IERROR = ', ierror
      return
    end if

    call i4mat_print ( m, n, a, '  The rowcolsum matrix:' )

  end do

  return
end
