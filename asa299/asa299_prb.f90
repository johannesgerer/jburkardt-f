program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA299_PRB.
!
!  Discussion:
!
!    ASA299_PRB calls the ASA299 test routines.
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA299_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA299 library.'

  call test01 ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA299_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests SIMPLEX_LATTICE_POINT_NEXT.
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 4

  integer i
  logical more
  integer, parameter :: t = 4
  integer x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  SIMPLEX_LATTICE_POINT_NEXT generates lattice points'
  write ( *, '(a)' ) '  in the simplex'
  write ( *, '(a)' ) '    0 <= X'
  write ( *, '(a)' ) '    sum ( X(1:N) ) <= T'
  write ( *, '(a,i8)' ) '  Here N = ', n
  write ( *, '(a,i8)' ) '  and T =  ', t
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Index        X(1)      X(2)      X(3)      X(4)'
  write ( *, '(a)' ) ' '

  more = .false.

  i = 0

  do

    call simplex_lattice_point_next ( n, t, more, x )

    i = i + 1

    write ( *, '(2x,i8,2x,4(2x,i8))' ) i, x(1:n)

    if ( .not. more )  then
      exit
    end if

  end do
 
  return
end
