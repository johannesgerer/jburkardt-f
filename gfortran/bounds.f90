program main

!*****************************************************************************80
!
!! MAIN is the main program for BOUNDS.
!
!  Discussion:
!
!    BOUNDS is a FORTRAN90 program in which an illegal array reference is made.
!
!    The GFORTRAN compiler switch "-fbounds-check" will catch this.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BOUNDS:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program uses an illegal memory reference.'
  write ( *, '(a)' ) '  In this case, an array element A(11) is read,'
  write ( *, '(a)' ) '  although the array is only dimensioned for size 10.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compilation with the GFORTRAN switch -fbounds-check'
  write ( *, '(a)' ) '  will generate a run-time warning.'
!
!  Initialize.
!
  do i = 1, 10
    a(i) = i
  end do
!
!  Add neighbor.
!  (And accidentally invoke nonexistent A(11)!)
!
  do i = 1, 10
    a(i) = a(i) + a(i+1)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BOUNDS:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
