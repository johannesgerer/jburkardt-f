program main

!*****************************************************************************80
!
!  Purpose:
!
!    MAIN is the main program for TEST01.
!
!  Discussion:
!
!    TEST02 has some uninitialized data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 May 2011
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A sample code for analysis by VALGRIND.'

  call junk_data ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine junk_data ( )

!*****************************************************************************80
!
!  Purpose:
!
!    JUNK_DATA has some uninitialized variables.
!
!  Discussion:
!
!    VALGRIND's MEMCHECK program monitors uninitialized variables, but does
!    not complain unless such a variable is used in a way that means its
!    value affects the program's results, that is, the value is printed,
!    or computed with.  Simply copying the unitialized data to another variable
!    is of no concern.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2011
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable :: x(:)

  allocate ( x(1:10) )
!
!  X = { 0, 1, 2, 3, 4, ?a, ?b, ?c, ?d, ?e }.
!
  do i = 1, 5
    x(i) = i - 1
  end do
!
!  Copy some values.
!  X = { 0, 1, ?c, 3, 4, ?b, ?b, ?c, ?d, ?e }.
!
  x(3) = x(8)
  x(6) = x(7)
!
!  Modify some uninitialized entries.
!  Memcheck doesn't seem to care about this.
!
  do i = 1, 10
    x(i) = 2 * x(i)
  end do
!
!  Print X.
!
  do i = 1, 10
    write ( *, '(2x,i2,2x,i2)' ) i, x(i)
  end do

  deallocate ( x )

  return
end
