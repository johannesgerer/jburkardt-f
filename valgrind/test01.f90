program main

!*****************************************************************************80
!
!  Purpose:
!
!    MAIN is the main program for TEST01.
!
!  Discussion:
!
!    TEST01 calls F, which has a memory "leak".  This memory leak can be
!    detected by VALGRID.
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

  integer ( kind = 4 ) n

  n = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  A sample code for analysis by VALGRIND.'

  call f ( n )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine f ( n )

!*****************************************************************************80
!
!  Purpose:
!
!    F computes N+1 entries of the Fibonacci sequence.
!
!  Discussion:
!
!    Unfortunately, F only allocates space for N entries.  Hence, the
!    assignment of a value to the N+1 entry causes a memory leak.
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
  integer ( kind = 4 ) n
  integer ( kind = 4 ), allocatable :: x(:)

  allocate ( x(1:n) )

  x(1) = 1
  write ( *, '(2x,i2,2x,i2)' ) 1, x(1)

  x(2) = 1
  write ( *, '(2x,i2,2x,i2)' ) 2, x(2)

  do i = 3, n + 1
    x(i) = x(i-1) + x(i-2)
    write ( *, '(2x,i2,2x,i2)' ) i, x(i)
  end do

  deallocate ( x )

  return
end

