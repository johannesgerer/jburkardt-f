program main

!*******************************************************************************
!
!! MAIN is the main program for the F90_SIMPLE example.
!
!  Discussion:
!
!    F90_SIMPLE is a simple FORTRAN90 program intended to demonstrate
!    the use of a makefile.  The program is divided into three files,
!    containing the main program, a midpoint rule method for approximating
!    the integral of a function, and a routine to evaluate the function.
!
!    The makefile contains the "knowledge" necessary to appropriately
!    compile the files, join them together, run the program and
!    even discard all the auxilliary files afterwards.
!
!  Modified:
!
!    04 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Local, real A, B, the endpoints of the interval of integration.
!
!    Local, real external F, the name of the function to be integrated.
!
!    Local, integer INT_NUM, the number of intervals to be used.
!
!    Local, real QUAD, the approximate value of the integral.
!
  implicit none

  real a
  real b
  real, external :: f
  integer int_num
  real quad
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'F90_SIMPLE'
  write ( *, '(a)' ) '  A simple FORTRAN90 program to demonstrate'
  write ( *, '(a)' ) '  the use of makefiles.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate the integral from 0 to 1000, of'
  write ( *, '(a)' ) '  F(T) = (4+T/365+1/2 sin(pi*T/91)) * (2+exp(-sin(2*pi*T)))'
  write ( *, '(a)' ) '  a function which models daily power consumption.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  subroutine midpoint ( a, b, f, int_num, quad )'
  write ( *, '(a)' ) '  estimates the integral using the midpoint rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  function f ( t )'
  write ( *, '(a)' ) '  evaluates the integrand.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Intervals   Estimate'
  write ( *, '(a)' ) ' '

  a = 0.0
  b = 1000.0
  int_num = 100

  do test = 1, 3

    call midpoint ( a, b, f, int_num, quad )

    write ( *, '(2x,i8,2x,g14.6)' ) int_num, quad

    int_num = int_num * 100

  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'F90_SIMPLE:'
  write ( *, '(a)' ) '  Normal end of execution.'
 
  stop
end
