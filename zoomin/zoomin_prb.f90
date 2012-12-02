program main

!*****************************************************************************80
!
!! MAIN is the main program for ZOOMIN_PRB.
!
!  Discussion:
!
!    ZOOMIN_PRB calls ZOOMIN to find a root of a nonlinear equation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZOOMIN_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ZOOMIN library.'

  call test01 ( )

  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZOOMIN_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 runs the tests on a polynomial function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), external :: func01
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) mult
  integer ( kind = 4 ) nder
  integer ( kind = 4 ) nsub

  a = 1.0D+00
  b = 1.5D+00
  c = 4.0D+00
!
!  Set the error tolerance.
!
  abserr = 0.00001D+00
!
!  Tell the program how many derivatives are available.
!
  nder = 3
!
!  IPOLY is the order of the polynomial function
!  or -1 if the function is not a polynomial.
!
  ipoly = 3
!
!  KMAX is the maximum number of iterations.
!
  kmax = 30
!
!  MULT is the multiplicity of the root.
!  If not known, set MULT to 1.
!
  mult = 1
!
!  For modified Newton methods, set the number of substeps.
!
  nsub = 3
!
!  Call ZOOMIN.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  (Polynomial function F(X))'
  write ( *, '(a)' ) '  Find a root of F(X)=(X+3)*(X+3)*(X-2)=0'
  write ( *, '(a)' ) ' '

  call zoomin ( a, b, c, abserr, kmax, func01, ipoly, mult, nder, nsub )

  return
end
function func01 ( x, ider )

!*****************************************************************************80
!
!! FUNC01 computes the function value for the first test.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the point at which the evaluation is to take place.
!
!    Input, integer ( kind = 4 ) IDER, specifies what is to be evaluated:
!    0, evaluate the function.
!    1, evaluate the first derivative.
!    2, evaluate the second derivative.
!    3, evaluate the third derivative.
!
!    Output, real FUNC01, the value of the function or derivative.
!
  implicit none

  real ( kind = 8 ) func01
  integer ( kind = 4 ) ider
  real ( kind = 8 ) x

  if ( ider == 0 ) then
    func01 = ( x + 3.0D+00 )**2 * ( x - 2.0E+00 )
  else if ( ider == 1 ) then
    func01 = ( x + 3.0D+00 ) * ( 3.0E+00 * x - 1.0E+00 )
  else if ( ider == 2 ) then
    func01 = 6.0D+00 * x + 8.0E+00
  else if ( ider == 3 ) then
    func01 = 6.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FUNC01 - Fatal error!'
    write ( *, '(a,i8)' ) '  Derivative of order IDER = ', ider
    write ( *, '(a)' ) '  was requested.'
    stop
  end if

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 runs the tests on a nonpolynomial function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), external :: func02
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) mult
  integer ( kind = 4 ) nder
  integer ( kind = 4 ) nsub

  a = 0.9D+00
  b = 0.4D+00
  c = 0.5D+00
!
!  Set the error tolerance.
!
  abserr = 0.00001D+00
!
!  Tell the program how many derivatives are available.
!
  nder = 3
!
!  IPOLY is the order of the polynomial function
!  or -1 if the function is not a polynomial.
!
  ipoly = - 1
!
!  KMAX is the maximum number of iterations.
!
  kmax = 60
!
!  MULT is the multiplicity of the root.
!  If not known, just set MULT to 1.
!
  mult = 1
!
!  For modified Newton methods, set the number of substeps.
!
  nsub = 3
!
!  Call ZOOMIN.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  (Nonpolynomial function F(X))'
  write ( *, '(a)' ) '  Find a root of F(X) = COS(X) - X'
  write ( *, '(a)' ) ' '

  call zoomin ( a, b, c, abserr, kmax, func02, ipoly, mult, nder, nsub )

  return
end
function func02 ( x, ider )

!*****************************************************************************80
!
!! FUNC02 computes the function value for the second test.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the point at which the evaluation is to take place.
!
!    Input, integer ( kind = 4 ) IDER, specifies what is to be evaluated:
!    0, evaluate the function.
!    1, evaluate the first derivative.
!    2, evaluate the second derivative.
!    3, evaluate the third derivative.
!
!    Output, real FUNC02, the value of the function or derivative.
!
  implicit none

  real ( kind = 8 ) func02
  integer ( kind = 4 ) ider
  real ( kind = 8 ) x

  if ( ider == 0 ) then
    func02 = cos ( x ) - x
  else if ( ider == 1 ) then
    func02 = - sin ( x ) - 1.0D+00
  else if ( ider == 2 ) then
    func02 = - cos ( x )
  else if ( ider == 3 ) then
    func02 = sin ( x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FUNC02 - Fatal error!'
    write ( *, '(a,i8)' ) '  Derivative of order IDER = ', ider
    write ( *, '(a)' ) '  was requested.'
    stop
  end if

  return
end
