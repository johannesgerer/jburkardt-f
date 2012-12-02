program main

!*****************************************************************************80
!
!! NL2SOL_PRB1 tests NL2SOL and NL2SNO on the Madsen example.
!
!  Modified:
!
!    03 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NL2SOL_PRB1:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the NL2SOL library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NL2SOL_PRB1:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests NL2SOL on the Madsen example.
!
  implicit none

  integer, parameter :: meqn = 3
  integer, parameter :: nvar = 2

  integer iv(60+nvar)
  external madj
  external madr
  external ufparm
  integer uiparm(1)
  real urparm(1)
  real v(93+meqn*nvar+3*meqn+nvar*(3*nvar+33)/2)
  real x(nvar)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Test the NL2SOL routine,'
  write ( *, '(a)' ) '  which requires a user residual and jacobian.'
!
!  Set the initial solution estimate.
!
  x(1:2) = (/ 3.0E+00, 1.0E+00 /)

  iv(1) = 0

  call nl2sol ( meqn, nvar, x, madr, madj, iv, v, uiparm, urparm, ufparm )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests NL2SNO on the Madsen example.
!
  implicit none

  integer, parameter :: meqn = 3
  integer, parameter :: nvar = 2

  integer iv(60+nvar)
  external madr
  external ufparm
  integer uiparm(1)
  real urparm(1)
  real v(147)
  real x(nvar)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Test the NL2SNO routine,'
  write ( *, '(a)' ) '  which requires only a user residual.'
  write ( *, '(a)' ) '  The jacobian is approximated internally.'
!
!  Set the initial solution estimate.
!
  x(1:2) = (/ 3.0E+00, 1.0E+00 /)

  iv(1) = 0

  call nl2sno ( meqn, nvar, x, madr, iv, v, uiparm, urparm, ufparm )

  return
end
subroutine madr ( meqn, nvar, x, nf, r, uiparm, urparm, ufparm )

!*****************************************************************************80
!
!! MADR computes the Madsen residual.
!
!  Discussion:
!
!    Given the value of the vector X, this routine computes the
!    value of F(X), the vector function whose norm we are trying
!    to minimize.
!
!  Modified:
!
!    25 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MEQN, the number of functions.
!
!    Input, integer NVAR, the number of variables.
!
!    Input, real X(NVAR), the current value of the variables.
!
!    Input, integer NF, the number of times the residual routine
!    has been called so far.
!
!    Output, real R(MEQN), the residual vector, that is, the
!    value of the functions for the given input value of the variables.
!
!    Input, integer UIPARM(*), a user array.
!
!    Input, real URPARM(*), a user array.
!
!    Input, external UFPARM, an external reference to a user subroutine
!    or function.
!
  implicit none

  integer meqn
  integer nvar

  integer nf
  real r(meqn)
  external ufparm
  integer uiparm(*)
  real urparm(*)
  real x(nvar)

  r(1) = x(1)**2 + x(2)**2 + x(1) * x(2)
  r(2) = sin ( x(1) )
  r(3) = cos ( x(2) )

  return
end
subroutine madj ( meqn, nvar, x, nf, jac, uiparm, urparm, ufparm )

!*****************************************************************************80
!
!! MADJ computes the Madsen jacobian.
!
!  Discussion:
!
!    Given the value of the vector X, this routine computes the
!    jacobian matrix of the vector function F(X), which has the
!    form
!
!      Jac(I,J) = d Fi / d Xj
!
!  Modified:
!
!    25 November 2002
!
!  Parameters:
!
!    Input, integer MEQN, the number of functions.
!
!    Input, integer NVAR, the number of variables.
!
!    Input, real X(NVAR), the current value of the variables.
!
!    Input, integer NF, the number of times the residual routine
!    has been called so far.
!
!    Output, real JAC(MEQN,NVAR), the jacobian matrix.  JAC(I,J) is
!    the derivative of function I with respect to variable J.
!
!    Input, integer UIPARM(*), a user array.
!
!    Input, real URPARM(*), a user array.
!
!    Input, external UFPARM, an external reference to a user subroutine
!    or function.
!
  implicit none

  integer meqn
  integer nvar

  real jac(meqn,nvar)
  integer nf
  external ufparm
  integer uiparm(*)
  real urparm(*)
  real x(nvar)

  jac(1,1) = 2.0E+00 * x(1) + x(2)
  jac(1,2) = 2.0E+00 * x(2) + x(1)

  jac(2,1) = cos ( x(1) )
  jac(2,2) = 0.0E+00

  jac(3,1) = 0.0E+00
  jac(3,2) = -sin ( x(2) )

  return
end
subroutine ufparm ( meqn, nvar, x )

!*****************************************************************************80
!
!! UFPARM is a user-supplied external routine.
!
!  Discussion:
!
!    The name of the routine, the argument list, and even whether
!    it is a function or subroutine, are left to the user.
!
!    NL2SOL simply passes the external reference from the calling
!    program through to the residual and jacobian routines.
!
!    If the user has no need for this facility, then a dummy
!    routine like this one may be used.
!
!  Modified:
!
!    07 February 2003
!
!  Parameters:
!
!    Input, integer MEQN, the number of functions.
!
!    Input, integer NVAR, the number of variables.
!
!    Input, real X(NVAR), the current value of the variables.
!
  implicit none

  integer meqn
  integer nvar

  real x(nvar)

  return
end

