program main

!*****************************************************************************80
!
!! MAIN is the main program for DQED_PRB1.
!
!  Discussion:
!
!    DQED_PRB1 tests DQED.
!
!    This routine illustrates the use of DQED, the Hanson-Krogh nonlinear least
!    squares solver, by solving a certain heart dipole moment equation.
!
!  Modified:
!
!    11 September 2002
!
!  Reference: 
!
!    John Dennis, David Gay, Phuong Vu.
!    A new nonlinear equations test problem, 
!    Rice University Department of Mathematics 
!    Report 83-16, 6/83, 6/85
!
  implicit none

  integer ( kind = 4 ), parameter :: ldfj = 8
  integer ( kind = 4 ), parameter :: liwa = 150
  integer ( kind = 4 ), parameter :: lwa = 2000
  integer ( kind = 4 ), parameter :: nvars = 8

  real ( kind = 8 ) bl(nvars+2)
  real ( kind = 8 ) bu(nvars+2)
  real ( kind = 8 ) df(nvars)
  external dqedhd
  real ( kind = 8 ) fj(ldfj,nvars+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(nvars+2)
  integer ( kind = 4 ) iopt(28)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iwa(liwa)
  integer ( kind = 4 ) mode
  real ( kind = 8 ) ropt(1)
  real ( kind = 8 ) sigma(nvars)
  character ( len = 20 ) title
  real ( kind = 8 ) wa(lwa)
!
!  To have the jacobian checked, set WANT_PR to true.
!
  logical, parameter :: want_pr = .false.
  real ( kind = 8 ) x(nvars)
  real ( kind = 8 ) xsave(nvars)
  real ( kind = 8 ) y(nvars)

  common /sigma/ sigma
  common /mode/  mode

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQED_PRB1'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A set of tests for DQED, which can solve'
  write ( *, '(a)' ) '  bounded and constrained linear least squares problems'
  write ( *, '(a)' ) '  and systems of nonlinear equations.'

  do
!
!  Read the problem title.
!
    read ( *, '(a)', iostat = ios ) title

    if ( ios /= 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
!
!  Read in the sums of values and initial parameter values.
!
    read ( *, *, iostat = ios ) sigma(1:nvars), xsave(1:nvars)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Unexpected end of file or other error.'
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'The input SIGMA vector:'
    write ( *, '(a)' ) ' '
    do i = 1, nvars
      write ( *, '(i6,g14.6)' ) i, sigma(i) 
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'The initial estimate for X:'
    write ( *, '(a)' ) ' '
    do i = 1, nvars
      write ( *, '(i6,g14.6)' ) i, xsave(i) 
    end do
!
!  Test the partial derivative computation.
!
    if ( want_pr ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Test the partial derivative computation:'
      write ( *, '(a)' ) ' '

      call dpchek ( df, dqedhd, fj, iopt, ldfj, nvars, ropt, xsave, y )

    end if
!
!  Test 1, with all equations normal.
!
    mode = 0

    x(1:nvars) = xsave(1:nvars)

    call test01 ( bl, bu, fj, ind, iopt, iwa, ldfj, liwa, lwa, mode, nvars, &
      ropt, sigma, wa, x )
!
!  Test 1, with first two linear equations used as constraints.
!
    mode = 1

    x(1:nvars) = xsave(1:nvars)

    call test01 ( bl, bu, fj, ind, iopt, iwa, ldfj, liwa, lwa, mode, nvars, &
      ropt, sigma, wa, x )
!
!  Test 2, with all equations normal.
!
    mode = 0

    x(1:nvars) = xsave(1:nvars)

    call test02 ( bl, bu, fj, ind, iopt, iwa, ldfj, liwa, lwa, mode, nvars, &
      ropt, sigma, wa, x )
!
!  Test 2, with first two linear equations used as constraints.
!
    mode = 1

    x(1:nvars) = xsave(1:nvars)

    call test02 ( bl, bu, fj, ind, iopt, iwa, ldfj, liwa, lwa, mode, nvars, &
      ropt, sigma, wa, x )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQED_PRB1'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( bl, bu, fj, ind, iopt, iwa, ldfj, liwa, lwa, mode, &
  nvars, ropt, sigma, wa, x )

!*****************************************************************************80
!
!! TEST01 uses an analytic jacobian.
!
!  Discussion:
!
!    The linear equations will not be constraints if MODE = 0.
!    Convergence is normally faster if the linear equations are used
!    as constraints, MODE = 1.
!
  implicit none

  integer ( kind = 4 ) ldfj
  integer ( kind = 4 ) liwa
  integer ( kind = 4 ) lwa
  integer ( kind = 4 ) nvars

  real ( kind = 8 ) bl(nvars+2)
  real ( kind = 8 ) bu(nvars+2)
  external dqedhd
  real ( kind = 8 ) fj(ldfj,nvars+1)
  real ( kind = 8 ) fnorm
  integer ( kind = 4 ) igo
  integer ( kind = 4 ) ind(nvars+2)
  integer ( kind = 4 ) iopt(28)
  integer ( kind = 4 ) iwa(liwa)
  integer ( kind = 4 ) mcon
  integer ( kind = 4 ) mequa
  integer ( kind = 4 ) mode
  real ( kind = 8 ) ropt(1)
  real ( kind = 8 ) sigma(nvars)
  real ( kind = 8 ) wa(lwa)
  real ( kind = 8 ) x(nvars)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use an analytic jacobian.'
  write ( *, '(a,i6)' ) '  MODE = ', mode
!
!  Tell how much storage the solver has.
!
  iwa(1) = lwa
  iwa(2) = liwa
!
!  Set up print option.  (not used here).
!
  iopt(1) = -1
  iopt(2) = 1
!
!  Set up for linear model without any quadratic terms. (not used).
!
  iopt(3) = -14
  iopt(4) = 1
!
!  Do not allow convergence to be claimed on small steps.
!
  iopt(5) = 17
  iopt(6) = 1
  iopt(7) = 15
!
!  Allow up to NVARS quadratic model terms.
!
  iopt(8) = nvars
!
!  Change condition number for quadratic model degree.
!
  iopt(9) = 10
  iopt(10) = 1
  ropt(1) = 1.0D+04
!
!  No more options.
!
  iopt(11) = 99
!
!  MODE = 0, no constraints.
!
  if ( mode == 0 ) then

    mcon = 0
!
!  MODE = 1, there are constraints.
!
!  (The first two equations are linear, and can be used as constraints).
!
  else

    mcon = 2
    ind(nvars+1) = 3
    ind(nvars+2) = 3
    bl(nvars+1) = sigma(1)
    bu(nvars+1) = sigma(1)
    bl(nvars+2) = sigma(2)
    bu(nvars+2) = sigma(2)

  end if
 
  mequa = nvars - mcon
!
!  All variables are otherwise free.
!
  ind(1:nvars) = 4

  call dqed ( dqedhd, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, fnorm, &
    igo, iopt, ropt, iwa, wa )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed minimizing X:'
  write ( *, '(a)' ) ' '
  write ( *, '(g14.6)' ) x(1:nvars)
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Residual after the fit = ',fnorm
  write ( *, '(a,i6)' ) '  DQED output flag IGO = ',igo
 
  return
end
subroutine test02 ( bl, bu, fj, ind, iopt, iwa, ldfj, liwa, lwa, mode, &
  nvars, ropt, sigma, wa, x )

!*****************************************************************************80
!
!! TEST02 uses a finite-difference approximated jacobian.
!
!  Discussion:
!
!    The linear equations will not be constraints if MODE = 0.
!    Convergence is normally faster if the linear equations are used
!    as constraints, MODE = 1.
!
  implicit none

  integer ( kind = 4 ) ldfj
  integer ( kind = 4 ) liwa
  integer ( kind = 4 ) lwa
  integer ( kind = 4 ) nvars

  real ( kind = 8 ) bl(nvars+2)
  real ( kind = 8 ) bu(nvars+2)
  real ( kind = 8 ) fj(ldfj,nvars+1)
  external fjaprx
  real ( kind = 8 ) fnorm
  integer ( kind = 4 ) igo
  integer ( kind = 4 ) ind(nvars+2)
  integer ( kind = 4 ) iopt(28)
  integer ( kind = 4 ) iwa(liwa)
  integer ( kind = 4 ) mcon
  integer ( kind = 4 ) mequa
  integer ( kind = 4 ) mode
  real ( kind = 8 ) ropt(1)
  real ( kind = 8 ) sigma(nvars)
  real ( kind = 8 ) wa(lwa)
  real ( kind = 8 ) x(nvars)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Use an approximate jacobian.'
  write ( *, '(a,i6)' ) '  MODE = ', mode
!
!  Tell how much storage the solver has.
!
  iwa(1) = lwa
  iwa(2) = liwa
!
!  Set up print option.  (not used here).
!
  iopt(1) = -1
  iopt(2) = 1
!
!  Set up for linear model without any quadratic terms. (not used).
!
  iopt(3) = -14
  iopt(4) = 1
!
!  Do not allow convergence to be claimed on small steps.
!
  iopt(5) = 17
  iopt(6) = 1
  iopt(7) = 15
!
!  Allow up to NVARS quadratic model terms.
!
  iopt(8) = nvars
!
!  Change condition number for quadratic model degree.
!
  iopt(9) = 10
  iopt(10) = 1
  ropt(1) = 10000.0D+00
!
!  No more options.
!
  iopt(11) = 99
!
!  MODE = 0, no constraints.
!
  if ( mode == 0 ) then

    mcon = 0
!
!  MODE = 1, there are constraints.
!  (The first two equations are linear, and can be used as
!  constraints instead).
!
  else

    mcon = 2
    ind(nvars+1) = 3
    ind(nvars+2) = 3
    bl(nvars+1) = sigma(1)
    bu(nvars+1) = sigma(1)
    bl(nvars+2) = sigma(2)
    bu(nvars+2) = sigma(2)

  end if
 
  mequa = nvars - mcon
!
!  All variables are otherwise free.
!
  ind(1:nvars) = 4

  call dqed ( fjaprx, mequa, nvars, mcon, ind, bl, bu, x, fj, ldfj, fnorm, &
    igo, iopt, ropt, iwa, wa )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Computed minimizing X:'
  write ( *, '(a)' ) ' '
  write ( *, '(g14.6)' ) x(1:nvars)
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Residual after the fit = ', fnorm
  write ( *, '(a,i6)' ) '  DQED output flag IGO = ', igo
 
  return
end
subroutine dqedhd ( x, fj, ldfj, igo, iopt, ropt )

!*****************************************************************************80
!
!! DQEDHD evaluates functions and derivatives for DQED.
!
!  Discussion:
!
!    The user problem has MCON constraint functions,
!    MEQUA least squares equations, and involves NVARS
!    unknown variables.
!
!    When this subprogram is entered, the general (near)
!    linear constraint partial derivatives, the derivatives
!    for the least squares equations, and the associated
!    function values are placed into the array FJ.
!    all partials and functions are evaluated at the point
!    in X.  then the subprogram returns to the calling
!    program unit. 
!
!    Typically one could do the following steps:
!
!    if ( igo /= 0 ) then
!      place the partials of the i-th constraint function with respect to 
!      variable j in the array fj(i,j), i = 1,...,mcon, j=1,...,nvars.
!
!    place the values of the i-th constraint equation into fj(i,nvars+1).
!
!    if ( igo /= 0 ) then
!      place the partials of the i-th least squares equation with respect 
!      to variable j in the array fj(i,j), i = 1,...,mequa, j = 1,...,nvars.
!
!    place the value of the i-th least squares equation into fj(i,nvars+1).
!
  implicit none

  integer ( kind = 4 ) ldfj
  integer ( kind = 4 ), parameter :: nvars = 8

  real ( kind = 8 ) fj(ldfj,nvars+1)
  integer ( kind = 4 ) igo
  integer ( kind = 4 ) iopt(*)
  integer ( kind = 4 ) mcon
  integer ( kind = 4 ) mequa
  integer ( kind = 4 ) mode
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) x(nvars)

  common /mode/  mode

  if ( mode == 0 ) then
    mcon = 0
  else
    mcon = 2
  end if

  mequa = nvars - mcon

  if ( igo /= 0 ) then
    call jack ( fj, ldfj, nvars, x )
  end if

  call func ( fj(1,9), iopt, mcon, mequa, nvars, ropt, x )

  return
end
subroutine fjaprx ( x, fj, ldfj, igo, iopt, ropt )

!*****************************************************************************80
!
!! FJAPRX uses DIFCEN to approximate the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldfj
  integer ( kind = 4 ), parameter :: nvars = 8

  real ( kind = 8 ) fj(ldfj,nvars+1)
  external func
  real ( kind = 8 ) fx(8)
  integer ( kind = 4 ) igo
  integer ( kind = 4 ) iopt(*)
  integer ( kind = 4 ) mcon
  integer ( kind = 4 ) mequa
  integer ( kind = 4 ) mode
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) x(nvars)

  common /mode/  mode

  if ( mode == 0 ) then
    mcon = 0
  else
    mcon = 2
  end if

  mequa = nvars - mcon

  call difcen ( fj, func, fx, iopt, ldfj, mcon, mequa, nvars, ropt, x )

  call func ( fj(1,nvars+1), iopt, mcon, mequa, nvars, ropt, x )

  return
end
subroutine func ( fx, iopt, mcon, mequa, nvars, ropt, x )

!*****************************************************************************80
!
!! FUNC evaluates the functions.
!
!  Discussion:
!
!    The system of equations has the form:
!
!    x(1)+x(2) = sigma(1)
!
!    x(3)+x(4) = sigma(2)
!
!    x(1)*x(5)+x(2)*x(6)-x(3)*x(7)-x(4)*x(8) = sigma(3)
!
!    x(1)*x(7)+x(2)*x(8)+x(3)*x(5)+x(4)*x(6) = sigma(4)
!
!    x(1)*(x(5)**2-x(7)**2)+x(2)*(x(6)**2-x(8)**2)+x(3)*(-2.0*x(5)*x(7))
!      +x(4)*(-2.0*x(6)*x(8)) = sigma(5)
!
!    x(1)*(2.0*x(5)*x(7))+x(2)*(2.0*x(6)*x(8))+x(3)*(x(5)**2-x(7)**2)
!      +x(4)*(x(6)**2-x(8)**2) = sigma(6)
!
!    x(1)*(x(5)*(x(5)**2-3.0*x(7)**2))+x(2)*(x(6)*(x(6)**2-3.0*x(8)**2))
!      +x(3)*(x(7)*(x(7)**2-3.0*x(5)**2))+x(4)*(x(8)*(x(8)**2-3.0*x(6)**2))
!       = sigma(7)
!
!    x(1)*(-x(7)*(x(7)**2-3.0*x(5)**2))+x(2)*(-x(8)*(x(8)**2-3.0*x(6)**2))
!      +x(3)*(x(5)*(x(5)**2-3.0*x(7)**2))+x(4)*(x(6)*(x(6)**2-3.0*x(8)**2))
!       = sigma(8)
!
  implicit none

  integer ( kind = 4 ) mcon
  integer ( kind = 4 ) mequa
  integer ( kind = 4 ) nvars

  real ( kind = 8 ) fx(mequa+mcon)
  integer ( kind = 4 ) iopt(*)
  integer ( kind = 4 ) mode
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) sigma(8)
  real ( kind = 8 ) x(nvars)

  common /sigma/ sigma
  common /mode/ mode
!
  if ( mode == 0 ) then
    fx(1) = x(1) + x(2) - sigma(1)
  else
    fx(1) = x(1) + x(2)
  end if

  if ( mode == 0 ) then
    fx(2) = x(3) + x(4) - sigma(2)
  else
    fx(2) = x(3) + x(4)
  end if

  fx(3) = x(1)*x(5) + x(2)*x(6) - x(3)*x(7) - x(4)*x(8) - sigma(3)

  fx(4) = x(1)*x(7) + x(2)*x(8) + x(3)*x(5) + x(4)*x(6) - sigma(4)

  fx(5) = x(1)*(x(5)**2-x(7)**2) + x(2)*(x(6)**2-x(8)**2) &
    +x(3)*(-2.0D+00*x(5)*x(7)) + x(4)*(-2.0D+00*x(6)*x(8)) - sigma(5)

  fx(6) = x(1)*(2.0D+00*x(5)*x(7)) + x(2)*(2.0D+00*x(6)*x(8)) &
    + x(3)*(x(5)**2-x(7)**2) &
    + x(4)*(x(6)**2-x(8)**2) - sigma(6)

  fx(7) = x(1)*(x(5)*(x(5)**2-3.0*x(7)**2)) &
    + x(2)*(x(6)*(x(6)**2-3.0D+00*x(8)**2)) &
    + x(3)*(x(7)*(x(7)**2-3.0D+00*x(5)**2)) &
    + x(4)*(x(8)*(x(8)**2-3.0D+00*x(6)**2)) &
    - sigma(7)

  fx(8) = x(1)*(-x(7)*(x(7)**2-3.0D+00*x(5)**2)) &
    + x(2)*(-x(8)*(x(8)**2-3.0D+00*x(6)**2)) &
    + x(3)*(x(5)*(x(5)**2-3.0D+00*x(7)**2)) &
    + x(4)*(x(6)*(x(6)**2-3.0D+00*x(8)**2)) - sigma(8)

  return
end
subroutine jack ( fj, ldfj, nvars, x )

!*****************************************************************************80
!
!! JACK evaluates the partial derivatives of the functions.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) FJ(LDFJ,NVARS), the MCON+MEQUA by NVARS
!    matrix of partial derivatives.
!
!    Input, integer ( kind = 4 ) LDFJ, the leading dimension of FJ, which must be
!    at least MCON+MEQUA.
!
!    Input, integer ( kind = 4 ) NVARS, the number of variables.
!
!    Input, real X(NVARS), the point at which the partial derivatives
!    are to be evaluated.
!
  implicit none

  integer ( kind = 4 ) ldfj
  integer ( kind = 4 ) nvars

  real ( kind = 8 ) fj(ldfj,nvars)
  real ( kind = 8 ) x(nvars)
!
!  Equation #1: 
!
!    x(1)+x(2) = sigma(1)
!
  fj(1,1) = 1.0D+00
  fj(1,2) = 1.0D+00
  fj(1,3) = 0.0D+00
  fj(1,4) = 0.0D+00
  fj(1,5) = 0.0D+00
  fj(1,6) = 0.0D+00
  fj(1,7) = 0.0D+00
  fj(1,8) = 0.0D+00
!
!  Equation #2:
!
!    x(3)+x(4) = sigma(2)
!
  fj(2,1) = 0.0D+00
  fj(2,2) = 0.0D+00
  fj(2,3) = 1.0D+00
  fj(2,4) = 1.0D+00
  fj(2,5) = 0.0D+00
  fj(2,6) = 0.0D+00
  fj(2,7) = 0.0D+00
  fj(2,8) = 0.0D+00
!
!  Equation #3:
!    x(1)*x(5)+x(2)*x(6)-x(3)*x(7)-x(4)*x(8) = sigma(3)
!
  fj(3,1) = x(5)
  fj(3,2) = x(6)
  fj(3,3) = -x(7)
  fj(3,4) = -x(8)
  fj(3,5) = x(1)
  fj(3,6) = x(2)
  fj(3,7) = -x(3)
  fj(3,8) = -x(4)
!
!  Equation #4:
!    x(1)*x(7)+x(2)*x(8)+x(3)*x(5)+x(4)*x(6) = sigma(4)
!
  fj(4,1) = x(7)
  fj(4,2) = x(8)
  fj(4,3) = x(5)
  fj(4,4) = x(6)
  fj(4,5) = x(3)
  fj(4,6) = x(4)
  fj(4,7) = x(1)
  fj(4,8) = x(2)
!
!  Equation #5:
!    x(1)*(x(5)**2-x(7)**2)+x(2)*(x(6)**2-x(8)**2)+x(3)*(-2.0*x(5)*x(7))
!      +x(4)*(-2.0*x(6)*x(8)) = sigma(5)
!
  fj(5,1) = x(5)**2-x(7)**2
  fj(5,2) = x(6)**2-x(8)**2
  fj(5,3) = -2.0D+00 *x(5)*x(7)
  fj(5,4) = -2.0D+00 *x(6)*x(8)
  fj(5,5) = 2.0D+00 *(x(1)*x(5)-x(3)*x(7))
  fj(5,6) = 2.0D+00 *(x(2)*x(6)-x(4)*x(8))
  fj(5,7) = -2.0D+00 *(x(1)*x(7)+x(3)*x(5))
  fj(5,8) = -2.0D+00 *(x(2)*x(8)+x(4)*x(6))
!
!  Equation #6:
!    x(1)*(2.0*x(5)*x(7))+x(2)*(2.0*x(6)*x(8))+x(3)*(x(5)**2-x(7)**2)
!      +x(4)*(x(6)**2-x(8)**2) = sigma(6)
!
  fj(6,1) = 2.0D+00 *x(5)*x(7)
  fj(6,2) = 2.0D+00 *x(6)*x(8)
  fj(6,3) = x(5)**2-x(7)**2
  fj(6,4) = x(6)**2-x(8)**2
  fj(6,5) = 2.0D+00 *(x(1)*x(7)+x(3)*x(5))
  fj(6,6) = 2.0D+00 *(x(2)*x(8)+x(4)*x(6))
  fj(6,7) = 2.0D+00 *(x(1)*x(5)-x(3)*x(7))
  fj(6,8) = 2.0D+00 *(x(2)*x(6)-x(4)*x(8))
!
!  Equation #7:
!    x(1)*(x(5)*(x(5)**2-3.0*x(7)**2))+x(2)*(x(6)*(x(6)**2-3.0*x(8)**2))
!      +x(3)*(x(7)*(x(7)**2-3.0*x(5)**2))+x(4)*(x(8)*(x(8)**2-3.0*x(6)**2))
!       = sigma(7)
!
  fj(7,1) = x(5)*(x(5)**2-3.0D+00*x(7)**2)
  fj(7,2) = x(6)*(x(6)**2-3.0D+00*x(8)**2)
  fj(7,3) = x(7)*(x(7)**2-3.0D+00*x(5)**2)
  fj(7,4) = x(8)*(x(8)**2-3.0D+00*x(6)**2)
  fj(7,5) = 3.0D+00*(x(1)*(x(5)**2-x(7)**2)+x(3)*(-2.0D+00*x(5)*x(7)))
  fj(7,6) = 3.0D+00*(x(2)*(x(6)**2-x(8)**2)+x(4)*(-2.0D+00*x(6)*x(8)))
  fj(7,7) = -3.0D+00*(x(1)*(2.0*x(5)*x(7))+x(3)*(x(5)**2-x(7)**2))
  fj(7,8) = -3.0D+00*(x(2)*(2.0*x(6)*x(8))+x(4)*(x(6)**2-x(8)**2))
!
!  Equation #8:
!    x(1)*(-x(7)*(x(7)**2-3.0*x(5)**2))+x(2)*(-x(8)*(x(8)**2-3.0*x(6)**2))
!      +x(3)*(x(5)*(x(5)**2-3.0*x(7)**2))+x(4)*(x(6)*(x(6)**2-3.0*x(8)**2))
!       = sigma(8)
!
  fj(8,1) = -x(7)*(x(7)**2-3.0D+00*x(5)**2)
  fj(8,2) = -x(8)*(x(8)**2-3.0D+00*x(6)**2)
  fj(8,3) = x(5)*(x(5)**2-3.0D+00*x(7)**2)
  fj(8,4) = x(6)*(x(6)**2-3.0D+00*x(8)**2)
  fj(8,5) = 3.0D+00*(x(1)*(2.0D+00*x(5)*x(7))+x(3)*(x(5)**2-x(7)**2))
  fj(8,6) = 3.0D+00*(x(2)*(2.0D+00*x(6)*x(8))+x(4)*(x(6)**2-x(8)**2))
  fj(8,7) = 3.0D+00*(x(1)*(x(5)**2-x(7)**2)+x(3)*(-2.0D+00*x(5)*x(7)))
  fj(8,8) = 3.0D+00*(x(2)*(x(6)**2-x(8)**2)+x(4)*(-2.0D+00*x(6)*x(8)))

  return
end
