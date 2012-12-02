subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function error_f ( x )

!*****************************************************************************80
!
!! ERROR_F evaluates the error function ERF(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    Rational Chebyshev approximations for the error function,
!    Mathematics of Computation,
!    1969, pages 631-638.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of ERF.
!
!    Output, real ( kind = 8 ) ERROR_F, the value of ERF(X).
!
  implicit none

  real    ( kind = 8 ), save, dimension ( 5 ) :: a = (/ &
    3.16112374387056560D+00, &
    1.13864154151050156D+02, &
    3.77485237685302021D+02, &
    3.20937758913846947D+03, &
    1.85777706184603153D-01 /)
  real    ( kind = 8 ), save, dimension ( 4 ) :: b = (/ &
    2.36012909523441209D+01, &
    2.44024637934444173D+02, &
    1.28261652607737228D+03, &
    2.84423683343917062D+03 /)
  real    ( kind = 8 ), save, dimension ( 9 ) :: c = (/ &
    5.64188496988670089D-01, &
    8.88314979438837594D+00, &
    6.61191906371416295D+01, &
    2.98635138197400131D+02, &
    8.81952221241769090D+02, &
    1.71204761263407058D+03, &
    2.05107837782607147D+03, &
    1.23033935479799725D+03, &
    2.15311535474403846D-08 /)
  real    ( kind = 8 ), save, dimension ( 8 ) :: d = (/ &
    1.57449261107098347D+01, &
    1.17693950891312499D+02, &
    5.37181101862009858D+02, &
    1.62138957456669019D+03, &
    3.29079923573345963D+03, &
    4.36261909014324716D+03, &
    3.43936767414372164D+03, &
    1.23033935480374942D+03 /)
  real    ( kind = 8 ) del
  real    ( kind = 8 ) error_f
  integer ( kind = 4 ) i
  real    ( kind = 8 ), save, dimension ( 6 ) :: p = (/ &
    3.05326634961232344D-01, &
    3.60344899949804439D-01, &
    1.25781726111229246D-01, &
    1.60837851487422766D-02, &
    6.58749161529837803D-04, &
    1.63153871373020978D-02 /)
  real    ( kind = 8 ), save, dimension ( 5 ) :: q = (/ &
    2.56852019228982242D+00, &
    1.87295284992346047D+00, &
    5.27905102951428412D-01, &
    6.05183413124413191D-02, &
    2.33520497626869185D-03 /)
  real    ( kind = 8 ), parameter :: sqrpi = 0.56418958354775628695D+00
  real    ( kind = 8 ), parameter :: thresh = 0.46875D+00
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xabs
  real    ( kind = 8 ), parameter :: xbig = 26.543D+00
  real    ( kind = 8 ) xden
  real    ( kind = 8 ) xnum
  real    ( kind = 8 ) xsq

  xabs = abs ( x )
!
!  Evaluate ERF(X) for |X| <= 0.46875.
!
  if ( xabs <= thresh ) then

    if ( epsilon ( xabs ) < xabs ) then
      xsq = xabs * xabs
    else
      xsq = 0.0D+00
    end if

    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do

    error_f = x * ( xnum + a(4) ) / ( xden + b(4) )
!
!  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
!
  else if ( xabs <= 4.0D+00 ) then

    xnum = c(9) * xabs
    xden = xabs
    do i = 1, 7
      xnum = ( xnum + c(i) ) * xabs
      xden = ( xden + d(i) ) * xabs
    end do

    error_f = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = aint ( xabs * 16.0D+00 ) / 16.0D+00
    del = ( xabs - xsq ) * ( xabs + xsq )
    error_f = exp ( - xsq * xsq ) * exp ( - del ) * error_f

    error_f = ( 0.5D+00 - error_f ) + 0.5D+00

    if ( x < 0.0D+00 ) then
      error_f = - error_f
    end if
!
!  Evaluate ERFC(X) for 4.0 < |X|.
!
  else

    if ( xbig <= xabs ) then

      if ( 0.0D+00 < x ) then
        error_f = 1.0D+00
      else
        error_f = -1.0D+00
      end if

    else

      xsq = 1.0D+00 / ( xabs * xabs )

      xnum = p(6) * xsq
      xden = xsq
      do i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
      end do

      error_f = xsq * ( xnum + p(5) ) / ( xden + q(5) )
      error_f = ( sqrpi - error_f ) / xabs
      xsq = aint ( xabs * 16.0D+00 ) / 16.0D+00
      del = ( xabs - xsq ) * ( xabs + xsq )
      error_f = exp ( - xsq * xsq ) * exp ( - del ) * error_f

      error_f = ( 0.5D+00 - error_f ) + 0.5D+00
      if ( x < 0.0D+00 ) then
        error_f = - error_f
      end if

    end if

  end if

  return
end
function p00_autonomous ( test )

!*****************************************************************************80
!
!! P00_AUTONOMOUS reports whether a given problem is autonomous.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the problem number.
!
!    Output, logical P00_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  integer ( kind = 4 ) test
  logical p00_autonomous
  logical p01_autonomous
  logical p02_autonomous
  logical p03_autonomous
  logical p04_autonomous
  logical p05_autonomous
  logical p06_autonomous
  logical p07_autonomous
  logical p08_autonomous
  logical p09_autonomous
  logical p10_autonomous
  logical p11_autonomous
  logical p12_autonomous
  logical p13_autonomous
  logical p14_autonomous
  logical p15_autonomous
  logical p16_autonomous
  logical p17_autonomous
  logical p18_autonomous
  logical p19_autonomous
  logical p20_autonomous
  logical p21_autonomous
  logical p22_autonomous
  logical p23_autonomous
  logical p24_autonomous
  logical p25_autonomous
  logical p26_autonomous
  logical p27_autonomous
  logical p28_autonomous
  logical p29_autonomous
  logical p30_autonomous
  logical p31_autonomous
  logical p32_autonomous
  logical p33_autonomous
  logical p34_autonomous
  logical p35_autonomous
  logical p36_autonomous
  logical p37_autonomous
  logical p38_autonomous
  logical p39_autonomous
  logical p40_autonomous

  if ( test == 1 ) then
    p00_autonomous = p01_autonomous ( )
  else if ( test == 2 ) then
    p00_autonomous = p02_autonomous ( )
  else if ( test == 3 ) then
    p00_autonomous = p03_autonomous ( )
  else if ( test == 4 ) then
    p00_autonomous = p04_autonomous ( )
  else if ( test == 5 ) then
    p00_autonomous = p05_autonomous ( )
  else if ( test == 6 ) then
    p00_autonomous = p06_autonomous ( )
  else if ( test == 7 ) then
    p00_autonomous = p07_autonomous ( )
  else if ( test == 8 ) then
    p00_autonomous = p08_autonomous ( )
  else if ( test == 9 ) then
    p00_autonomous = p09_autonomous ( )
  else if ( test == 10 ) then
    p00_autonomous = p10_autonomous ( )
  else if ( test == 11 ) then
    p00_autonomous = p11_autonomous ( )
  else if ( test == 12 ) then
    p00_autonomous = p12_autonomous ( )
  else if ( test == 13 ) then
    p00_autonomous = p13_autonomous ( )
  else if ( test == 14 ) then
    p00_autonomous = p14_autonomous ( )
  else if ( test == 15 ) then
    p00_autonomous = p15_autonomous ( )
  else if ( test == 16 ) then
    p00_autonomous = p16_autonomous ( )
  else if ( test == 17 ) then
    p00_autonomous = p17_autonomous ( )
  else if ( test == 18 ) then
    p00_autonomous = p18_autonomous ( )
  else if ( test == 19 ) then
    p00_autonomous = p19_autonomous ( )
  else if ( test == 20 ) then
    p00_autonomous = p20_autonomous ( )
  else if ( test == 21 ) then
    p00_autonomous = p21_autonomous ( )
  else if ( test == 22 ) then
    p00_autonomous = p22_autonomous ( )
  else if ( test == 23 ) then
    p00_autonomous = p23_autonomous ( )
  else if ( test == 24 ) then
    p00_autonomous = p24_autonomous ( )
  else if ( test == 25 ) then
    p00_autonomous = p25_autonomous ( )
  else if ( test == 26 ) then
    p00_autonomous = p26_autonomous ( )
  else if ( test == 27 ) then
    p00_autonomous = p27_autonomous ( )
  else if ( test == 28 ) then
    p00_autonomous = p28_autonomous ( )
  else if ( test == 29 ) then
    p00_autonomous = p29_autonomous ( )
  else if ( test == 30 ) then
    p00_autonomous = p30_autonomous ( )
  else if ( test == 31 ) then
    p00_autonomous = p31_autonomous ( )
  else if ( test == 32 ) then
    p00_autonomous = p32_autonomous ( )
  else if ( test == 33 ) then
    p00_autonomous = p33_autonomous ( )
  else if ( test == 34 ) then
    p00_autonomous = p34_autonomous ( )
  else if ( test == 35 ) then
    p00_autonomous = p35_autonomous ( )
  else if ( test == 36 ) then
    p00_autonomous = p36_autonomous ( )
  else if ( test == 37 ) then
    p00_autonomous = p37_autonomous ( )
  else if ( test == 38 ) then
    p00_autonomous = p38_autonomous ( )
  else if ( test == 39 ) then
    p00_autonomous = p39_autonomous ( )
  else if ( test == 40 ) then
    p00_autonomous = p40_autonomous ( )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_AUTONOMOUS - Fatal error!'
    write ( *, '(a,i8)' ) '  Unexpected value of TEST = ', test
    stop
  end if

  return
end
subroutine p00_equil ( test, neqn, y, next )

!*****************************************************************************80
!
!! P00_EQUIL returns equilibrium solutions of any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the problem number.
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  integer ( kind = 4 ) test
  real ( kind = 8 ) y(neqn)

  if ( test == 1 ) then
    call p01_equil ( neqn, y, next )
  else if ( test == 2 ) then
    call p02_equil ( neqn, y, next )
  else if ( test == 3 ) then
    call p03_equil ( neqn, y, next )
  else if ( test == 4 ) then
    call p04_equil ( neqn, y, next )
  else if ( test == 5 ) then
    call p05_equil ( neqn, y, next )
  else if ( test == 6 ) then
    call p06_equil ( neqn, y, next )
  else if ( test == 7 ) then
    call p07_equil ( neqn, y, next )
  else if ( test == 8 ) then
    call p08_equil ( neqn, y, next )
  else if ( test == 9 ) then
    call p09_equil ( neqn, y, next )
  else if ( test == 10 ) then
    call p10_equil ( neqn, y, next )
  else if ( test == 11 ) then
    call p11_equil ( neqn, y, next )
  else if ( test == 12 ) then
    call p12_equil ( neqn, y, next )
  else if ( test == 13 ) then
    call p13_equil ( neqn, y, next )
  else if ( test == 14 ) then
    call p14_equil ( neqn, y, next )
  else if ( test == 15 ) then
    call p15_equil ( neqn, y, next )
  else if ( test == 16 ) then
    call p16_equil ( neqn, y, next )
  else if ( test == 17 ) then
    call p17_equil ( neqn, y, next )
  else if ( test == 18 ) then
    call p18_equil ( neqn, y, next )
  else if ( test == 19 ) then
    call p19_equil ( neqn, y, next )
  else if ( test == 20 ) then
    call p20_equil ( neqn, y, next )
  else if ( test == 21 ) then
    call p21_equil ( neqn, y, next )
  else if ( test == 22 ) then
    call p22_equil ( neqn, y, next )
  else if ( test == 23 ) then
    call p23_equil ( neqn, y, next )
  else if ( test == 24 ) then
    call p24_equil ( neqn, y, next )
  else if ( test == 25 ) then
    call p25_equil ( neqn, y, next )
  else if ( test == 26 ) then
    call p26_equil ( neqn, y, next )
  else if ( test == 27 ) then
    call p27_equil ( neqn, y, next )
  else if ( test == 28 ) then
    call p28_equil ( neqn, y, next )
  else if ( test == 29 ) then
    call p29_equil ( neqn, y, next )
  else if ( test == 30 ) then
    call p30_equil ( neqn, y, next )
  else if ( test == 31 ) then
    call p31_equil ( neqn, y, next )
  else if ( test == 32 ) then
    call p32_equil ( neqn, y, next )
  else if ( test == 33 ) then
    call p33_equil ( neqn, y, next )
  else if ( test == 34 ) then
    call p34_equil ( neqn, y, next )
  else if ( test == 35 ) then
    call p35_equil ( neqn, y, next )
  else if ( test == 36 ) then
    call p36_equil ( neqn, y, next )
  else if ( test == 37 ) then
    call p37_equil ( neqn, y, next )
  else if ( test == 38 ) then
    call p38_equil ( neqn, y, next )
  else if ( test == 39 ) then
    call p39_equil ( neqn, y, next )
  else if ( test == 40 ) then
    call p40_equil ( neqn, y, next )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_EQUIL - Fatal error!'
    write ( *, '(a,i6)' ) '  Unrecognized problem number = ', test
    stop
  end if

  return
end
subroutine p00_euler_step ( test, neqn, t0, y0, t1, y1 )

!*****************************************************************************80
!
!! P00_EULER_STEP takes a single Euler step from (T0,Y0) to (T1,Y1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the problem number.
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T0, Y0(NEQN), the arguments of the derivative
!    function.
!
!    Input, real ( kind = 8 ) T1, the point at which an estimate of the solution
!    is desired.
!
!    Output, real ( kind = 8 ) Y1(NEQN), the estimated solution at T1.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) dt
  integer ( kind = 4 ) test
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  real ( kind = 8 ) y0(neqn)
  real ( kind = 8 ) y1(neqn)
  real ( kind = 8 ) yp0(neqn)

  dt = t1 - t0

  call p00_fun ( test, neqn, t0, y0, yp0 )

  y1(1:neqn) = y0(1:neqn) + dt * yp0(1:neqn)

  return
end
subroutine p00_fun ( test, neqn, t, y, yp )

!*****************************************************************************80
!
!! P00_FUN evaluates the derivative function for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the problem number.
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) test
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  if ( test == 1 ) then
    call p01_fun ( neqn, t, y, yp )
  else if ( test == 2 ) then
    call p02_fun ( neqn, t, y, yp )
  else if ( test == 3 ) then
    call p03_fun ( neqn, t, y, yp )
  else if ( test == 4 ) then
    call p04_fun ( neqn, t, y, yp )
  else if ( test == 5 ) then
    call p05_fun ( neqn, t, y, yp )
  else if ( test == 6 ) then
    call p06_fun ( neqn, t, y, yp )
  else if ( test == 7 ) then
    call p07_fun ( neqn, t, y, yp )
  else if ( test == 8 ) then
    call p08_fun ( neqn, t, y, yp )
  else if ( test == 9 ) then
    call p09_fun ( neqn, t, y, yp )
  else if ( test == 10 ) then
    call p10_fun ( neqn, t, y, yp )
  else if ( test == 11 ) then
    call p11_fun ( neqn, t, y, yp )
  else if ( test == 12 ) then
    call p12_fun ( neqn, t, y, yp )
  else if ( test == 13 ) then
    call p13_fun ( neqn, t, y, yp )
  else if ( test == 14 ) then
    call p14_fun ( neqn, t, y, yp )
  else if ( test == 15 ) then
    call p15_fun ( neqn, t, y, yp )
  else if ( test == 16 ) then
    call p16_fun ( neqn, t, y, yp )
  else if ( test == 17 ) then
    call p17_fun ( neqn, t, y, yp )
  else if ( test == 18 ) then
    call p18_fun ( neqn, t, y, yp )
  else if ( test == 19 ) then
    call p19_fun ( neqn, t, y, yp )
  else if ( test == 20 ) then
    call p20_fun ( neqn, t, y, yp )
  else if ( test == 21 ) then
    call p21_fun ( neqn, t, y, yp )
  else if ( test == 22 ) then
    call p22_fun ( neqn, t, y, yp )
  else if ( test == 23 ) then
    call p23_fun ( neqn, t, y, yp )
  else if ( test == 24 ) then
    call p24_fun ( neqn, t, y, yp )
  else if ( test == 25 ) then
    call p25_fun ( neqn, t, y, yp )
  else if ( test == 26 ) then
    call p26_fun ( neqn, t, y, yp )
  else if ( test == 27 ) then
    call p27_fun ( neqn, t, y, yp )
  else if ( test == 28 ) then
    call p28_fun ( neqn, t, y, yp )
  else if ( test == 29 ) then
    call p29_fun ( neqn, t, y, yp )
  else if ( test == 30 ) then
    call p30_fun ( neqn, t, y, yp )
  else if ( test == 31 ) then
    call p31_fun ( neqn, t, y, yp )
  else if ( test == 32 ) then
    call p32_fun ( neqn, t, y, yp )
  else if ( test == 33 ) then
    call p33_fun ( neqn, t, y, yp )
  else if ( test == 34 ) then
    call p34_fun ( neqn, t, y, yp )
  else if ( test == 35 ) then
    call p35_fun ( neqn, t, y, yp )
  else if ( test == 36 ) then
    call p36_fun ( neqn, t, y, yp )
  else if ( test == 37 ) then
    call p37_fun ( neqn, t, y, yp )
  else if ( test == 38 ) then
    call p38_fun ( neqn, t, y, yp )
  else if ( test == 39 ) then
    call p39_fun ( neqn, t, y, yp )
  else if ( test == 40 ) then
    call p40_fun ( neqn, t, y, yp )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FUN - Fatal error!'
    write ( *, '(a,i6)' ) '  Unrecognized problem number = ', test
    stop
  end if

  return
end
subroutine p00_jac ( test, neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P00_JAC evaluates the jacobian for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the problem number.
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  integer ( kind = 4 ) test
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  if ( test == 1 ) then
    call p01_jac ( neqn, ldj, t, y, jac )
  else if ( test == 2 ) then
    call p02_jac ( neqn, ldj, t, y, jac )
  else if ( test == 3 ) then
    call p03_jac ( neqn, ldj, t, y, jac )
  else if ( test == 4 ) then
    call p04_jac ( neqn, ldj, t, y, jac )
  else if ( test == 5 ) then
    call p05_jac ( neqn, ldj, t, y, jac )
  else if ( test == 6 ) then
    call p06_jac ( neqn, ldj, t, y, jac )
  else if ( test == 7 ) then
    call p07_jac ( neqn, ldj, t, y, jac )
  else if ( test == 8 ) then
    call p08_jac ( neqn, ldj, t, y, jac )
  else if ( test == 9 ) then
    call p09_jac ( neqn, ldj, t, y, jac )
  else if ( test == 10 ) then
    call p10_jac ( neqn, ldj, t, y, jac )
  else if ( test == 11 ) then
    call p11_jac ( neqn, ldj, t, y, jac )
  else if ( test == 12 ) then
    call p12_jac ( neqn, ldj, t, y, jac )
  else if ( test == 13 ) then
    call p13_jac ( neqn, ldj, t, y, jac )
  else if ( test == 14 ) then
    call p14_jac ( neqn, ldj, t, y, jac )
  else if ( test == 15 ) then
    call p15_jac ( neqn, ldj, t, y, jac )
  else if ( test == 16 ) then
    call p16_jac ( neqn, ldj, t, y, jac )
  else if ( test == 17 ) then
    call p17_jac ( neqn, ldj, t, y, jac )
  else if ( test == 18 ) then
    call p18_jac ( neqn, ldj, t, y, jac )
  else if ( test == 19 ) then
    call p19_jac ( neqn, ldj, t, y, jac )
  else if ( test == 20 ) then
    call p20_jac ( neqn, ldj, t, y, jac )
  else if ( test == 21 ) then
    call p21_jac ( neqn, ldj, t, y, jac )
  else if ( test == 22 ) then
    call p22_jac ( neqn, ldj, t, y, jac )
  else if ( test == 23 ) then
    call p23_jac ( neqn, ldj, t, y, jac )
  else if ( test == 24 ) then
    call p24_jac ( neqn, ldj, t, y, jac )
  else if ( test == 25 ) then
    call p25_jac ( neqn, ldj, t, y, jac )
  else if ( test == 26 ) then
    call p26_jac ( neqn, ldj, t, y, jac )
  else if ( test == 27 ) then
    call p27_jac ( neqn, ldj, t, y, jac )
  else if ( test == 28 ) then
    call p28_jac ( neqn, ldj, t, y, jac )
  else if ( test == 29 ) then
    call p29_jac ( neqn, ldj, t, y, jac )
  else if ( test == 30 ) then
    call p30_jac ( neqn, ldj, t, y, jac )
  else if ( test == 31 ) then
    call p31_jac ( neqn, ldj, t, y, jac )
  else if ( test == 32 ) then
    call p32_jac ( neqn, ldj, t, y, jac )
  else if ( test == 33 ) then
    call p33_jac ( neqn, ldj, t, y, jac )
  else if ( test == 34 ) then
    call p34_jac ( neqn, ldj, t, y, jac )
  else if ( test == 35 ) then
    call p35_jac ( neqn, ldj, t, y, jac )
  else if ( test == 36 ) then
    call p36_jac ( neqn, ldj, t, y, jac )
  else if ( test == 37 ) then
    call p37_jac ( neqn, ldj, t, y, jac )
  else if ( test == 38 ) then
    call p38_jac ( neqn, ldj, t, y, jac )
  else if ( test == 39 ) then
    call p39_jac ( neqn, ldj, t, y, jac )
  else if ( test == 40 ) then
    call p40_jac ( neqn, ldj, t, y, jac )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_JAC - Fatal error!'
    write ( *, '(a,i6)' ) '  Unrecognized problem number = ', test
    stop
  end if

  return
end
subroutine p00_neqn ( test, neqn )

!*****************************************************************************80
!
!! P00_NEQN returns the number of equations for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the problem number.
!
!    Output, integer ( kind = 4 ) NEQN, the number of variables.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) test

  if (  test == 1 ) then
    call p01_neqn ( neqn )
  else if ( test == 2 ) then
    call p02_neqn ( neqn )
  else if ( test == 3 ) then
    call p03_neqn ( neqn )
  else if ( test == 4 ) then
    call p04_neqn ( neqn )
  else if ( test == 5 ) then
    call p05_neqn ( neqn )
  else if ( test == 6 ) then
    call p06_neqn ( neqn )
  else if ( test == 7 ) then
    call p07_neqn ( neqn )
  else if ( test == 8 ) then
    call p08_neqn ( neqn )
  else if ( test == 9 ) then
    call p09_neqn ( neqn )
  else if ( test == 10 ) then
    call p10_neqn ( neqn )
  else if ( test == 11 ) then
    call p11_neqn ( neqn )
  else if ( test == 12 ) then
    call p12_neqn ( neqn )
  else if ( test == 13 ) then
    call p13_neqn ( neqn )
  else if ( test == 14 ) then
    call p14_neqn ( neqn )
  else if ( test == 15 ) then
    call p15_neqn ( neqn )
  else if ( test == 16 ) then
    call p16_neqn ( neqn )
  else if ( test == 17 ) then
    call p17_neqn ( neqn )
  else if ( test == 18 ) then
    call p18_neqn ( neqn )
  else if ( test == 19 ) then
    call p19_neqn ( neqn )
  else if ( test == 20 ) then
    call p20_neqn ( neqn )
  else if ( test == 21 ) then
    call p21_neqn ( neqn )
  else if ( test == 22 ) then
    call p22_neqn ( neqn )
  else if ( test == 23 ) then
    call p23_neqn ( neqn )
  else if ( test == 24 ) then
    call p24_neqn ( neqn )
  else if ( test == 25 ) then
    call p25_neqn ( neqn )
  else if ( test == 26 ) then
    call p26_neqn ( neqn )
  else if ( test == 27 ) then
    call p27_neqn ( neqn )
  else if ( test == 28 ) then
    call p28_neqn ( neqn )
  else if ( test == 29 ) then
    call p29_neqn ( neqn )
  else if ( test == 30 ) then
    call p30_neqn ( neqn )
  else if ( test == 31 ) then
    call p31_neqn ( neqn )
  else if ( test == 32 ) then
    call p32_neqn ( neqn )
  else if ( test == 33 ) then
    call p33_neqn ( neqn )
  else if ( test == 34 ) then
    call p34_neqn ( neqn )
  else if ( test == 35 ) then
    call p35_neqn ( neqn )
  else if ( test == 36 ) then
    call p36_neqn ( neqn )
  else if ( test == 37 ) then
    call p37_neqn ( neqn )
  else if ( test == 38 ) then
    call p38_neqn ( neqn )
  else if ( test == 39 ) then
    call p39_neqn ( neqn )
  else if ( test == 40 ) then
    call p40_neqn ( neqn )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_NEQN - Fatal error!'
    write ( *, '(a,i6)' ) '  Unrecognized problem number  = ', test
    stop
  end if

  return
end
subroutine p00_rk_step ( test, neqn, order, t0, y0, t1, y1 )

!*****************************************************************************80
!
!! P00_RK_STEP takes a single Runge-Kutta step from (T0,Y0) to (T1,Y1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the problem number.
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the Runge-Kutta method to
!    be employed.  Legal values are 1 through 5.
!
!    Input, real ( kind = 8 ) T0, Y0(NEQN), the arguments of the derivative
!    function.
!
!    Input, real ( kind = 8 ) T1, the point at which an estimate of the solution
!    is desired.
!
!    Output, real ( kind = 8 ) Y1(NEQN), the estimated solution at T1.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) dt
  integer ( kind = 4 ) test
  integer ( kind = 4 ) order
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  real ( kind = 8 ) tk1
  real ( kind = 8 ) tk2
  real ( kind = 8 ) tk3
  real ( kind = 8 ) tk4
  real ( kind = 8 ) tk5
  real ( kind = 8 ) y0(neqn)
  real ( kind = 8 ) y1(neqn)
  real ( kind = 8 ) yk1(neqn)
  real ( kind = 8 ) yk2(neqn)
  real ( kind = 8 ) yk3(neqn)
  real ( kind = 8 ) yk4(neqn)
  real ( kind = 8 ) yk5(neqn)
  real ( kind = 8 ) yp0(neqn)
  real ( kind = 8 ) ypk1(neqn)
  real ( kind = 8 ) ypk2(neqn)
  real ( kind = 8 ) ypk3(neqn)
  real ( kind = 8 ) ypk4(neqn)
  real ( kind = 8 ) ypk5(neqn)

  dt = t1 - t0

  if ( order == 1 ) then

    call p00_fun ( test, neqn, t0, y0, yp0 )

    y1(1:neqn) = y0(1:neqn) + dt * yp0(1:neqn)

  else if ( order == 2 ) then

    call p00_fun ( test, neqn, t0, y0, yp0 )

    yk1(1:neqn) = y0(1:neqn) + dt * yp0(1:neqn)

    tk1 = t0 + dt

    call p00_fun ( test, neqn, tk1, yk1, ypk1 )

    y1(1:neqn) = y0(1:neqn) + dt * ( yp0(1:neqn) + ypk1(1:neqn) ) / 2.0D+00

  else if ( order == 3 ) then

    call p00_fun ( test, neqn, t0, y0, yp0 )

    yk1(1:neqn) = y0(1:neqn) + 0.5D+00 * dt * yp0(1:neqn)

    tk1 = t0 + 0.5D+00 * dt
    call p00_fun ( test, neqn, tk1, yk1, ypk1 )

    yk2(1:neqn) = y0(1:neqn) + dt * ( 2.0D+00 * ypk1(1:neqn) - yp0(1:neqn) )

    tk2 = t0 + dt
    call p00_fun ( test, neqn, tk2, yk2, ypk2 )

    y1(1:neqn) = y0(1:neqn) + ( dt / 6.0D+00 ) &
      * ( yp0(1:neqn) + 4.0D+00 * ypk1(1:neqn) + ypk2(1:neqn) )

  else if ( order == 4 ) then

    call p00_fun ( test, neqn, t0, y0, yp0 )

    yk1(1:neqn) = y0(1:neqn) + 0.5D+00 * dt * yp0(1:neqn)

    tk1 = t0 + 0.5D+00 * dt

    call p00_fun ( test, neqn, tk1, yk1, ypk1 )

    yk2(1:neqn) = y0(1:neqn) + 0.5D+00 * dt * ypk1(1:neqn)

    tk2 = t0 + 0.5D+00 * dt
    call p00_fun ( test, neqn, tk2, yk2, ypk2 )

    tk3 = t0 + dt

    yk3(1:neqn) = y0(1:neqn) + dt * ypk2(1:neqn)

    call p00_fun ( test, neqn, tk3, yk3, ypk3 )

    y1(1:neqn) = y0(1:neqn) + ( dt / 6.0D+00 ) * ( &
                  yp0(1:neqn) &
      + 2.0D+00 * ypk1(1:neqn) &
      + 2.0D+00 * ypk2(1:neqn) &
      +           ypk3(1:neqn) )

  else if ( order == 5 ) then

    call p00_fun ( test, neqn, t0, y0, yp0 )

    yk1(1:neqn) = y0(1:neqn) + 0.25D+00 * dt * yp0(1:neqn)

    tk1 = t0 + 0.25D+00 * dt

    call p00_fun ( test, neqn, tk1, yk1, ypk1 )

    yk2(1:neqn) = y0(1:neqn) + dt * ( &
        3.0D+00 * yp0(1:neqn) &
      + 9.0D+00 * ypk1(1:neqn) ) / 32.0D+00

    tk2 = t0 + 3.0D+00 * dt / 8.0D+00

    call p00_fun ( test, neqn, tk2, yk2, ypk2 )

    yk3(1:neqn) = y0(1:neqn) + dt * ( &
        1932.0D+00 * yp0(1:neqn) &
      - 7200.0D+00 * ypk1(1:neqn) &
      + 7296.0D+00 * ypk2(1:neqn) ) / 2197.0D+00

    tk3 = t0 + 12.0D+00 * dt / 13.0D+00

    call p00_fun ( test, neqn, tk3, yk3, ypk3 )

    yk4(1:neqn) = y0(1:neqn) + dt * ( &
      + (  439.0D+00 /  216.0D+00 ) * yp0(1:neqn) &
      -      8.0D+00                * ypk1(1:neqn) &
      + ( 3680.0D+00 /  513.0D+00 ) * ypk2(1:neqn) &
      - (  845.0D+00 / 4104.0D+00 ) * ypk3(1:neqn) )

    tk4 = t0 + dt

    call p00_fun ( test, neqn, tk4, yk4, ypk4 )

    yk5(1:neqn) = y0(1:neqn) + dt * ( &
      - (    8.0D+00 /   27.0D+00 ) * yp0(1:neqn) &
      + (    2.0D+00              ) * ypk1(1:neqn) &
      - ( 3544.0D+00 / 2565.0D+00 ) * ypk2(1:neqn) &
      + ( 1859.0D+00 / 4104.0D+00 ) * ypk3(1:neqn) &
      - (   11.0D+00 /   40.0D+00 ) * ypk4(1:neqn) )

    tk5 = t0 + 0.5D+00 * dt

    call p00_fun ( test, neqn, tk5, yk5, ypk5 )

    y1(1:neqn)  = y0(1:neqn) + dt * ( &
        (    16.0D+00 / 135.0D+00   ) * yp0(1:neqn) &
      + (  6656.0D+00 / 12825.0D+00 ) * ypk2(1:neqn) &
      + ( 28561.0D+00 / 56430.0D+00 ) * ypk3(1:neqn) &
      - (     9.0D+00 / 50.0D+00    ) * ypk4(1:neqn) &
      + (     2.0D+00 / 55.0D+00    ) * ypk5(1:neqn) )

  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_RK_STEP - Fatal error!'
    write ( *, '(a,i6)' ) '  Unavailable Runge Kutta order = ', order
    stop
  end if

  return
end
subroutine p00_scale ( test, neqn, scale )

!*****************************************************************************80
!
!! P00_SCALE returns scale factors for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the problem number.
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) test
  real ( kind = 8 ) scale(neqn)

  if ( test == 1 ) then
    call p01_scale ( neqn, scale )
  else if ( test == 2 ) then
    call p02_scale ( neqn, scale )
  else if ( test == 3 ) then
    call p03_scale ( neqn, scale )
  else if ( test == 4 ) then
    call p04_scale ( neqn, scale )
  else if ( test == 5 ) then
    call p05_scale ( neqn, scale )
  else if ( test == 6 ) then
    call p06_scale ( neqn, scale )
  else if ( test == 7 ) then
    call p07_scale ( neqn, scale )
  else if ( test == 8 ) then
    call p08_scale ( neqn, scale )
  else if ( test == 9 ) then
    call p09_scale ( neqn, scale )
  else if ( test == 10 ) then
    call p10_scale ( neqn, scale )
  else if ( test == 11 ) then
    call p11_scale ( neqn, scale )
  else if ( test == 12 ) then
    call p12_scale ( neqn, scale )
  else if ( test == 13 ) then
    call p13_scale ( neqn, scale )
  else if ( test == 14 ) then
    call p14_scale ( neqn, scale )
  else if ( test == 15 ) then
    call p15_scale ( neqn, scale )
  else if ( test == 16 ) then
    call p16_scale ( neqn, scale )
  else if ( test == 17 ) then
    call p17_scale ( neqn, scale )
  else if ( test == 18 ) then
    call p18_scale ( neqn, scale )
  else if ( test == 19 ) then
    call p19_scale ( neqn, scale )
  else if ( test == 20 ) then
    call p20_scale ( neqn, scale )
  else if ( test == 21 ) then
    call p21_scale ( neqn, scale )
  else if ( test == 22 ) then
    call p22_scale ( neqn, scale )
  else if ( test == 23 ) then
    call p23_scale ( neqn, scale )
  else if ( test == 24 ) then
    call p24_scale ( neqn, scale )
  else if ( test == 25 ) then
    call p25_scale ( neqn, scale )
  else if ( test == 26 ) then
    call p26_scale ( neqn, scale )
  else if ( test == 27 ) then
    call p27_scale ( neqn, scale )
  else if ( test == 28 ) then
    call p28_scale ( neqn, scale )
  else if ( test == 29 ) then
    call p29_scale ( neqn, scale )
  else if ( test == 30 ) then
    call p30_scale ( neqn, scale )
  else if ( test == 31 ) then
    call p31_scale ( neqn, scale )
  else if ( test == 32 ) then
    call p32_scale ( neqn, scale )
  else if ( test == 33 ) then
    call p33_scale ( neqn, scale )
  else if ( test == 34 ) then
    call p34_scale ( neqn, scale )
  else if ( test == 35 ) then
    call p35_scale ( neqn, scale )
  else if ( test == 36 ) then
    call p36_scale ( neqn, scale )
  else if ( test == 37 ) then
    call p37_scale ( neqn, scale )
  else if ( test == 38 ) then
    call p38_scale ( neqn, scale )
  else if ( test == 39 ) then
    call p39_scale ( neqn, scale )
  else if ( test == 40 ) then
    call p40_scale ( neqn, scale )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SCALE - Fatal error!'
    write ( *, '(a,i6)' ) '  Unrecognized problem = ', test
  end if

  return
end
subroutine p00_start ( test, neqn, t_start, y_start )

!*****************************************************************************80
!
!! P00_START returns the starting point for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the problem number.
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) test
  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  if ( test == 1 ) then
    call p01_start ( neqn, t_start, y_start )
  else if ( test == 2 ) then
    call p02_start ( neqn, t_start, y_start )
  else if ( test == 3 ) then
    call p03_start ( neqn, t_start, y_start )
  else if ( test == 4 ) then
    call p04_start ( neqn, t_start, y_start )
  else if ( test == 5 ) then
    call p05_start ( neqn, t_start, y_start )
  else if ( test == 6 ) then
    call p06_start ( neqn, t_start, y_start )
  else if ( test == 7 ) then
    call p07_start ( neqn, t_start, y_start )
  else if ( test == 8 ) then
    call p08_start ( neqn, t_start, y_start )
  else if ( test == 9 ) then
    call p09_start ( neqn, t_start, y_start )
  else if ( test == 10 ) then
    call p10_start ( neqn, t_start, y_start )
  else if ( test == 11 ) then
    call p11_start ( neqn, t_start, y_start )
  else if ( test == 12 ) then
    call p12_start ( neqn, t_start, y_start )
  else if ( test == 13 ) then
    call p13_start ( neqn, t_start, y_start )
  else if ( test == 14 ) then
    call p14_start ( neqn, t_start, y_start )
  else if ( test == 15 ) then
    call p15_start ( neqn, t_start, y_start )
  else if ( test == 16 ) then
    call p16_start ( neqn, t_start, y_start )
  else if ( test == 17 ) then
    call p17_start ( neqn, t_start, y_start )
  else if ( test == 18 ) then
    call p18_start ( neqn, t_start, y_start )
  else if ( test == 19 ) then
    call p19_start ( neqn, t_start, y_start )
  else if ( test == 20 ) then
    call p20_start ( neqn, t_start, y_start )
  else if ( test == 21 ) then
    call p21_start ( neqn, t_start, y_start )
  else if ( test == 22 ) then
    call p22_start ( neqn, t_start, y_start )
  else if ( test == 23 ) then
    call p23_start ( neqn, t_start, y_start )
  else if ( test == 24 ) then
    call p24_start ( neqn, t_start, y_start )
  else if ( test == 25 ) then
    call p25_start ( neqn, t_start, y_start )
  else if ( test == 26 ) then
    call p26_start ( neqn, t_start, y_start )
  else if ( test == 27 ) then
    call p27_start ( neqn, t_start, y_start )
  else if ( test == 28 ) then
    call p28_start ( neqn, t_start, y_start )
  else if ( test == 29 ) then
    call p29_start ( neqn, t_start, y_start )
  else if ( test == 30 ) then
    call p30_start ( neqn, t_start, y_start )
  else if ( test == 31 ) then
    call p31_start ( neqn, t_start, y_start )
  else if ( test == 32 ) then
    call p32_start ( neqn, t_start, y_start )
  else if ( test == 33 ) then
    call p33_start ( neqn, t_start, y_start )
  else if ( test == 34 ) then
    call p34_start ( neqn, t_start, y_start )
  else if ( test == 35 ) then
    call p35_start ( neqn, t_start, y_start )
  else if ( test == 36 ) then
    call p36_start ( neqn, t_start, y_start )
  else if ( test == 37 ) then
    call p37_start ( neqn, t_start, y_start )
  else if ( test == 38 ) then
    call p38_start ( neqn, t_start, y_start )
  else if ( test == 39 ) then
    call p39_start ( neqn, t_start, y_start )
  else if ( test == 40 ) then
    call p40_start ( neqn, t_start, y_start )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_START - Fatal error!'
    write ( *, '(a,i6)' ) '  Unrecognized problem = ', test
  end if

  return
end
subroutine p00_stop ( test, neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P00_STOP returns the stopping point for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the problem number.
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) test
  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  if ( test == 1 ) then
    call p01_stop ( neqn, t_stop, y_stop )
  else if ( test == 2 ) then
    call p02_stop ( neqn, t_stop, y_stop )
  else if ( test == 3 ) then
    call p03_stop ( neqn, t_stop, y_stop )
  else if ( test == 4 ) then
    call p04_stop ( neqn, t_stop, y_stop )
  else if ( test == 5 ) then
    call p05_stop ( neqn, t_stop, y_stop )
  else if ( test == 6 ) then
    call p06_stop ( neqn, t_stop, y_stop )
  else if ( test == 7 ) then
    call p07_stop ( neqn, t_stop, y_stop )
  else if ( test == 8 ) then
    call p08_stop ( neqn, t_stop, y_stop )
  else if ( test == 9 ) then
    call p09_stop ( neqn, t_stop, y_stop )
  else if ( test == 10 ) then
    call p10_stop ( neqn, t_stop, y_stop )
  else if ( test == 11 ) then
    call p11_stop ( neqn, t_stop, y_stop )
  else if ( test == 12 ) then
    call p12_stop ( neqn, t_stop, y_stop )
  else if ( test == 13 ) then
    call p13_stop ( neqn, t_stop, y_stop )
  else if ( test == 14 ) then
    call p14_stop ( neqn, t_stop, y_stop )
  else if ( test == 15 ) then
    call p15_stop ( neqn, t_stop, y_stop )
  else if ( test == 16 ) then
    call p16_stop ( neqn, t_stop, y_stop )
  else if ( test == 17 ) then
    call p17_stop ( neqn, t_stop, y_stop )
  else if ( test == 18 ) then
    call p18_stop ( neqn, t_stop, y_stop )
  else if ( test == 19 ) then
    call p19_stop ( neqn, t_stop, y_stop )
  else if ( test == 20 ) then
    call p20_stop ( neqn, t_stop, y_stop )
  else if ( test == 21 ) then
    call p21_stop ( neqn, t_stop, y_stop )
  else if ( test == 22 ) then
    call p22_stop ( neqn, t_stop, y_stop )
  else if ( test == 23 ) then
    call p23_stop ( neqn, t_stop, y_stop )
  else if ( test == 24 ) then
    call p24_stop ( neqn, t_stop, y_stop )
  else if ( test == 25 ) then
    call p25_stop ( neqn, t_stop, y_stop )
  else if ( test == 26 ) then
    call p26_stop ( neqn, t_stop, y_stop )
  else if ( test == 27 ) then
    call p27_stop ( neqn, t_stop, y_stop )
  else if ( test == 28 ) then
    call p28_stop ( neqn, t_stop, y_stop )
  else if ( test == 29 ) then
    call p29_stop ( neqn, t_stop, y_stop )
  else if ( test == 30 ) then
    call p30_stop ( neqn, t_stop, y_stop )
  else if ( test == 31 ) then
    call p31_stop ( neqn, t_stop, y_stop )
  else if ( test == 32 ) then
    call p32_stop ( neqn, t_stop, y_stop )
  else if ( test == 33 ) then
    call p33_stop ( neqn, t_stop, y_stop )
  else if ( test == 34 ) then
    call p34_stop ( neqn, t_stop, y_stop )
  else if ( test == 35 ) then
    call p35_stop ( neqn, t_stop, y_stop )
  else if ( test == 36 ) then
    call p36_stop ( neqn, t_stop, y_stop )
  else if ( test == 37 ) then
    call p37_stop ( neqn, t_stop, y_stop )
  else if ( test == 38 ) then
    call p38_stop ( neqn, t_stop, y_stop )
  else if ( test == 39 ) then
    call p39_stop ( neqn, t_stop, y_stop )
  else if ( test == 40 ) then
    call p40_stop ( neqn, t_stop, y_stop )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_STOP - Fatal error!'
    write ( *, '(a,i6)' ) '  Unrecognized problem = ', test
  end if

  return
end
subroutine p00_test_num ( test_num )

!*****************************************************************************80
!
!! P00_TEST_NUM returns the number of problems available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) TEST_NUM, the number of test problems.
!
  implicit none

  integer ( kind = 4 ) test_num

  test_num = 40

  return
end
subroutine p00_title ( test, title )

!*****************************************************************************80
!
!! P00_TITLE returns the title for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TEST, the problem number.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer ( kind = 4 ) test
  character ( len = * ) title

  if ( test == 1 ) then
    call p01_title ( title )
  else if ( test == 2 ) then
    call p02_title ( title )
  else if ( test == 3 ) then
    call p03_title ( title )
  else if ( test == 4 ) then
    call p04_title ( title )
  else if ( test == 5 ) then
    call p05_title ( title )
  else if ( test == 6 ) then
    call p06_title ( title )
  else if ( test == 7 ) then
    call p07_title ( title )
  else if ( test == 8 ) then
    call p08_title ( title )
  else if ( test == 9 ) then
    call p09_title ( title )
  else if ( test == 10 ) then
    call p10_title ( title )
  else if ( test == 11 ) then
    call p11_title ( title )
  else if ( test == 12 ) then
    call p12_title ( title )
  else if ( test == 13 ) then
    call p13_title ( title )
  else if ( test == 14 ) then
    call p14_title ( title )
  else if ( test == 15 ) then
    call p15_title ( title )
  else if ( test == 16 ) then
    call p16_title ( title )
  else if ( test == 17 ) then
    call p17_title ( title )
  else if ( test == 18 ) then
    call p18_title ( title )
  else if ( test == 19 ) then
    call p19_title ( title )
  else if ( test == 20 ) then
    call p20_title ( title )
  else if ( test == 21 ) then
    call p21_title ( title )
  else if ( test == 22 ) then
    call p22_title ( title )
  else if ( test == 23 ) then
    call p23_title ( title )
  else if ( test == 24 ) then
    call p24_title ( title )
  else if ( test == 25 ) then
    call p25_title ( title )
  else if ( test == 26 ) then
    call p26_title ( title )
  else if ( test == 27 ) then
    call p27_title ( title )
  else if ( test == 28 ) then
    call p28_title ( title )
  else if ( test == 29 ) then
    call p29_title ( title )
  else if ( test == 30 ) then
    call p30_title ( title )
  else if ( test == 31 ) then
    call p31_title ( title )
  else if ( test == 32 ) then
    call p32_title ( title )
  else if ( test == 33 ) then
    call p33_title ( title )
  else if ( test == 34 ) then
    call p34_title ( title )
  else if ( test == 35 ) then
    call p35_title ( title )
  else if ( test == 36 ) then
    call p36_title ( title )
  else if ( test == 37 ) then
    call p37_title ( title )
  else if ( test == 38 ) then
    call p38_title ( title )
  else if ( test == 39 ) then
    call p39_title ( title )
  else if ( test == 40 ) then
    call p40_title ( title )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i6)' ) '  Unrecognized problem number = ', test
    stop
  end if

  return
end
function p01_autonomous ( )

!*****************************************************************************80
!
!! P01_AUTONOMOUS reports whether problem 1 is autonomous.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P01_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p01_autonomous

  p01_autonomous = .true.

  return
end
subroutine p01_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P01_EQUIL returns equilibrium solutions of problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p01_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P01_FUN evaluates the function for problem 1.
!
!  Discussion:
!
!    y' = -y
!    y(0) = 1
!
!    1 equation.
!    Enright and Pryce nonstiff problem #A1.
!    Autonomous.
!
!    Exact solution:
!
!      y(t) = exp(-t)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = - y(1)

  return
end
subroutine p01_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P01_JAC evaluates the jacobian for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = - 1.0D+00

  return
end
subroutine p01_neqn ( neqn )

!*****************************************************************************80
!
!! P01_NEQN returns the number of equations for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 1

  return
end
subroutine p01_scale ( neqn, scale )

!*****************************************************************************80
!
!! P01_SCALE returns scale factors for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1) = 1.0D+00

  return
end
subroutine p01_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P01_START returns the starting point for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = 1.0D+00

  return
end
subroutine p01_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P01_STOP returns the stopping point for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00
  y_stop(1) = 2.061153353012535D-09

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns the title of problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 1, Enright and Pryce #A1'

  return
end
function p02_autonomous ( )

!*****************************************************************************80
!
!! P02_AUTONOMOUS reports whether problem 2 is autonomous.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P02_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p02_autonomous

  p02_autonomous = .true.

  return
end
subroutine p02_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P02_EQUIL returns equilibrium solutions of problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p02_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P02_FUN evaluates the function for problem 2.
!
!  Discussion:
!
!    y' = -(y^3)/2
!    y(0) = 1
!
!    1 equation.
!    Enright and Pryce nonstiff problem #A2.
!    Autonomous.
!
!    Exact solution:
!
!      y(t) = 1 / sqrt ( t + 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = - 0.5D+00 * y(1)**3

  return
end
subroutine p02_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P02_JAC evaluates the jacobian for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = - 1.5D+00 * y(1)**2

  return
end
subroutine p02_neqn ( neqn )

!*****************************************************************************80
!
!! P02_NEQN returns the number of equations for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 1

  return
end
subroutine p02_scale ( neqn, scale )

!*****************************************************************************80
!
!! P02_SCALE returns scale factors for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1) = 1.0D+00

  return
end
subroutine p02_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P02_START returns the starting point for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = 1.0D+00

  return
end
subroutine p02_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P02_STOP returns the stopping point for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00
  y_stop(1) = 0.2182178902359887D+00

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns the title of problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 2, Enright and Pryce #A2'

  return
end
function p03_autonomous ( )

!*****************************************************************************80
!
!! P03_AUTONOMOUS reports whether problem 3 is autonomous.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P03_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p03_autonomous

  p03_autonomous = .false.

  return
end
subroutine p03_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P03_EQUIL returns equilibrium solutions of problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p03_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P03_FUN evaluates the function for problem 3.
!
!  Discussion:
!
!    y' = cos(t) * y
!    y(0) = 1
!
!    1 equation.
!    Enright and Pryce nonstiff problem #A3.
!    Not autonomous.
!
!    Exact solution:
!
!      y(t) = exp ( sin ( t ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(1) * cos ( t )

  return
end
subroutine p03_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P03_JAC evaluates the jacobian for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = cos ( t )

  return
end
subroutine p03_neqn ( neqn )

!*****************************************************************************80
!
!! P03_NEQN returns the number of equations for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 1

  return
end
subroutine p03_scale ( neqn, scale )

!*****************************************************************************80
!
!! P03_SCALE returns scale factors for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1) = 2.71D+00

  return
end
subroutine p03_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P03_START returns the starting point for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = 1.0D+00

  return
end
subroutine p03_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P03_STOP returns the stopping point for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00
  y_stop(1) = 2.491650271850414D+00

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns the title of problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 3, Enright and Pryce #A3'

  return
end
function p04_autonomous ( )

!*****************************************************************************80
!
!! P04_AUTONOMOUS reports whether problem 4 is autonomous.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P04_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p04_autonomous

  p04_autonomous = .true.

  return
end
subroutine p04_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P04_EQUIL returns equilibrium solutions of problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1) = 0.0D+00
  else if ( next == 1 ) then
    next = 2
    y(1) = 20.0D+00
  else
    next = 0
  end if

  return
end
subroutine p04_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P04_FUN evaluates the function for problem 4.
!
!  Discussion:
!
!    y' = y*(20-y)/80
!    y(0) = 1
!
!    1 equation.
!    Enright and Pryce nonstiff problem #A4.
!    Autonomous.
!
!    Exact solution:
!
!      y(t) = 20 / ( 1 + 19 * exp ( -t / 4 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(1) * ( 20.0D+00 - y(1) ) / 80.0D+00

  return
end
subroutine p04_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P04_JAC evaluates the jacobian for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = ( 10.0D+00 - y(1) ) / 40.0D+00

  return
end
subroutine p04_neqn ( neqn )

!*****************************************************************************80
!
!! P04_NEQN returns the number of equations for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 1

  return
end
subroutine p04_scale ( neqn, scale )

!*****************************************************************************80
!
!! P04_SCALE returns scale factors for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1) = 17.7D+00

  return
end
subroutine p04_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P04_START returns the starting point for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = 1.0D+00

  return
end
subroutine p04_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P04_STOP returns the stopping point for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00
  y_stop(1) = 17.73016648131483D+00

  return
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns the title of problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 4, Enright and Pryce #A4'

  return
end
function p05_autonomous ( )

!*****************************************************************************80
!
!! P05_AUTONOMOUS reports whether problem 5 is autonomous.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P05_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p05_autonomous

  p05_autonomous = .false.

  return
end
subroutine p05_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P05_EQUIL returns equilibrium solutions of problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  next = 0

  return
end
subroutine p05_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P05_FUN evaluates the function for problem 5.
!
!  Discussion:
!
!    y' = (y-t)/(y+t)
!    y(0) = 4
!
!    1 equation.
!    Enright and Pryce nonstiff problem #A5.
!    Not autonomous.
!
!    Exact solution:
!
!      r = sqrt ( t + y(t)**2 )
!      theta = atan ( y(t) / t )
!
!      r = 4 * exp ( pi/2 - theta )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = ( y(1) - t ) / ( y(1) + t )

  return
end
subroutine p05_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P05_JAC evaluates the jacobian for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = 2.0D+00 * t / ( y(1) + t )**2

  return
end
subroutine p05_neqn ( neqn )

!*****************************************************************************80
!
!! P05_NEQN returns the number of equations for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 1

  return
end
subroutine p05_scale ( neqn, scale )

!*****************************************************************************80
!
!! P05_SCALE returns scale factors for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1) = 6.2D+00

  return
end
subroutine p05_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P05_START returns the starting point for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = 4.0D+00

  return
end
subroutine p05_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P05_STOP returns the stopping point for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00
  y_stop(1) = - 0.788782668896419D+00

  return
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns the title of problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 5, Enright and Pryce #A5'

  return
end
function p06_autonomous ( )

!*****************************************************************************80
!
!! P06_AUTONOMOUS reports whether problem 6 is autonomous.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P06_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p06_autonomous

  p06_autonomous = .true.

  return
end
subroutine p06_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P06_EQUIL returns equilibrium solutions of problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:2) = (/ 0.0D+00, 0.0D+00 /)
  else if ( next == 1 ) then
    next = 2
    y(1:2) = (/ 1.0D+00, 1.0D+00 /)
  else
    next = 0
  end if

  return
end
subroutine p06_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P06_FUN evaluates the function for problem 6.
!
!  Discussion:
!
!    y1' = 2 y1 * ( 1 - y2 )
!    y2' = - y2 * ( 1 - y1 )
!    y1(0) = 1
!    y2(0) = 3
!
!    2 equations.
!    Enright and Pryce nonstiff problem #B1.
!    Autonomous.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = 2.0D+00 * y(1) * ( 1.0D+00 - y(2) )
  yp(2) =         - y(2) * ( 1.0D+00 - y(1) )

  return
end
subroutine p06_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P06_JAC evaluates the jacobian for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = 2.0D+00 * ( 1.0D+00 - y(2) )
  jac(1,2) = 2.0D+00 * y(1)

  jac(2,1) = y(2)
  jac(2,2) = - ( 1.0D+00 - y(1) )

  return
end
subroutine p06_neqn ( neqn )

!*****************************************************************************80
!
!! P06_NEQN returns the number of equations for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p06_scale ( neqn, scale )

!*****************************************************************************80
!
!! P06_SCALE returns scale factors for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:2) = (/ 4.25D+00, 3.00D+00 /)

  return
end
subroutine p06_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P06_START returns the starting point for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = (/ 1.0D+00, 3.0D+00 /)

  return
end
subroutine p06_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P06_STOP returns the stopping point for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00
  y_stop(1:2) = (/ 0.6761876008576667D+00, 0.1860816099640036D+00 /)

  return
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns the title of problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 6, Enright and Pryce #B1'

  return
end
function p07_autonomous ( )

!*****************************************************************************80
!
!! P07_AUTONOMOUS reports whether problem 7 is autonomous.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P07_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p07_autonomous

  p07_autonomous = .true.

  return
end
subroutine p07_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P07_EQUIL returns equilibrium solutions of problem 7.
!
!  Discussion:
!
!    Any solution with Y(1) = Y(2) = Y(3) is an equilibrium
!    solution.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:3) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  else
    next = 0
  end if

  return
end
subroutine p07_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P07_FUN evaluates the function for problem 7.
!
!  Discussion:
!
!    y1' = -y1 +   y2
!    y2' =  y1 - 2 y2 + y3
!    y3' =         y2 - y3
!    y1(0) = 2
!    y2(0) = 0
!    y3(0) = 1
!
!    3 equations.
!    Enright and Pryce nonstiff problem #B2.
!    Autonomous.
!
!    Note that the quantity (y1+y2+y3) is conserved by the exact solution.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = - y(1) +           y(2)
  yp(2) =   y(1) - 2.0D+00 * y(2) + y(3)
  yp(3) =                    y(2) - y(3)

  return
end
subroutine p07_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P07_JAC evaluates the jacobian for problem 7.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = - 1.0D+00
  jac(1,2) = + 1.0D+00
  jac(1,3) =   0.0D+00

  jac(2,1) = + 1.0D+00
  jac(2,2) = - 2.0D+00
  jac(2,3) = + 1.0D+00

  jac(3,1) =   0.0D+00
  jac(3,2) = + 1.0D+00
  jac(3,3) = - 1.0D+00

  return
end
subroutine p07_neqn ( neqn )

!*****************************************************************************80
!
!! P07_NEQN returns the number of equations for problem 7.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 3

  return
end
subroutine p07_scale ( neqn, scale )

!*****************************************************************************80
!
!! P07_SCALE returns scale factors for problem 7.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:3) = (/ 2.0D+00, 1.0D+00, 1.0D+00 /)

  return
end
subroutine p07_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P07_START returns the starting point for problem 7.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn) = (/ 2.0D+00, 0.0D+00, 1.0D+00 /)

  return
end
subroutine p07_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P07_STOP returns the stopping point for problem 7.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:3) = (/ 1.000000001030576D+00, &
                   1.000000000000000D+00, &
                   0.9999999989694235D+00 /)

  return
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns the title of problem 7.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 7, Enright and Pryce #B2'

  return
end
function p08_autonomous ( )

!*****************************************************************************80
!
!! P08_AUTONOMOUS reports whether problem 8 is autonomous.
!
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P08_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p08_autonomous

  p08_autonomous = .true.

  return
end
subroutine p08_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P08_EQUIL returns equilibrium solutions of problem 8.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:3) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  else
    next = 0
  end if

  return
end
subroutine p08_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P08_FUN evaluates the function for problem 8.
!
!  Discussion:
!
!    y1' = - y1
!    y2' =   y1 - y2^2
!    y3' =        y2^2
!    y1(0) = 1
!    y2(0) = 0
!    y3(0) = 0
!
!    3 equations.
!    Enright and Pryce nonstiff problem #B3.
!    Autonomous.
!
!    Notice that the quantity (y1+y2+y3) is conserved by the exact solution.
!
!  Modified:
!
!    07 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1:3) = (/ -y(1), y(1) - y(2)**2, y(2)**2 /)

  return
end
subroutine p08_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P08_JAC evaluates the jacobian for problem 8.
!
!  Modified:
!
!    07 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = - 1.0D+00
  jac(1,2) =   0.0D+00
  jac(1,3) =   0.0D+00

  jac(2,1) =   1.0D+00
  jac(2,2) = - 2.0D+00 * y(2)
  jac(2,3) =   0.0D+00

  jac(3,1) =   0.0D+00
  jac(3,2) = - 2.0D+00 * y(2)
  jac(3,3) =   0.0D+00

  return
end
subroutine p08_neqn ( neqn )

!*****************************************************************************80
!
!! P08_NEQN returns the number of equations for problem 8.
!
!  Modified:
!
!    07 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 3

  return
end
subroutine p08_scale ( neqn, scale )

!*****************************************************************************80
!
!! P08_SCALE returns scale factors for problem 8.
!
!  Modified:
!
!    07 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:3) = (/ 1.0D+00, 0.519D+00, 0.947D+00 /)

  return
end
subroutine p08_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P08_START returns the starting point for problem 8.
!
!  Modified:
!
!    07 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)

  return
end
subroutine p08_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P08_STOP returns the stopping point for problem 8.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:3) = (/ 2.061153488557776D-09, &
                   5.257228022048349D-02, &
                   9.474277177183630D-01 /)

  return
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns the title of problem 8.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 8, Enright and Pryce #B3'

  return
end
function p09_autonomous ( )

!*****************************************************************************80
!
!! P09_AUTONOMOUS reports whether problem 9 is autonomous.
!
!  Modified:
!
!    02 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P09_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p09_autonomous

  p09_autonomous = .true.

  return
end
subroutine p09_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P09_EQUIL returns equilibrium solutions of problem 9.
!
!  Modified:
!
!    02 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:3) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  else
    next = 0
  end if

  return
end
subroutine p09_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P09_FUN evaluates the function for problem 9.
!
!  Discussion:
!
!    y1' =  ( - y2 - y1 * y3 ) / sqrt ( y1^2 + y2^2 )
!    y2' =  (   y1 - y2 * y3 ) / sqrt ( y1^2 + y2^2 )
!    y3' =      y1             / sqrt ( y1^2 + y2^2 )
!    y1(0) = 3
!    y2(0) = 0
!    y3(0) = 0
!
!    3 equations.
!    Enright and Pryce nonstiff problem #B4.
!    Autonomous.
!
!  Discussion:
!
!    Enright and Pryce nonstiff problem #B4.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) norm
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  norm = sqrt ( y(1)**2 + y(2)**2 )

  if ( 0.0D+00 < norm ) then
    yp(1) = - y(2) - y(1) * y(3) / norm
    yp(2) =   y(1) - y(2) * y(3) / norm
    yp(3) =                 y(1) / norm
  else
    yp(1:3) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  end if

  return
end
subroutine p09_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P09_JAC evaluates the jacobian for problem 9.
!
!  Modified:
!
!    02 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) norm
  real ( kind = 8 ) norm3
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  norm = sqrt ( y(1)**2 + y(2)**2 )
  norm3 = norm**3

  jac(1,1) = - y(3) * y(2)**2 / norm3
  jac(1,2) = - 1.0D+00 + y(1) * y(2) * y(3) / norm3
  jac(1,3) = - y(1) / norm

  jac(2,1) = 1.0D+00 + y(1) * y(2) * y(3) / norm3
  jac(2,2) = - y(1)**2 * y(3) / norm3
  jac(2,3) = - y(2) / norm

  jac(3,1) =   y(2)**2 / norm3
  jac(3,2) = - y(1) * y(2) / norm3
  jac(3,3) = 0.0D+00

  return
end
subroutine p09_neqn ( neqn )

!*****************************************************************************80
!
!! P09_NEQN returns the number of equations for problem 9.
!
!  Modified:
!
!    02 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 3

  return
end
subroutine p09_scale ( neqn, scale )

!*****************************************************************************80
!
!! P09_SCALE returns scale factors for problem 9.
!
!  Modified:
!
!    02 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:3) = (/ 1.0D+00, 1.0D+00, 1.0D+00 /)

  return
end
subroutine p09_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P09_START returns the starting point for problem 9.
!
!  Modified:
!
!    02 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = (/ 3.0D+00, 0.0D+00, 0.0D+00 /)

  return
end
subroutine p09_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P09_STOP returns the stopping point for problem 9.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:3) = (/ &
    9.826950928005993D-01, &
    2.198447081694832D+00, &
    9.129452507276399D-01 /)

  return
end
subroutine p09_title ( title )

!*****************************************************************************80
!
!! P09_TITLE returns the title of problem 9.
!
!  Modified:
!
!    02 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 9, Enright and Pryce #B4'

  return
end
function p10_autonomous ( )

!*****************************************************************************80
!
!! P10_AUTONOMOUS reports whether problem 10 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P10_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p10_autonomous

  p10_autonomous = .true.

  return
end
subroutine p10_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P10_EQUIL returns equilibrium solutions of problem 10.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:3) = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
  else if ( next == 1 ) then
    next = 2
    y(1:3) = (/ 0.0D+00, 1.0D+00, 0.0D+00 /)
  else if ( next == 2 ) then
    next = 3
    y(1:3) = (/ 0.0D+00, 0.0D+00, 1.0D+00 /)
  else
    next = 0
  end if

  return
end
subroutine p10_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P10_FUN evaluates the function for problem 10.
!
!  Discussion:
!
!    y1' =   y2 * y3
!    y2' = - y1 * y3
!    y3' = - 0.51 * y1 * y2
!    y1(0) = 0
!    y2(0) = 1
!    y3(0) = 1
!
!    3 equations.
!    Enright and Pryce nonstiff problem #B5.
!    Autonomous.
!
!  Discussion:
!
!    Enright and Pryce nonstiff problem #B5.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) =   y(2) * y(3)
  yp(2) = - y(1) * y(3)
  yp(3) = - 0.51D+00 * y(1) * y(2)

  return
end
subroutine p10_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P10_JAC evaluates the jacobian for problem 10.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1:3) = (/  0.00D+00,  1.00D+00,  1.00D+00 /)
  jac(2,1:3) = (/ -1.00D+00,  0.00D+00, -1.00D+00 /)
  jac(3,1:3) = (/ -0.51D+00, -0.51D+00,  0.00D+00 /)

  return
end
subroutine p10_neqn ( neqn )

!*****************************************************************************80
!
!! P10_NEQN returns the number of equations for problem 10.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 3

  return
end
subroutine p10_scale ( neqn, scale )

!*****************************************************************************80
!
!! P10_SCALE returns scale factors for problem 10.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:3) = (/ 1.0D+00, 1.0D+00, 1.0D+00 /)

  return
end
subroutine p10_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P10_START returns the starting point for problem 10.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = (/ 0.0D+00, 1.0D+00, 1.0D+00 /)

  return
end
subroutine p10_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P10_STOP returns the stopping point for problem 10.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:3) = (/ &
    -9.396570798729192D-01, &
    -3.421177754000779D-01, &
     7.414126596199957D-01 /)

  return
end
subroutine p10_title ( title )

!*****************************************************************************80
!
!! P10_TITLE returns the title of problem 10.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 10, Enright and Pryce #B5'

  return
end
function p11_autonomous ( )

!*****************************************************************************80
!
!! P11_AUTONOMOUS reports whether problem 11 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P11_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p11_autonomous

  p11_autonomous = .true.

  return
end
subroutine p11_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P11_EQUIL returns equilibrium solutions of problem 11.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p11_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P11_FUN evaluates the function for problem 11.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1)        =             - y(1)
  yp(2:neqn-1) = y(1:neqn-2) - y(2:neqn-1)
  yp(neqn)     = y(neqn-1)

  return
end
subroutine p11_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P11_JAC evaluates the jacobian for problem 11.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) i
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1:neqn,1:neqn) = 0.0D+00

  do i = 1, neqn-1
    jac(i,i) = -1.0D+00
  end do

  do i = 2, neqn
    jac(i,i-1) = 1.0D+00
  end do

  return
end
subroutine p11_neqn ( neqn )

!*****************************************************************************80
!
!! P11_NEQN returns the number of equations for problem 11.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 10

  return
end
subroutine p11_scale ( neqn, scale )

!*****************************************************************************80
!
!! P11_SCALE returns scale factors for problem 11.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p11_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P11_START returns the starting point for problem 11.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1) = 1.0D+00
  y_start(2:neqn) = 0.0D+00

  return
end
subroutine p11_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P11_STOP returns the stopping point for problem 11.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ &
    2.061153622240064D-09, &
    4.122307244619555D-08, &
    4.122307244716968D-07, &
    2.748204829855288D-06, &
    1.374102414941961D-05, &
    5.496409659803266D-05, &
    1.832136553274552D-04, &
    5.234675866508716D-04, &
    1.308668966628220D-03, &
    9.979127409508656D-01 /)

  return
end
subroutine p11_title ( title )

!*****************************************************************************80
!
!! P11_TITLE returns the title of problem 11.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 11, Enright and Pryce #C1'

  return
end
function p12_autonomous ( )

!*****************************************************************************80
!
!! P12_AUTONOMOUS reports whether problem 12 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P12_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p12_autonomous

  p12_autonomous = .true.

  return
end
subroutine p12_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P12_EQUIL returns equilibrium solutions of problem 12.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p12_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P12_FUN evaluates the function for problem 12.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) i
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1)        =             - y(1)
  do i = 2, neqn-1
    yp(i) = real ( i - 1, kind = 8 ) * y(i-1) - real ( i, kind = 8 ) * y(i)
  end do
  yp(neqn) = real ( neqn - 1, kind = 8 ) * y(neqn-1)

  return
end
subroutine p12_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P12_JAC evaluates the jacobian for problem 12.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) i
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1:neqn,1:neqn) = 0.0D+00

  do i = 1, neqn-1
    jac(i,i) = - real ( i, kind = 8 )
  end do

  do i = 2, neqn
    jac(i,i-1) = real ( i - 1, kind = 8 )
  end do

  return
end
subroutine p12_neqn ( neqn )

!*****************************************************************************80
!
!! P12_NEQN returns the number of equations for problem 12.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 10

  return
end
subroutine p12_scale ( neqn, scale )

!*****************************************************************************80
!
!! P12_SCALE returns scale factors for problem 12.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p12_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P12_START returns the starting point for problem 12.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1) = 1.0D+00
  y_start(2:neqn) = 0.0D+00

  return
end
subroutine p12_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P12_STOP returns the stopping point for problem 12.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ &
    2.061153577984930D-09, &
    2.061153573736588D-09, &
    2.061153569488245D-09, &
    2.061153565239902D-09, &
    2.061153560991560D-09, &
    2.061153556743217D-09, &
    2.061153552494874D-09, &
    2.061153548246532D-09, &
    2.061153543998189D-09, &
    9.999999814496180D-01 /)

  return
end
subroutine p12_title ( title )

!*****************************************************************************80
!
!! P12_TITLE returns the title of problem 12.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 12, Enright and Pryce #C2'

  return
end
function p13_autonomous ( )

!*****************************************************************************80
!
!! P13_AUTONOMOUS reports whether problem 13 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P13_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p13_autonomous

  p13_autonomous = .true.

  return
end
subroutine p13_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P13_EQUIL returns equilibrium solutions of problem 13.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p13_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P13_FUN evaluates the function for problem 13.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1)        =             - 2.0D+00 * y(1) + y(2)
  yp(2:neqn-1) = y(1:neqn-2) - 2.0D+00 * y(2:neqn-1) + y(3:neqn)
  yp(neqn)     = y(neqn-1)   - 2.0D+00 * y(neqn)

  return
end
subroutine p13_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P13_JAC evaluates the jacobian for problem 13.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) i
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1:neqn,1:neqn) = 0.0D+00

  do i = 2, neqn
    jac(i,i-1) = 1.0D+00
  end do

  do i = 1, neqn
    jac(i,i) = -2.0D+00
  end do

  do i = 1, neqn-1
    jac(i,i+1) = 1.0D+00
  end do

  return
end
subroutine p13_neqn ( neqn )

!*****************************************************************************80
!
!! P13_NEQN returns the number of equations for problem 13.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 10

  return
end
subroutine p13_scale ( neqn, scale )

!*****************************************************************************80
!
!! P13_SCALE returns scale factors for problem 13.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p13_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P13_START returns the starting point for problem 13.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1) = 1.0D+00
  y_start(2:neqn) = 0.0D+00

  return
end
subroutine p13_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P13_STOP returns the stopping point for problem 13.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ &
    2.948119211022058D-03, &
    5.635380154844266D-03, &
    7.829072515926013D-03, &
    9.348257908594937D-03, &
    1.007943610301970D-02, &
    9.982674171429909D-03, &
    9.088693332766085D-03, &
    7.489115195185912D-03, &
    5.322964130953349D-03, &
    2.762434379029886D-03 /)

  return
end
subroutine p13_title ( title )

!*****************************************************************************80
!
!! P13_TITLE returns the title of problem 13.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 13, Enright and Pryce #C3'

  return
end
function p14_autonomous ( )

!*****************************************************************************80
!
!! P14_AUTONOMOUS reports whether problem 14 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P14_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p14_autonomous

  p14_autonomous = .true.

  return
end
subroutine p14_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P14_EQUIL returns equilibrium solutions of problem 14.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p14_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P14_FUN evaluates the function for problem 14.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1)        =             - 2.0D+00 * y(1) + y(2)
  yp(2:neqn-1) = y(1:neqn-2) - 2.0D+00 * y(2:neqn-1) + y(3:neqn)
  yp(neqn)     = y(neqn-1)   - 2.0D+00 * y(neqn)

  return
end
subroutine p14_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P14_JAC evaluates the jacobian for problem 14.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) i
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1:neqn,1:neqn) = 0.0D+00

  do i = 2, neqn
    jac(i,i-1) = 1.0D+00
  end do

  do i = 1, neqn
    jac(i,i) = -2.0D+00
  end do

  do i = 1, neqn-1
    jac(i,i+1) = 1.0D+00
  end do

  return
end
subroutine p14_neqn ( neqn )

!*****************************************************************************80
!
!! P14_NEQN returns the number of equations for problem 14.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 51

  return
end
subroutine p14_scale ( neqn, scale )

!*****************************************************************************80
!
!! P14_SCALE returns scale factors for problem 14.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p14_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P14_START returns the starting point for problem 14.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1) = 1.0D+00
  y_start(2:neqn) = 0.0D+00

  return
end
subroutine p14_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P14_STOP returns the stopping point for problem 14.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ &
    3.124111453721466D-03, &
    6.015416842150318D-03, &
    8.470021834842650D-03, &
    1.033682931733337D-02, &
    1.153249572873923D-02, &
    1.204549525737964D-02, &
    1.192957068015293D-02, &
    1.128883207111195D-02, &
    1.025804501391024D-02, &
    8.982017581934167D-03, &
    7.597500902492453D-03, &
    6.219920556824985D-03, &
    4.935916341009131D-03, &
    3.801432544256119D-03, &
    2.844213677587894D-03, &
    2.069123394222672D-03, &
    1.464687282843915D-03, &
    1.009545263941126D-03, &
    6.779354330227017D-04, &
    4.437815269118510D-04, &
    2.833264542938954D-04, &
    1.765005798796805D-04, &
    1.073342592697238D-04, &
    6.374497601777217D-05, &
    3.698645309704183D-05, &
    2.097466832643746D-05, &
    1.162956710412555D-05, &
    6.306710405783322D-06, &
    3.346286430868515D-06, &
    1.737760074184334D-06, &
    8.835366904275847D-07, &
    4.399520411127637D-07, &
    2.146181897152360D-07, &
    1.025981211654928D-07, &
    4.807864068784215D-08, &
    2.209175152474847D-08, &
    9.956251263138180D-09, &
    4.402193653748924D-09, &
    1.910149382204028D-09, &
    8.135892921473050D-10, &
    3.402477118549235D-10, &
    1.397485617545782D-10, &
    5.638575303049199D-11, &
    2.235459707956947D-11, &
    8.710498036398032D-12, &
    3.336554275346643D-12, &
    1.256679567784939D-12, &
    4.654359053128788D-13, &
    1.693559145599857D-13, &
    5.996593816663054D-14, &
    1.891330702629865D-14 /)

  return
end
subroutine p14_title ( title )

!*****************************************************************************80
!
!! P14_TITLE returns the title of problem 14.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 14, Enright and Pryce #C4'

  return
end
function p15_autonomous ( )

!*****************************************************************************80
!
!! P15_AUTONOMOUS reports whether problem 15 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P15_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p15_autonomous

  p15_autonomous = .true.

  return
end
subroutine p15_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P15_EQUIL returns equilibrium solutions of problem 15.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p15_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P15_FUN evaluates the function for problem 15.
!
!  Discussion:
!
!    This system models the motion of the five outer planets of the
!    solar system.
!
!  Modified:
!
!    19 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j
  real ( kind = 8 ) k2
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  real ( kind = 8 ) m0
  real ( kind = 8 ) m(5)
  integer ( kind = 4 ) mm
  real ( kind = 8 ) p
  real ( kind = 8 ) q(5,5)
  real ( kind = 8 ) r(5)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  call p15_param ( 'GET', 'K2', k2 )
  call p15_param ( 'GET', 'M0', m0 )
  call p15_param ( 'GET', 'M1', m(1) )
  call p15_param ( 'GET', 'M2', m(2) )
  call p15_param ( 'GET', 'M3', m(3) )
  call p15_param ( 'GET', 'M4', m(4) )
  call p15_param ( 'GET', 'M5', m(5) )

  i = 0
  do l = 3, 15, 3
    i = i + 1
    p = y(l-2)**2 + y(l-1)**2 + y(l)**2
    r(i) = 1.0D+00 / ( p * sqrt ( p ) )
    j = 0
    do ll = 3, 15, 3
      j = j + 1
      if ( ll /= l ) then
        p = ( y(l-2) - y(ll-2) )**2 + ( y(l-1) - y(ll-1) )**2 &
          + ( y(l) - y(ll) )**2
        q(i,j) = 1.0D+00 / ( p * sqrt ( p ) )
        q(j,i) = q(i,j)
      end if
    end do
  end do

  i3 = 0
  do i = 1, 5
    i3 = i3 + 3
    do ll = i3-2, i3
      mm = ll - i3
      yp(ll) = y(ll+15)
      p = 0.0D+00
      do j = 1, 5
        mm = mm + 3
        if ( j /= i ) then
          p = p + m(j) &
             * ( y(mm) * ( q(i,j) - r(j) ) - y(ll) * q(i,j) )
        end if
      end do
      yp(ll+15) = k2 * ( - ( m0 + m(i) ) * y(ll) * r(i) + p )
    end do
  end do

  return
end
subroutine p15_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P15_JAC evaluates the jacobian for problem 15.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) k2
  real ( kind = 8 ) m0
  real ( kind = 8 ) m1
  real ( kind = 8 ) m2
  real ( kind = 8 ) m3
  real ( kind = 8 ) m4
  real ( kind = 8 ) m5
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  call p15_param ( 'GET', 'K2', k2 )
  call p15_param ( 'GET', 'M0', m0 )
  call p15_param ( 'GET', 'M1', m1 )
  call p15_param ( 'GET', 'M2', m2 )
  call p15_param ( 'GET', 'M3', m3 )
  call p15_param ( 'GET', 'M4', m4 )
  call p15_param ( 'GET', 'M5', m5 )

  jac(1:neqn,1:neqn) = 0.0D+00

  return
end
subroutine p15_neqn ( neqn )

!*****************************************************************************80
!
!! P15_NEQN returns the number of equations for problem 15.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 30

  return
end
subroutine p15_param ( action, name, value )

!*****************************************************************************80
!
!! P15_PARAM handles the parameters for problem 15.
!
!  Discussion:
!
!    K2 is the gravitational constant;
!    M0 is the lumped mass of the sun and four inner planets;
!    M1 is the mass of Jupiter;
!    M2 is the mass of Saturn;
!    M3 is the mass of Uranus;
!    M4 is the mass of Neptune;
!    M5 is the mass of Pluto.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'K2', 'M0', 'M1', 'M2', 'M3', 'M4' or 'M5'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: k2 = 2.95912208286D+00
  real ( kind = 8 ), save :: m0 = 1.00000597682D+00
  real ( kind = 8 ), save :: m1 = 0.954786104043D-03
  real ( kind = 8 ), save :: m2 = 0.285583733151D-03
  real ( kind = 8 ), save :: m3 = 0.437273164546D-04
  real ( kind = 8 ), save :: m4 = 0.517759138449D-04
  real ( kind = 8 ), save :: m5 = 0.277777777778D-05
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'K2' ) ) then
      value = k2
    else if ( s_eqi ( name, 'M0' ) ) then
      value = m0
    else if ( s_eqi ( name, 'M1' ) ) then
      value = m1
    else if ( s_eqi ( name, 'M2' ) ) then
      value = m2
    else if ( s_eqi ( name, 'M3' ) ) then
      value = m3
    else if ( s_eqi ( name, 'M4' ) ) then
      value = m4
    else if ( s_eqi ( name, 'M5' ) ) then
      value = m5
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P15_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'K2' ) ) then
      k2 = value
    else if ( s_eqi ( name, 'M0' ) ) then
      m0 = value
    else if ( s_eqi ( name, 'M1' ) ) then
      m1 = value
    else if ( s_eqi ( name, 'M2' ) ) then
      m2 = value
    else if ( s_eqi ( name, 'M3' ) ) then
      m3 = value
    else if ( s_eqi ( name, 'M4' ) ) then
      m4 = value
    else if ( s_eqi ( name, 'M5' ) ) then
      m5 = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P15_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P15_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p15_scale ( neqn, scale )

!*****************************************************************************80
!
!! P15_SCALE returns scale factors for problem 15.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p15_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P15_START returns the starting point for problem 15.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn)= (/ &
     3.42947415189D+00, &
     3.35386959711D+00, &
     1.35494901715D+00, &
     6.64145542550D+00, &
     5.97156957878D+00, &
     2.18231499728D+00, &
    11.2630437207D+00, &
    14.6952576794D+00, &
     6.27960525067D+00, &
   -30.1552268759D+00, &
     1.65699966404D+00, &
     1.43785752721D+00, &
   -21.1238353380D+00, &
    28.4465098142D+00, &
    15.3882659679D+00, &
    -0.557160570446D+00, &
     0.505696783289D+00, &
     0.230578543901D+00, &
    -0.415570776342D+00, &
     0.365682722812D+00, &
     0.169143213293D+00, &
    -0.325325669158D+00, &
     0.189706021964D+00, &
     0.0877265322780D+00, &
    -0.0240476254170D+00, &
    -0.287659532608D+00, &
    -0.117219543175D+00, &
    -0.176860753121D+00, &
    -0.216393453025D+00, &
    -0.0148647893090D+00 /)

  return
end
subroutine p15_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P15_STOP returns the stopping point for problem 15.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ &
    -4.792730224323733D+00, &
    -2.420550725448973D+00, &
    -9.212509306014886D-01, &
    -4.217310404035213D+00, &
     7.356202947498970D+00, &
     3.223785985421212D+00, &
     4.035559443262270D+00, &
     1.719865528670555D+01, &
     7.478910794233703D+00, &
    -2.998759326324844D+01, &
    -4.107310937550929D+00, &
    -9.277008321754407D-01, &
    -2.442125302518482D+01, &
     2.381459045746554D+01, &
     1.492096306951359D+01, &
     3.499208963063806D-01, &
    -5.748487687912825D-01, &
    -2.551694020879149D-01, &
    -5.237040978903326D-01, &
    -2.493000463579661D-01, &
    -8.045341642044464D-02, &
    -3.875289237334110D-01, &
     5.648603288767891D-02, &
     3.023606472143342D-02, &
     4.133856546712445D-02, &
    -2.862393029841379D-01, &
    -1.183032405136207D-01, &
    -1.511986457359206D-01, &
    -2.460068894318766D-01, &
    -3.189687411323877D-02 /)

  return
end
subroutine p15_title ( title )

!*****************************************************************************80
!
!! P15_TITLE returns the title of problem 15.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 15, Enright and Pryce #C5'

  return
end
function p16_autonomous ( )

!*****************************************************************************80
!
!! P16_AUTONOMOUS reports whether problem 16 is autonomous.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P16_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p16_autonomous

  p16_autonomous = .true.

  return
end
subroutine p16_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P16_EQUIL returns equilibrium solutions of problem 16.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p16_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P16_FUN evaluates the function for problem 16.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(3)
  yp(2) = y(4)
  yp(3) = - y(1) / ( sqrt ( ( y(1)**2 + y(2)**2 ) ) )**3
  yp(4) = - y(2) / ( sqrt ( ( y(1)**2 + y(2)**2 ) ) )**3

  return
end
subroutine p16_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P16_JAC evaluates the jacobian for problem 16.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1:neqn,1:neqn) = 0.0D+00

  return
end
subroutine p16_neqn ( neqn )

!*****************************************************************************80
!
!! P16_NEQN returns the number of equations for problem 16.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 4

  return
end
subroutine p16_param ( action, name, value )

!*****************************************************************************80
!
!! P16_PARAM handles the parameters for problem 16.
!
!  Modified:
!
!    20 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'DELTA'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: delta = 0.1D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      value = delta
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P16_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      delta = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P16_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P16_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p16_scale ( neqn, scale )

!*****************************************************************************80
!
!! P16_SCALE returns scale factors for problem 16.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p16_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P16_START returns the starting point for problem 16.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) delta
  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  call p16_param ( 'GET', 'DELTA', delta )

  y_start(1) = 1.0D+00 - delta
  y_start(2) = 0.0D+00
  y_start(3) = 0.0D+00
  y_start(4) = sqrt ( ( 1.0D+00 + delta ) / ( 1.0D+00 - delta ) )

  return
end
subroutine p16_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P16_STOP returns the stopping point for problem 16.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ &
    2.198835352008397D-01, &
    9.427076846341813D-01, &
   -9.787659841058176D-01, &
    3.287977990962036D-01 /)

  return
end
subroutine p16_title ( title )

!*****************************************************************************80
!
!! P16_TITLE returns the title of problem 16.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 16, Enright and Pryce #D1'

  return
end
function p17_autonomous ( )

!*****************************************************************************80
!
!! P17_AUTONOMOUS reports whether problem 17 is autonomous.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P17_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p17_autonomous

  p17_autonomous = .true.

  return
end
subroutine p17_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P17_EQUIL returns equilibrium solutions of problem 17.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p17_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P17_FUN evaluates the function for problem 17.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(3)
  yp(2) = y(4)
  yp(3) = - y(1) / ( sqrt ( ( y(1)**2 + y(2)**2 ) ) )**3
  yp(4) = - y(2) / ( sqrt ( ( y(1)**2 + y(2)**2 ) ) )**3

  return
end
subroutine p17_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P17_JAC evaluates the jacobian for problem 17.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1:neqn,1:neqn) = 0.0D+00

  return
end
subroutine p17_neqn ( neqn )

!*****************************************************************************80
!
!! P17_NEQN returns the number of equations for problem 17.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 4

  return
end
subroutine p17_param ( action, name, value )

!*****************************************************************************80
!
!! P17_PARAM handles the parameters for problem 17.
!
!  Modified:
!
!    20 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'DELTA'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: delta = 0.3D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      value = delta
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P17_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      delta = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P17_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P17_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p17_scale ( neqn, scale )

!*****************************************************************************80
!
!! P17_SCALE returns scale factors for problem 17.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p17_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P17_START returns the starting point for problem 17.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) delta
  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  call p17_param ( 'GET', 'DELTA', delta )

  y_start(1) = 1.0D+00 - delta
  y_start(2) = 0.0D+00
  y_start(3) = 0.0D+00
  y_start(4) = sqrt ( ( 1.0D+00 + delta ) / ( 1.0D+00 - delta ) )

  return
end
subroutine p17_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P17_STOP returns the stopping point for problem 17.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ &
   -1.777027357140412D-01, &
    9.467784719905892D-01, &
   -1.030294163192969D+00, &
    1.211074890053952D-01 /)

  return
end
subroutine p17_title ( title )

!*****************************************************************************80
!
!! P17_TITLE returns the title of problem 17.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 17, Enright and Pryce #D2'

  return
end
function p18_autonomous ( )

!*****************************************************************************80
!
!! P18_AUTONOMOUS reports whether problem 18 is autonomous.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P18_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p18_autonomous

  p18_autonomous = .true.

  return
end
subroutine p18_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P18_EQUIL returns equilibrium solutions of problem 18.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p18_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P18_FUN evaluates the function for problem 18.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(3)
  yp(2) = y(4)
  yp(3) = - y(1) / ( sqrt ( ( y(1)**2 + y(2)**2 ) ) )**3
  yp(4) = - y(2) / ( sqrt ( ( y(1)**2 + y(2)**2 ) ) )**3

  return
end
subroutine p18_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P18_JAC evaluates the jacobian for problem 18.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1:neqn,1:neqn) = 0.0D+00

  return
end
subroutine p18_neqn ( neqn )

!*****************************************************************************80
!
!! P18_NEQN returns the number of equations for problem 18.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 4

  return
end
subroutine p18_param ( action, name, value )

!*****************************************************************************80
!
!! P18_PARAM handles the parameters for problem 18.
!
!  Modified:
!
!    20 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'DELTA'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: delta = 0.5D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      value = delta
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P18_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      delta = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P18_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P18_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p18_scale ( neqn, scale )

!*****************************************************************************80
!
!! P18_SCALE returns scale factors for problem 18.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p18_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P18_START returns the starting point for problem 18.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) delta
  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  delta = 0.5D+00

  y_start(1) = 1.0D+00 - delta
  y_start(2) = 0.0D+00
  y_start(3) = 0.0D+00
  y_start(4) = sqrt ( ( 1.0D+00 + delta ) / ( 1.0D+00 - delta ) )

  return
end
subroutine p18_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P18_STOP returns the stopping point for problem 18.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ &
   -5.780432953035361D-01, &
    8.633840009194193D-01, &
   -9.595083730380727D-01, &
   -6.504915126712089D-02 /)

  return
end
subroutine p18_title ( title )

!*****************************************************************************80
!
!! P18_TITLE returns the title of problem 18.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 18, Enright and Pryce #D3'

  return
end
function p19_autonomous ( )

!*****************************************************************************80
!
!! P19_AUTONOMOUS reports whether problem 19 is autonomous.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P19_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p19_autonomous

  p19_autonomous = .true.

  return
end
subroutine p19_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P19_EQUIL returns equilibrium solutions of problem 19.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p19_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P19_FUN evaluates the function for problem 19.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(3)
  yp(2) = y(4)
  yp(3) = - y(1) / ( sqrt ( ( y(1)**2 + y(2)**2 ) ) )**3
  yp(4) = - y(2) / ( sqrt ( ( y(1)**2 + y(2)**2 ) ) )**3

  return
end
subroutine p19_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P19_JAC evaluates the jacobian for problem 19.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1:neqn,1:neqn) = 0.0D+00

  return
end
subroutine p19_neqn ( neqn )

!*****************************************************************************80
!
!! P19_NEQN returns the number of equations for problem 19.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 4

  return
end
subroutine p19_param ( action, name, value )

!*****************************************************************************80
!
!! P19_PARAM handles the parameters for problem 19.
!
!  Modified:
!
!    20 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'DELTA'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: delta = 0.7D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      value = delta
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P19_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      delta = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P19_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P19_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p19_scale ( neqn, scale )

!*****************************************************************************80
!
!! P19_SCALE returns scale factors for problem 19.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p19_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P19_START returns the starting point for problem 19.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) delta
  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  call p19_param ( 'GET', 'DELTA', delta )

  y_start(1) = 1.0D+00 - delta
  y_start(2) = 0.0D+00
  y_start(3) = 0.0D+00
  y_start(4) = sqrt ( ( 1.0D+00 + delta ) / ( 1.0D+00 - delta ) )

  return
end
subroutine p19_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P19_STOP returns the stopping point for problem 19.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ &
   -9.538990293416394D-01, &
    6.907409024219432D-01, &
   -8.212674270877433D-01, &
   -1.539574259125825D-01 /)

  return
end
subroutine p19_title ( title )

!*****************************************************************************80
!
!! P19_TITLE returns the title of problem 19.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 19, Enright and Pryce #D4'

  return
end
function p20_autonomous ( )

!*****************************************************************************80
!
!! P20_AUTONOMOUS reports whether problem 20 is autonomous.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P20_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p20_autonomous

  p20_autonomous = .true.

  return
end
subroutine p20_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P20_EQUIL returns equilibrium solutions of problem 20.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p20_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P20_FUN evaluates the function for problem 20.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(3)
  yp(2) = y(4)
  yp(3) = - y(1) / ( sqrt ( ( y(1)**2 + y(2)**2 ) ) )**3
  yp(4) = - y(2) / ( sqrt ( ( y(1)**2 + y(2)**2 ) ) )**3

  return
end
subroutine p20_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P20_JAC evaluates the jacobian for problem 20.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1:neqn,1:neqn) = 0.0D+00

  return
end
subroutine p20_neqn ( neqn )

!*****************************************************************************80
!
!! P20_NEQN returns the number of equations for problem 20.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 4

  return
end
subroutine p20_param ( action, name, value )

!*****************************************************************************80
!
!! P20_PARAM handles the parameters for problem 20.
!
!  Modified:
!
!    20 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'DELTA'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: delta = 0.9D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      value = delta
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      delta = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P20_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p20_scale ( neqn, scale )

!*****************************************************************************80
!
!! P20_SCALE returns scale factors for problem 20.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p20_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P20_START returns the starting point for problem 20.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) delta
  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  delta = 0.9D+00

  y_start(1) = 1.0D+00 - delta
  y_start(2) = 0.0D+00
  y_start(3) = 0.0D+00
  y_start(4) = sqrt ( ( 1.0D+00 + delta ) / ( 1.0D+00 - delta ) )

  return
end
subroutine p20_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P20_STOP returns the stopping point for problem 20.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) delta
  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ &
   -1.295266250987574D+00, &
    4.003938963792321D-01, &
   -6.775390924707566D-01, &
   -1.270838154278686D-01 /)

  return
end
subroutine p20_title ( title )

!*****************************************************************************80
!
!! P20_TITLE returns the title of problem 20.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 20, Enright and Pryce #D5'

  return
end
function p21_autonomous ( )

!*****************************************************************************80
!
!! P21_AUTONOMOUS reports whether problem 21 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P21_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p21_autonomous

  p21_autonomous = .true.

  return
end
subroutine p21_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P21_EQUIL returns equilibrium solutions of problem 21.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p21_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P21_FUN evaluates the function for problem 21.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(2)
  yp(2) = - ( 1.0D+00 - 0.25D+00 / ( t + 1.0D+00 )**2 ) * y(1) &
    - 1.0D+00 / ( t + 1.0D+00 ) * y(2)

  return
end
subroutine p21_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P21_JAC evaluates the jacobian for problem 21.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = 0.0D+00
  jac(1,2) = 1.0D+00

  jac(2,1) = 1.0 - 0.25D+00 / ( t + 1.0D+00 )**2
  jac(2,2) = - 1.0D+00 / ( t + 1.0D+00 )

  return
end
subroutine p21_neqn ( neqn )

!*****************************************************************************80
!
!! P21_NEQN returns the number of equations for problem 21.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p21_scale ( neqn, scale )

!*****************************************************************************80
!
!! P21_SCALE returns scale factors for problem 21.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p21_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P21_START returns the starting point for problem 21.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn) = (/ 0.6713967071418030D+00, 0.09540051444747446D+00 /)

  return
end
subroutine p21_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P21_STOP returns the stopping point for problem 21.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ 1.456723600728308D-01, -9.883500195574063D-02 /)

  return
end
subroutine p21_title ( title )

!*****************************************************************************80
!
!! P21_TITLE returns the title of problem 21.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 21, Enright and Pryce #E1'

  return
end
function p22_autonomous ( )

!*****************************************************************************80
!
!! P22_AUTONOMOUS reports whether problem 22 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P22_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p22_autonomous

  p22_autonomous = .true.

  return
end
subroutine p22_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P22_EQUIL returns equilibrium solutions of problem 22.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p22_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P22_FUN evaluates the function for problem 22.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(2)
  yp(2) = ( 1.0D+00 - y(1)**2 ) * y(2) - y(1)

  return
end
subroutine p22_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P22_JAC evaluates the jacobian for problem 22.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = 0.0D+00
  jac(1,2) = 1.0D+00

  jac(2,1) = - 2.0D+00 * y(1) * y(2) - 1.0D+00
  jac(2,2) = 1.0D+00 - y(1)**2

  return
end
subroutine p22_neqn ( neqn )

!*****************************************************************************80
!
!! P22_NEQN returns the number of equations for problem 22.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p22_scale ( neqn, scale )

!*****************************************************************************80
!
!! P22_SCALE returns scale factors for problem 22.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p22_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P22_START returns the starting point for problem 22.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn) = (/ 2.0D+00, 0.0D+00 /)

  return
end
subroutine p22_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P22_STOP returns the stopping point for problem 22.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ 2.008149762174948D+00, -4.250887527320057D-02 /)

  return
end
subroutine p22_title ( title )

!*****************************************************************************80
!
!! P22_TITLE returns the title of problem 22.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 22, Enright and Pryce #E2'

  return
end
function p23_autonomous ( )

!*****************************************************************************80
!
!! P23_AUTONOMOUS reports whether problem 23 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P23_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p23_autonomous

  p23_autonomous = .true.

  return
end
subroutine p23_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P23_EQUIL returns equilibrium solutions of problem 23.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none
!
  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p23_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P23_FUN evaluates the function for problem 23.
!
!  Modified:
!
!    06 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(2)
  yp(2) = y(1)**3 / 6.0D+00 - y(1) + 2.0D+00 * sin ( 2.78535D+00 * t )

  return
end
subroutine p23_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P23_JAC evaluates the jacobian for problem 23.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = 0.0D+00
  jac(1,2) = 1.0D+00

  jac(2,1) = 0.5D+00 * y(1)**2 - 1.0D+00
  jac(2,2) = 0.0D+00

  return
end
subroutine p23_neqn ( neqn )

!*****************************************************************************80
!
!! P23_NEQN returns the number of equations for problem 23.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p23_scale ( neqn, scale )

!*****************************************************************************80
!
!! P23_SCALE returns scale factors for problem 23.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p23_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P23_START returns the starting point for problem 23.
!
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn) = (/ 0.0D+00, 0.0D+00 /)

  return
end
subroutine p23_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P23_STOP returns the stopping point for problem 23.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ -1.004178858647128D-01, 2.411400132095954D-01 /)

  return
end
subroutine p23_title ( title )

!*****************************************************************************80
!
!! P23_TITLE returns the title of problem 23.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 23, Enright and Pryce #E3'

  return
end
function p24_autonomous ( )

!*****************************************************************************80
!
!! P24_AUTONOMOUS reports whether problem 24 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P24_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p24_autonomous

  p24_autonomous = .true.

  return
end
subroutine p24_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P24_EQUIL returns equilibrium solutions of problem 24.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p24_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P24_FUN evaluates the function for problem 24.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(2)
  yp(2) = 0.032D+00 - 0.4D+00 * y(2)**2

  return
end
subroutine p24_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P24_JAC evaluates the jacobian for problem 24.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = 0.0D+00
  jac(1,2) = 1.0D+00

  jac(2,1) = 0.0D+00
  jac(2,2) = -0.8D+00 * y(2)

  return
end
subroutine p24_neqn ( neqn )

!*****************************************************************************80
!
!! P24_NEQN returns the number of equations for problem 24.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p24_scale ( neqn, scale )

!*****************************************************************************80
!
!! P24_SCALE returns scale factors for problem 24.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p24_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P24_START returns the starting point for problem 24.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn) = (/ 30.0D+00, 0.0D+00 /)

  return
end
subroutine p24_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P24_STOP returns the stopping point for problem 24.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ 3.395091444646555D+01, 2.767822659672869D-01 /)

  return
end
subroutine p24_title ( title )

!*****************************************************************************80
!
!! P24_TITLE returns the title of problem 24.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 24, Enright and Pryce #E4'

  return
end
function p25_autonomous ( )

!*****************************************************************************80
!
!! P25_AUTONOMOUS reports whether problem 25 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P25_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p25_autonomous

  p25_autonomous = .false.

  return
end
subroutine p25_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P25_EQUIL returns equilibrium solutions of problem 25.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p25_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P25_FUN evaluates the function for problem 25.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(2)
  yp(2) = sqrt ( 1.0D+00 + y(2)**2 ) / ( 25.0D+00 - t )

  return
end
subroutine p25_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P25_JAC evaluates the jacobian for problem 25.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = 0.0D+00
  jac(1,2) = 1.0D+00

  jac(2,1) = 0.0D+00
  jac(2,2) = y(2) / ( sqrt ( 1.0D+00 + y(2)**2 ) * ( 25.0D+00 - t ) )

  return
end
subroutine p25_neqn ( neqn )

!*****************************************************************************80
!
!! P25_NEQN returns the number of equations for problem 25.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p25_scale ( neqn, scale )

!*****************************************************************************80
!
!! P25_SCALE returns scale factors for problem 25.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p25_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P25_START returns the starting point for problem 25.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn) = (/ 0.0D+00, 0.0D+00 /)

  return
end
subroutine p25_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P25_STOP returns the stopping point for problem 25.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ 1.411797390542629D+01, 2.400000000000002D+00 /)

  return
end
subroutine p25_title ( title )

!*****************************************************************************80
!
!! P25_TITLE returns the title of problem 25.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 25, Enright and Pryce #E5'

  return
end
function p26_autonomous ( )

!*****************************************************************************80
!
!! P26_AUTONOMOUS reports whether problem 26 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P26_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p26_autonomous

  p26_autonomous = .true.

  return
end
subroutine p26_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P26_EQUIL returns equilibrium solutions of problem 26.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p26_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P26_FUN evaluates the function for problem 26.
!
!  Modified:
!
!    19 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) a
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  call p26_param ( 'GET', 'A', a )

  yp(1) = y(2)

  if ( mod ( int ( t ), 2 ) == 0 ) then
    yp(2) = 2.0D+00 * a * y(2) - ( pi**2 + a**2 ) * y(1) + 1.0D+00
  else
    yp(2) = 2.0D+00 * a * y(2) - ( pi**2 + a**2 ) * y(1) - 1.0D+00
  end if

  return
end
subroutine p26_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P26_JAC evaluates the jacobian for problem 26.
!
!  Modified:
!
!    19 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) a
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  call p26_param ( 'GET', 'A', a )

  jac(1,1) = 0.0D+00
  jac(1,2) = 1.0D+00

  jac(2,1) = - ( pi**2 + a**2 )
  jac(2,2) = 2.0D+00 * a

  return
end
subroutine p26_neqn ( neqn )

!*****************************************************************************80
!
!! P26_NEQN returns the number of equations for problem 26.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p26_param ( action, name, value )

!*****************************************************************************80
!
!! P26_PARAM handles the parameters for problem 26.
!
!  Discussion:
!
!    A is the only parameter, and has a default value of 0.1.
!
!  Modified:
!
!    19 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'A'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: a = 0.1D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'A' ) ) then
      value = a
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P26_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'A' ) ) then
      a = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P26_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P26_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p26_scale ( neqn, scale )

!*****************************************************************************80
!
!! P26_SCALE returns scale factors for problem 26.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p26_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P26_START returns the starting point for problem 26.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn) = (/ 0.0D+00, 0.0D+00 /)

  return
end
subroutine p26_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P26_STOP returns the stopping point for problem 26.
!
!  Discussion:
!
!    The value of Y_STOP is only valid for the default parameter
!    value A = 0.1.
!
!  Modified:
!
!    19 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ -1.294460621213470D+01, -2.208575158908672D-15 /)

  return
end
subroutine p26_title ( title )

!*****************************************************************************80
!
!! P26_TITLE returns the title of problem 26.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 26, Enright and Pryce #F1'

  return
end
function p27_autonomous ( )

!*****************************************************************************80
!
!! P27_AUTONOMOUS reports whether problem 27 is autonomous.
!
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P27_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p27_autonomous

  p27_autonomous = .false.

  return
end
subroutine p27_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P27_EQUIL returns equilibrium solutions of problem 27.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next

  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p27_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P27_FUN evaluates the function for problem 27.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  if ( mod ( int ( t ), 2 ) == 0 ) then
    yp(1) = 55.0D+00 - 1.5D+00 * y(1)
  else
    yp(1) = 55.0D+00 - 0.5D+00 * y(1)
  end if

  return
end
subroutine p27_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P27_JAC evaluates the jacobian for problem 27.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  if ( mod ( int ( t ), 2 ) == 0 ) then
    jac(1,1) = - 1.5D+00
  else
    jac(1,1) = - 0.5D+00
  end if

  return
end
subroutine p27_neqn ( neqn )

!*****************************************************************************80
!
!! P27_NEQN returns the number of equations for problem 27.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 1

  return
end
subroutine p27_scale ( neqn, scale )

!*****************************************************************************80
!
!! P27_SCALE returns scale factors for problem 27.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p27_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P27_START returns the starting point for problem 27.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn) = (/ 110.0D+00 /)

  return
end
subroutine p27_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P27_STOP returns the stopping point for problem 27.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ 70.03731057008607D+00 /)

  return
end
subroutine p27_title ( title )

!*****************************************************************************80
!
!! P27_TITLE returns the title of problem 27.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 27, Enright and Pryce #F2'

  return
end
function p28_autonomous ( )

!*****************************************************************************80
!
!! P28_AUTONOMOUS reports whether problem 14 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P28_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p28_autonomous

  p28_autonomous = .true.

  return
end
subroutine p28_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P28_EQUIL returns equilibrium solutions of problem 28.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p28_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P28_FUN evaluates the function for problem 28.
!
!  Modified:
!
!    19 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(2)

  yp(2) = 0.01D+00 * y(2) * ( 1.0D+00 - y(1)**2 ) - y(1) &
    - abs ( sin ( pi * t ) )

  return
end
subroutine p28_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P28_JAC evaluates the jacobian for problem 28.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = 0.0D+00
  jac(1,2) = 1.0D+00

  jac(2,1) = - 0.02D+00 * y(2) * y(1) - 1.0D+00
  jac(2,2) = 0.01D+00 * ( 1.0D+00 - y(1)**2 )

  return
end
subroutine p28_neqn ( neqn )

!*****************************************************************************80
!
!! P28_NEQN returns the number of equations for problem 28.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p28_scale ( neqn, scale )

!*****************************************************************************80
!
!! P28_SCALE returns scale factors for problem 28.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p28_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P28_START returns the starting point for problem 28.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn) = (/ 0.0D+00, 0.0D+00 /)

  return
end
subroutine p28_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P28_STOP returns the stopping point for problem 28.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ -3.726957553088175D-01, -6.230137949234190D-01 /)

  return
end
subroutine p28_title ( title )

!*****************************************************************************80
!
!! P28_TITLE returns the title of problem 28.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 28, Enright and Pryce #F3'

  return
end
function p29_autonomous ( )

!*****************************************************************************80
!
!! P29_AUTONOMOUS reports whether problem 29 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P29_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p29_autonomous

  p29_autonomous = .false.

  return
end
subroutine p29_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P29_EQUIL returns equilibrium solutions of problem 29.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p29_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P29_FUN evaluates the function for problem 29.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  if ( t <= 10.0D+00 ) then
    yp(1) = - 2.0D+00 / 21.0D+00 &
      - 120.0D+00 * ( t - 5.0D+00 ) / ( 1.0D+00 + 4.0D+00 * ( t - 5.0D+00 )**2 )
  else
    yp(1) = - 2.0D+00 * y(1)
  end if

  return
end
subroutine p29_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P29_JAC evaluates the jacobian for problem 29.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  if ( t <= 10.0D+00 ) then
    jac(1,1) = 0.0D+00
  else
    jac(1,1) = -2.0D+00
  end if

  return
end
subroutine p29_neqn ( neqn )

!*****************************************************************************80
!
!! P29_NEQN returns the number of equations for problem 29.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 1

  return
end
subroutine p29_scale ( neqn, scale )

!*****************************************************************************80
!
!! P29_SCALE returns scale factors for problem 29.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p29_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P29_START returns the starting point for problem 29.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn) = 1.0D+00

  return
end
subroutine p29_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P29_STOP returns the stopping point for problem 29.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ 9.815017249707434D-11 /)

  return
end
subroutine p29_title ( title )

!*****************************************************************************80
!
!! P29_TITLE returns the title of problem 29.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 29, Enright and Pryce #F4'

  return
end
function p30_autonomous ( )

!*****************************************************************************80
!
!! P30_AUTONOMOUS reports whether problem 30 is autonomous.
!
!  Modified:
!
!    03 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P30_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p30_autonomous

  p30_autonomous = .false.

  return
end
subroutine p30_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P30_EQUIL returns equilibrium solutions of problem 30.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none
  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:neqn) = 0.0D+00
  else
    next = 0
  end if

  return
end
subroutine p30_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P30_FUN evaluates the function for problem 30.
!
!  Modified:
!
!    06 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wayne Enright, John Pryce,
!    Algorithm 648,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 1, pages 28-34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) csum
  integer ( kind = 4 ) i
  real ( kind = 8 ) pprime
  real ( kind = 8 ) r8_cube_root
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  csum = 0.0D+00
  do i = 1, 19
    csum = csum + real ( i, kind = 8 )**(4.0D+00/3.0D+00)
  end do

  pprime = 0.0D+00
  do i = 1, 19
    pprime = pprime + ( 4.0D+00 / 3.0D+00 ) &
      * r8_cube_root ( t - real ( i, kind = 8 ) )
  end do

  yp(1) = pprime * y(1) / csum

  return
end
subroutine p30_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P30_JAC evaluates the jacobian for problem 30.
!
!  Modified:
!
!    06 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) csum
  integer ( kind = 4 ) i
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) pprime
  real ( kind = 8 ) r8_cube_root
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  csum = 0.0D+00
  do i = 1, 19
    csum = csum + real ( i, kind = 8 )**(4.0D+00/3.0D+00)
  end do

  pprime = 0.0D+00
  do i = 1, 19
    pprime = pprime + ( 4.0D+00 / 3.0D+00 ) &
      * r8_cube_root ( t - real ( i, kind = 8 ) )
  end do

  jac(1,1) = pprime / csum

  return
end
subroutine p30_neqn ( neqn )

!*****************************************************************************80
!
!! P30_NEQN returns the number of equations for problem 30.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 1

  return
end
subroutine p30_scale ( neqn, scale )

!*****************************************************************************80
!
!! P30_SCALE returns scale factors for problem 30.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:neqn) = 1.0D+00

  return
end
subroutine p30_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P30_START returns the starting point for problem 30.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00

  y_start(1:neqn) = 1.0D+00

  return
end
subroutine p30_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P30_STOP returns the stopping point for problem 30.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:neqn) = (/ 1.0D+00 /)

  return
end
subroutine p30_title ( title )

!*****************************************************************************80
!
!! P30_TITLE returns the title of problem 30.
!
!  Modified:
!
!    04 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 30, Enright and Pryce #F5'

  return
end
function p31_autonomous ( )

!*****************************************************************************80
!
!! P31_AUTONOMOUS reports whether problem 31 is autonomous.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P31_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p31_autonomous

  p31_autonomous = .true.

  return
end
subroutine p31_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P31_EQUIL returns equilibrium solutions of problem 31.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  call p31_param ( 'GET', 'A', a )
  call p31_param ( 'GET', 'B', b )
  call p31_param ( 'GET', 'C', c )
  call p31_param ( 'GET', 'D', d )

  if ( next == 0 ) then
    next = 1
    y(1:2) = (/ 0.0D+00, 0.0D+00 /)
  else if ( next == 1 .and. c /= 0.0D+00 .and. d /= 0.0D+00 ) then
    next = 2
    y(1:2) = (/ d / c, a / b /)
  else
    next = 0
  end if

  return
end
subroutine p31_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P31_FUN evaluates the function for problem 31.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  call p31_param ( 'GET', 'A', a )
  call p31_param ( 'GET', 'B', b )
  call p31_param ( 'GET', 'C', c )
  call p31_param ( 'GET', 'D', d )

  yp(1) = ( a - b * y(2) ) * y(1)
  yp(2) = ( c * y(1) - d ) * y(2)

  return
end
subroutine p31_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P31_JAC evaluates the jacobian for problem 31.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  call p31_param ( 'GET', 'A', a )
  call p31_param ( 'GET', 'B', b )
  call p31_param ( 'GET', 'C', c )
  call p31_param ( 'GET', 'D', d )

  jac(1,1) =  a - b * y(2)
  jac(1,2) =    - b * y(1)

  jac(2,1) = c * y(2)
  jac(2,2) = c * y(1) - d

  return
end
subroutine p31_neqn ( neqn )

!*****************************************************************************80
!
!! P31_NEQN returns the number of equations for problem 31.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p31_param ( action, name, value )

!*****************************************************************************80
!
!! P31_PARAM handles the parameters for problem 31.
!
!  Discussion:
!
!    A, B, C, and D are the positive coefficients in the equations:
!
!    Y1' = ( A - B * Y2 ) * Y1
!    Y2' = ( C * Y1 - D ) * Y2
!
!  Modified:
!
!    19 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'A', 'B', 'C' or 'D'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: a = 5.0D+00
  real ( kind = 8 ), save :: b = 1.0D+00
  real ( kind = 8 ), save :: c = 0.5D+00
  real ( kind = 8 ), save :: d = 2.0D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'A' ) ) then
      value = a
    else if ( s_eqi ( name, 'B' ) ) then
      value = b
    else if ( s_eqi ( name, 'C' ) ) then
      value = c
    else if ( s_eqi ( name, 'D' ) ) then
      value = d
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P31_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'A' ) ) then
      a = value
    else if ( s_eqi ( name, 'B' ) ) then
      b = value
    else if ( s_eqi ( name, 'C' ) ) then
      c = value
    else if ( s_eqi ( name, 'D' ) ) then
      d = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P31_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P31_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p31_scale ( neqn, scale )

!*****************************************************************************80
!
!! P31_SCALE returns scale factors for problem 31.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) scale(neqn)

  call p31_param ( 'GET', 'A', a )
  call p31_param ( 'GET', 'B', b )
  call p31_param ( 'GET', 'C', c )
  call p31_param ( 'GET', 'D', d )

  if ( c /= 0.0D+00 .and. d /= 0.0D+00 ) then
    scale(1:2) = (/ d / c, a / b /)
  else
    scale(1:2) = (/ 1.0D+00, 1.0D+00 /)
  end if

  return
end
subroutine p31_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P31_START returns the starting point for problem 31.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = (/ 2.0D+00, 2.0D+00 /)

  return
end
subroutine p31_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P31_STOP returns the stopping point for problem 31.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 10.0D+00

  y_stop(1:2) = (/ 2.20050D+00, 10.2726D+00 /)

  return
end
subroutine p31_title ( title )

!*****************************************************************************80
!
!! P31_TITLE returns the title of problem 31.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 31, Lotka-Volterra Predator-Prey Equations.'

  return
end
function p32_autonomous ( )

!*****************************************************************************80
!
!! P32_AUTONOMOUS reports whether problem 32 is autonomous.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P32_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p32_autonomous

  p32_autonomous = .true.

  return
end
subroutine p32_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P32_EQUIL returns equilibrium solutions of problem 32.
!
!  Modified:
!
!    19 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) beta
  integer ( kind = 4 ) next
  real ( kind = 8 ) rho
  real ( kind = 8 ) s
  real ( kind = 8 ) sigma
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  call p32_param ( 'GET', 'BETA', beta )
  call p32_param ( 'GET', 'RHO', rho )
  call p32_param ( 'GET', 'SIGMA', sigma )

  if ( next == 0 ) then
    next = 1
    y(1:3) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  else if ( rho <= 1.0D+00 ) then
    next = 0
  else if ( next == 1 ) then
    next = 2
    t = rho - 1.0D+00
    s = sqrt ( beta * ( rho - 1.0D+00 ) )
    y(1:3) = (/ s, s, t /)
  else if ( next == 2 ) then
    next = 3
    t = rho - 1.0D+00
    s = sqrt ( beta * ( rho - 1.0D+00 ) )
    y(1:3) = (/ -s, -s, t /)
  else
    next = 0
  end if

  return
end
subroutine p32_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P32_FUN evaluates the function for problem 32.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) beta
  real ( kind = 8 ) rho
  real ( kind = 8 ) sigma
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  call p32_param ( 'GET', 'BETA', beta )
  call p32_param ( 'GET', 'RHO', rho )
  call p32_param ( 'GET', 'SIGMA', sigma )

  yp(1) = sigma * ( y(2) - y(1) )
  yp(2) = rho * y(1) - y(2) - y(1) * y(3)
  yp(3) = - beta * y(3) + y(1) * y(2)

  return
end
subroutine p32_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P32_JAC evaluates the jacobian for problem 32.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) beta
  real ( kind = 8 ) rho
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) sigma
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  call p32_param ( 'GET', 'BETA', beta )
  call p32_param ( 'GET', 'RHO', rho )
  call p32_param ( 'GET', 'SIGMA', sigma )

  jac(1,1) = - sigma
  jac(1,2) =   sigma
  jac(1,3) =   0.0D+00

  jac(2,1) =   rho - y(3)
  jac(2,2) = - 1.0D+00
  jac(2,3) = - y(1)

  jac(3,1) =   y(2)
  jac(3,2) =   y(1)
  jac(3,3) =   - beta

  return
end
subroutine p32_neqn ( neqn )

!*****************************************************************************80
!
!! P32_NEQN returns the number of equations for problem 32.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 3

  return
end
subroutine p32_param ( action, name, value )

!*****************************************************************************80
!
!! P32_PARAM handles the parameters for problem 32.
!
!  Discussion:
!
!    BETA, RHO and SIGMA should be positive.
!
!    For "interesting" behavior, RHO should be greater than 1.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'BETA', 'RHO', or 'SIGMA'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: beta = 8.0D+00 / 3.0D+00
  character ( len = * ) name
  real ( kind = 8 ), save :: rho = 28.0D+00
  logical s_eqi
  real ( kind = 8 ), save :: sigma = 10.0D+00
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'BETA' ) ) then
      value = beta
    else if ( s_eqi ( name, 'RHO' ) ) then
      value = rho
    else if ( s_eqi ( name, 'SIGMA' ) ) then
      value = sigma
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P32_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'BETA' ) ) then
      beta = value
    else if ( s_eqi ( name, 'RHO' ) ) then
      rho = value
    else if ( s_eqi ( name, 'SIGMA' ) ) then
      sigma = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P32_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P32_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p32_scale ( neqn, scale )

!*****************************************************************************80
!
!! P32_SCALE returns scale factors for problem 32.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:3) = (/ 1.0D+00, 1.0D+00, 1.0D+00 /)

  return
end
subroutine p32_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P32_START returns the starting point for problem 32.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = (/ 2.0D+00, 2.0D+00, 21.0D+00 /)

  return
end
subroutine p32_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P32_STOP returns the stopping point for problem 32.
!
!  Discussion:
!
!    The system is chaotic, and so a dummy stop value is put here.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:3) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)

  return
end
subroutine p32_title ( title )

!*****************************************************************************80
!
!! P32_TITLE returns the title of problem 32.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 32, The Lorenz System'

  return
end
function p33_autonomous ( )

!*****************************************************************************80
!
!! P33_AUTONOMOUS reports whether problem 33 is autonomous.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P33_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p33_autonomous

  p33_autonomous = .true.

  return
end
subroutine p33_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P33_EQUIL returns equilibrium solutions of problem 33.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:2) = (/ 0.0D+00, 0.0D+00 /)
  else
    next = 0
  end if

  return
end
subroutine p33_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P33_FUN evaluates the function for problem 33.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) delta
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  call p33_param ( 'GET', 'DELTA', delta )

  yp(1) = y(2)
  yp(2) = delta * ( 1.0D+00 - y(1)**2 ) * y(2) - y(1)

  return
end
subroutine p33_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P33_JAC evaluates the jacobian for problem 33.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) delta
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  call p33_param ( 'GET', 'DELTA', delta )

  jac(1,1) =  0.0D+00
  jac(1,2) =  1.0D+00

  jac(2,1) = - 1.0D+00 - 2.0D+00 * delta * y(1) * y(2)
  jac(2,2) = delta * ( 1.0D+00 - y(1)**2 )

  return
end
subroutine p33_neqn ( neqn )

!*****************************************************************************80
!
!! P33_NEQN returns the number of equations for problem 33.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p33_param ( action, name, value )

!*****************************************************************************80
!
!! P33_PARAM handles the parameters for problem 33.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'DELTA'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: delta = 1.0D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      value = delta
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P33_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      delta = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P33_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P33_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p33_scale ( neqn, scale )

!*****************************************************************************80
!
!! P33_SCALE returns scale factors for problem 33.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:2) = (/ 1.0D+00, 1.0D+00 /)

  return
end
subroutine p33_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P33_START returns the starting point for problem 33.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = (/ 2.0D+00, 2.0D+00 /)

  return
end
subroutine p33_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P33_STOP returns the stopping point for problem 33.
!
!  Modified:
!
!    21 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:2) = (/ 0.756245D+00, 2.67294D+00 /)

  return
end
subroutine p33_title ( title )

!*****************************************************************************80
!
!! P33_TITLE returns the title of problem 33.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 33, The Van der Pol equation'

  return
end
function p34_autonomous ( )

!*****************************************************************************80
!
!! P34_AUTONOMOUS reports whether problem 34 is autonomous.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P34_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p34_autonomous

  p34_autonomous = .true.

  return
end
subroutine p34_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P34_EQUIL returns equilibrium solutions of problem 34.
!
!  Modified:
!
!    19 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:2) = (/ 0.0D+00, 0.0D+00 /)
  else if ( next == 1 ) then
    next = 2
    y(1:2) = (/ pi, 0.0D+00 /)
  else
    next = 0
  end if

  return
end
subroutine p34_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P34_FUN evaluates the function for problem 34.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) g
  real ( kind = 8 ) k
  real ( kind = 8 ) l
  real ( kind = 8 ) m
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  call p34_param ( 'GET', 'G', g )
  call p34_param ( 'GET', 'K', k )
  call p34_param ( 'GET', 'L', l )
  call p34_param ( 'GET', 'M', m )

  yp(1) = y(2)
  yp(2) = - ( g / l ) * y(1) - ( k / m ) * y(2)

  return
end
subroutine p34_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P34_JAC evaluates the jacobian for problem 34.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) g
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) k
  real ( kind = 8 ) l
  real ( kind = 8 ) m
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  call p34_param ( 'GET', 'G', g )
  call p34_param ( 'GET', 'K', k )
  call p34_param ( 'GET', 'L', l )
  call p34_param ( 'GET', 'M', m )

  jac(1,1) =  0.0D+00
  jac(1,2) =  1.0D+00

  jac(2,1) = - ( g / l )
  jac(2,2) = - ( k / m )

  return
end
subroutine p34_neqn ( neqn )

!*****************************************************************************80
!
!! P34_NEQN returns the number of equations for problem 34.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p34_param ( action, name, value )

!*****************************************************************************80
!
!! P34_PARAM handles the parameters for problem 34.
!
!  Discussion:
!
!    G is the gravitational force;
!    K is the damping coefficient;
!    L is the length of the string;
!    M is the mass of the pendulum bob.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'G', 'K', 'L' or 'M'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: g = 32.0D+00
  real ( kind = 8 ), save :: k = 1.0D+00
  real ( kind = 8 ), save :: l = 1.0D+00
  real ( kind = 8 ), save :: m = 1.0D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'G' ) ) then
      value = g
    else if ( s_eqi ( name, 'K' ) ) then
      value = k
    else if ( s_eqi ( name, 'L' ) ) then
      value = l
    else if ( s_eqi ( name, 'M' ) ) then
      value = m
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P34_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'G' ) ) then
      g = value
    else if ( s_eqi ( name, 'K' ) ) then
      k = value
    else if ( s_eqi ( name, 'L' ) ) then
      l = value
    else if ( s_eqi ( name, 'M' ) ) then
      m = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P34_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P34_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p34_scale ( neqn, scale )

!*****************************************************************************80
!
!! P34_SCALE returns scale factors for problem 34.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:2) = (/ 1.0D+00, 1.0D+00 /)

  return
end
subroutine p34_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P34_START returns the starting point for problem 34.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = (/ 2.0D+00, 2.0D+00 /)

  return
end
subroutine p34_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P34_STOP returns the stopping point for problem 34.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:2) = (/ 0.695786D-04, 0.277616D-03 /)

  return
end
subroutine p34_title ( title )

!*****************************************************************************80
!
!! P34_TITLE returns the title of problem 34.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 34, The Linearized Damped Pendulum'

  return
end
function p35_autonomous ( )

!*****************************************************************************80
!
!! P35_AUTONOMOUS reports whether problem 35 is autonomous.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P35_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p35_autonomous

  p35_autonomous = .true.

  return
end
subroutine p35_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P35_EQUIL returns equilibrium solutions of problem 35.
!
!  Modified:
!
!    24 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:2) = (/ 0.0D+00, 0.0D+00 /)
  else
    next = 0
  end if

  return
end
subroutine p35_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P35_FUN evaluates the function for problem 35.
!
!  Modified:
!
!    24 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) g
  real ( kind = 8 ) k
  real ( kind = 8 ) l
  real ( kind = 8 ) m
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  call p35_param ( 'GET', 'G', g )
  call p35_param ( 'GET', 'K', k )
  call p35_param ( 'GET', 'L', l )
  call p35_param ( 'GET', 'M', m )

  yp(1) = y(2)
  yp(2) = - ( g / l ) * sin ( y(1) ) - ( k / m ) * y(2)

  return
end
subroutine p35_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P35_JAC evaluates the jacobian for problem 35.
!
!  Modified:
!
!    24 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) g
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) k
  real ( kind = 8 ) l
  real ( kind = 8 ) m
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  call p35_param ( 'GET', 'G', g )
  call p35_param ( 'GET', 'K', k )
  call p35_param ( 'GET', 'L', l )
  call p35_param ( 'GET', 'M', m )

  jac(1,1) =  0.0D+00
  jac(1,2) =  1.0D+00

  jac(2,1) = - ( g / l ) * cos ( y(1) )
  jac(2,2) = - ( k / m )

  return
end
subroutine p35_neqn ( neqn )

!*****************************************************************************80
!
!! P35_NEQN returns the number of equations for problem 35.
!
!  Modified:
!
!    24 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p35_param ( action, name, value )

!*****************************************************************************80
!
!! P35_PARAM handles the parameters for problem 35.
!
!  Discussion:
!
!    G is the gravitational force;
!    K is the damping coefficient;
!    L is the length of the string;
!    M is the mass of the pendulum bob.
!
!  Modified:
!
!    24 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'G', 'K', 'L' or 'M'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: g = 32.0D+00
  real ( kind = 8 ), save :: k = 1.0D+00
  real ( kind = 8 ), save :: l = 1.0D+00
  real ( kind = 8 ), save :: m = 1.0D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'G' ) ) then
      value = g
    else if ( s_eqi ( name, 'K' ) ) then
      value = k
    else if ( s_eqi ( name, 'L' ) ) then
      value = l
    else if ( s_eqi ( name, 'M' ) ) then
      value = m
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P35_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'G' ) ) then
      g = value
    else if ( s_eqi ( name, 'K' ) ) then
      k = value
    else if ( s_eqi ( name, 'L' ) ) then
      l = value
    else if ( s_eqi ( name, 'M' ) ) then
      m = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P35_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P35_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p35_scale ( neqn, scale )

!*****************************************************************************80
!
!! P35_SCALE returns scale factors for problem 35.
!
!  Modified:
!
!    24 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:2) = (/ 1.0D+00, 1.0D+00 /)

  return
end
subroutine p35_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P35_START returns the starting point for problem 35.
!
!  Modified:
!
!    24 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = (/ 2.0D+00, 2.0D+00 /)

  return
end
subroutine p35_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P35_STOP returns the stopping point for problem 35.
!
!  Modified:
!
!    24 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 20.0D+00

  y_stop(1:2) = (/ -0.584253D-04, 0.359969D-03 /)

  return
end
subroutine p35_title ( title )

!*****************************************************************************80
!
!! P35_TITLE returns the title of problem 35.
!
!  Modified:
!
!    24 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 35, The Nonlinear Damped Pendulum'

  return
end
function p36_autonomous ( )

!*****************************************************************************80
!
!! P36_AUTONOMOUS reports whether problem 36 is autonomous.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P36_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p36_autonomous

  p36_autonomous = .true.

  return
end
subroutine p36_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P36_EQUIL returns equilibrium solutions of problem 36.
!
!  Modified:
!
!    26 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:2) = (/ 0.0D+00, 0.0D+00 /)
  else if ( next == 1 ) then
    next = 2
    y(1:2) = (/ 1.0D+00, 0.0D+00 /)
  else if ( next == 2 ) then
    next = 3
    y(1:2) = (/ -1.0D+00, 0.0D+00 /)
  else
    next = 0
  end if

  return
end
subroutine p36_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P36_FUN evaluates the function for problem 36.
!
!  Modified:
!
!    26 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(2)
  yp(2) = y(1) * ( 1.0D+00 - y(1)**2 )

  return
end
subroutine p36_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P36_JAC evaluates the jacobian for problem 36.
!
!  Modified:
!
!    26 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) =  0.0D+00
  jac(1,2) =  1.0D+00

  jac(2,1) = 1.0D+00 - 3.0D+00 * y(1)**2
  jac(2,2) = 0.0D+00

  return
end
subroutine p36_neqn ( neqn )

!*****************************************************************************80
!
!! P36_NEQN returns the number of equations for problem 36.
!
!  Modified:
!
!    26 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p36_scale ( neqn, scale )

!*****************************************************************************80
!
!! P36_SCALE returns scale factors for problem 36.
!
!  Modified:
!
!    26 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:2) = (/ 1.0D+00, 1.0D+00 /)

  return
end
subroutine p36_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P36_START returns the starting point for problem 36.
!
!  Modified:
!
!    26 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = (/ 0.5D+00, 0.0D+00 /)

  return
end
subroutine p36_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P36_STOP returns the stopping point for problem 36.
!
!  Modified:
!
!    26 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 100.0D+00

  y_stop(1:2) = (/ 0.667726D+00, -0.254738D+00 /)

  return
end
subroutine p36_title ( title )

!*****************************************************************************80
!
!! P36_TITLE returns the title of problem 36.
!
!  Modified:
!
!    26 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 36, Duffing''s Equation'

  return
end
function p37_autonomous ( )

!*****************************************************************************80
!
!! P37_AUTONOMOUS reports whether problem 37 is autonomous.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P37_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p37_autonomous

  p37_autonomous = .false.

  return
end
subroutine p37_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P37_EQUIL returns equilibrium solutions of problem 37.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1:2) = (/ 0.0D+00, 0.0D+00 /)
  else
    next = 0
  end if

  return
end
subroutine p37_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P37_FUN evaluates the function for problem 37.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) a
  real ( kind = 8 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) w
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  call p37_param ( 'GET', 'A', a )
  call p37_param ( 'GET', 'K', k )
  call p37_param ( 'GET', 'W', w )

  yp(1) = y(2)
  yp(2) = y(1) * ( 1.0D+00 - y(1)**2 ) - k * y(2) + a * cos ( w * t )

  return
end
subroutine p37_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P37_JAC evaluates the jacobian for problem 37.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) a
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) w
  real ( kind = 8 ) y(neqn)

  call p37_param ( 'GET', 'A', a )
  call p37_param ( 'GET', 'K', k )
  call p37_param ( 'GET', 'W', w )

  jac(1,1) =  0.0D+00
  jac(1,2) =  1.0D+00

  jac(2,1) = 1.0D+00 - 3.0D+00 * y(1)**2
  jac(2,2) = - k

  return
end
subroutine p37_neqn ( neqn )

!*****************************************************************************80
!
!! P37_NEQN returns the number of equations for problem 37.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 2

  return
end
subroutine p37_param ( action, name, value )

!*****************************************************************************80
!
!! P37_PARAM handles the parameters for problem 37.
!
!  Discussion:
!
!    A is the amplitude on the forcing term;
!    K is the damping coefficient.
!    W is the frequency coefficient on the forcing term;
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which should
!    be 'A', 'K' or 'W'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: a = 0.3D+00
  real ( kind = 8 ), save :: k = 0.2D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value
  real ( kind = 8 ), save :: w = 1.0D+00

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'A' ) ) then
      value = a
    else if ( s_eqi ( name, 'K' ) ) then
      value = k
    else if ( s_eqi ( name, 'W' ) ) then
      value = w
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P37_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'A' ) ) then
      a = value
    else if ( s_eqi ( name, 'K' ) ) then
      k = value
    else if ( s_eqi ( name, 'W' ) ) then
      w = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P37_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P37_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p37_scale ( neqn, scale )

!*****************************************************************************80
!
!! P37_SCALE returns scale factors for problem 37.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1:2) = (/ 1.0D+00, 1.0D+00 /)

  return
end
subroutine p37_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P37_START returns the starting point for problem 37.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1:neqn) = (/ 0.5D+00, 0.0D+00 /)

  return
end
subroutine p37_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P37_STOP returns the stopping point for problem 37.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  t_stop = 100.0D+00

  y_stop(1:2) = (/ -1.21774D+00, -0.548248D+00 /)

  return
end
subroutine p37_title ( title )

!*****************************************************************************80
!
!! P37_TITLE returns the title of problem 37.
!
!  Modified:
!
!    27 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 37, Duffing''s Equation with Damping and Forcing'

  return
end
function p38_autonomous ( )

!*****************************************************************************80
!
!! P38_AUTONOMOUS reports whether problem 38 is autonomous.
!
!  Modified:
!
!    14 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P38_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p38_autonomous

  p38_autonomous = .true.

  return
end
subroutine p38_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P38_EQUIL returns equilibrium solutions of problem 38.
!
!  Modified:
!
!    14 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    next = 1
    y(1) = 0.0D+00
  else if ( next == 1 ) then
    next = 2
    y(1) = 1.0D+00
  else
    next = 0
  end if

  return
end
subroutine p38_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P38_FUN evaluates the function for problem 38.
!
!  Discussion:
!
!    Moler attributes this problem to Lawrence Shampine.
!
!    The equation describes the radius of a ball of flame that
!    begins, at time 0, at DELTA.
!
!      Y(0) = DELTA
!
!    The rate of fuel consumption is proportional to the volume, and
!    the rate of fuel intake is proportional to the area of the ball.
!    We take the constant of proportionality to be 1.
!
!      Y' = Y^2 - Y3
!
!    The data is normalized so that Y = 1 is the equilibrium solution.
!
!    The computation is to be made from T = 0 to T = 2/DELTA.
!
!    For values of DELTA close to 1, such as 0.01, the equation is
!    not stiff.  But for DELTA = 0.0001, the equation can become
!    stiff as the solution approaches the equilibrium solution Y = 1,
!    and computed solutions may be wildly inaccurate or cautious
!    solvers may take very small timesteps.
!
!    The exact solution involves the Lambert W function, which
!    is defined by
!
!      W(z) * exp ( W(z) ) = z
!
!    and if we set
!
!      A = ( 1 / DELTA - 1 )
!
!    then
!
!      Y(T) = 1 / ( W(A*exp(A-T)) + 1 )
!
!  Modified:
!
!    14 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cleve Moler,
!    Cleve's Corner: Stiff Differential Equations,
!    MATLAB News and Notes,
!    May 2003, pages 12-13.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  yp(1) = y(1)**2 * ( 1.0D+00 - y(1) )

  return
end
subroutine p38_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P38_JAC evaluates the jacobian for problem 38.
!
!  Modified:
!
!    14 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) =  y(1) * ( 2.0D+00 - 3.0D+00 * y(1) )

  return
end
subroutine p38_neqn ( neqn )

!*****************************************************************************80
!
!! P38_NEQN returns the number of equations for problem 38.
!
!  Modified:
!
!    14 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 1

  return
end
subroutine p38_param ( action, name, value )

!*****************************************************************************80
!
!! P38_PARAM handles the parameters for problem 38.
!
!  Discussion:
!
!    DELTA is the initial radius of the flame.
!
!  Modified:
!
!    14 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which
!    should be 'DELTA'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: delta = 0.01D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      value = delta
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P38_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'DELTA' ) ) then
      delta = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P38_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P38_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p38_scale ( neqn, scale )

!*****************************************************************************80
!
!! P38_SCALE returns scale factors for problem 38.
!
!  Modified:
!
!    14 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1) = 1.0D+00

  return
end
subroutine p38_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P38_START returns the starting point for problem 38.
!
!  Modified:
!
!    14 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) delta
  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  call p38_param ( 'GET', 'DELTA', delta )

  t_start = 0.0D+00
  y_start(1) = delta

  return
end
subroutine p38_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P38_STOP returns the stopping point for problem 38.
!
!  Modified:
!
!    14 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) delta
  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  call p38_param ( 'GET', 'DELTA', delta )

  t_stop = 2.0D+00 / delta
  y_stop(1) = 1.0D+00

  return
end
subroutine p38_title ( title )

!*****************************************************************************80
!
!! P38_TITLE returns the title of problem 38.
!
!  Modified:
!
!    14 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 38, Shampine''s Ball of Flame'

  return
end
function p39_autonomous ( )

!*****************************************************************************80
!
!! P39_AUTONOMOUS reports whether problem 39 is autonomous.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P39_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p39_autonomous

  p39_autonomous = .false.

  return
end
subroutine p39_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P39_EQUIL returns equilibrium solutions of problem 39.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  next = 0

  return
end
subroutine p39_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P39_FUN evaluates the function for problem 39.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Arnold, John Polking,
!    Ordinary Differential Equations using Matlab,
!    Prentice Hall, 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  call p39_param ( 'GET', 'A', a )
  call p39_param ( 'GET', 'B', b )

  yp(1) = y(1)**2 - a * t + b

  return
end
subroutine p39_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P39_JAC evaluates the jacobian for problem 39.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  jac(1,1) = 2.0D+00 * y(1)

  return
end
subroutine p39_neqn ( neqn )

!*****************************************************************************80
!
!! P39_NEQN returns the number of equations for problem 39.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 1

  return
end
subroutine p39_param ( action, name, value )

!*****************************************************************************80
!
!! P39_PARAM handles the parameters for problem 39.
!
!  Discussion:
!
!    The ODE is
!
!      y' = y**2 - a * t + b
!
!    The parameters A and B may be controlled by this routine.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which
!    should be 'A' or 'B'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: a = 1.00D+00
  real ( kind = 8 ), save :: b = 0.00D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'A' ) ) then
      value = a
    else if ( s_eqi ( name, 'B' ) ) then
      value = b
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P39_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'A' ) ) then
      a = value
    else if ( s_eqi ( name, 'B' ) ) then
      b = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P39_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P39_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p39_scale ( neqn, scale )

!*****************************************************************************80
!
!! P39_SCALE returns scale factors for problem 39.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1) = 1.0D+00

  return
end
subroutine p39_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P39_START returns the starting point for problem 39.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = 0.0D+00
  y_start(1) = 0.5D+00

  return
end
subroutine p39_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P39_STOP returns the stopping point for problem 39.
!
!  Discussion:
!
!    I need to sit down and work these out...
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  call p39_param ( 'GET', 'A', a )
  call p39_param ( 'GET', 'B', b )

  if ( a /= 0.0D+00 ) then
    t_stop = ( b + 9.0D+00 ) / a
    y_stop(1) = -3.0D+00
  else if ( b /= 0.0D+00 ) then
    t_stop = 0.0D+00
    y_stop(1) = 0.0D+00
  else
    t_stop = 10.0D+00
    y_stop(1) = 0.1D+00
  end if

  return
end
subroutine p39_title ( title )

!*****************************************************************************80
!
!! P39_TITLE returns the title of problem 39.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 39, Polking''s first order ODE'

  return
end
function p40_autonomous ( )

!*****************************************************************************80
!
!! P40_AUTONOMOUS reports whether problem 40 is autonomous.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, logical P40_AUTONOMOUS, is TRUE if the equation is autonomous.
!
  implicit none

  logical p40_autonomous

  p40_autonomous = .false.

  return
end
subroutine p40_equil ( neqn, y, next )

!*****************************************************************************80
!
!! P40_EQUIL returns equilibrium solutions of problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) Y(NEQN), the "next" equilibrium solution, if any.
!
!    Input/output, integer ( kind = 4 ) NEXT, on input the index of the previous
!    equilibrium, which should be 0 on first call.  On output, the index
!    of the current equilibrium, or 0 if there are no more.
!
  implicit none

  integer ( kind = 4 ) neqn

  integer ( kind = 4 ) next
  real ( kind = 8 ) y(neqn)

  if ( next == 0 ) then
    y(1:neqn) = 0.0D+00
    next = 1
  else
    next = 0
  end if

  return
end
subroutine p40_fun ( neqn, t, y, yp )

!*****************************************************************************80
!
!! P40_FUN evaluates the function for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the derivative
!    function.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative function.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) eps
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  call p40_param ( 'GET', 'EPS', eps )

  yp(1) = y(1) * ( y(1) - t ) / eps

  return
end
subroutine p40_jac ( neqn, ldj, t, y, jac )

!*****************************************************************************80
!
!! P40_JAC evaluates the jacobian for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) LDJ, the leading dimension of JAC, which
!    must be at least NEQN.
!
!    Input, real ( kind = 8 ) T, Y(NEQN), the arguments of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(LDJ,NEQN), the jacobian matrix.
!
  implicit none

  integer ( kind = 4 ) ldj
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) eps
  real ( kind = 8 ) jac(ldj,neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)

  call p40_param ( 'GET', 'EPS', eps )

  jac(1,1) = ( 2.0D+00 * y(1) - t ) / eps

  return
end
subroutine p40_neqn ( neqn )

!*****************************************************************************80
!
!! P40_NEQN returns the number of equations for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NEQN, the number of equations.
!
  implicit none

  integer ( kind = 4 ) neqn

  neqn = 1

  return
end
subroutine p40_param ( action, name, value )

!*****************************************************************************80
!
!! P40_PARAM handles the parameters for problem 40.
!
!  Discussion:
!
!    The ODE is
!
!      y' = y ( y - t ) / eps
!
!    The parameter EPS may be controlled by this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is 'SET' to set a value, or 'GET' to
!    get a value.
!
!    Input, character ( len = * ) NAME, the name of the variable, which
!    should be EPS'.
!
!    Input/output, real ( kind = 8 ) VALUE, the value to assign to the
!    variable, or to retrieve from the variable.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: eps = 0.01D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'GET' ) ) then

    if ( s_eqi ( name, 'EPS' ) ) then
      value = eps
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P40_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else if ( s_eqi ( action, 'SET' ) ) then

    if ( s_eqi ( name, 'EPS' ) ) then
      eps = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P40_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized variable: "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P40_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action: "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p40_scale ( neqn, scale )

!*****************************************************************************80
!
!! P40_SCALE returns scale factors for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) SCALE(NEQN), the scaling factors.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) scale(neqn)

  scale(1) = 1.0D+00

  return
end
subroutine p40_start ( neqn, t_start, y_start )

!*****************************************************************************80
!
!! P40_START returns the starting point for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_START, Y_START(NEQN), the initial data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) t_start
  real ( kind = 8 ) y_start(neqn)

  t_start = -1.0D+00
  y_start(1) = -1.0D+00

  return
end
subroutine p40_stop ( neqn, t_stop, y_stop )

!*****************************************************************************80
!
!! P40_STOP returns the stopping point for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Output, real ( kind = 8 ) T_STOP, Y_STOP(NEQN), the final data.
!
  implicit none

  integer ( kind = 4 ) neqn

  real ( kind = 8 ) eps
  real ( kind = 8 ) error_f
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y_stop(neqn)

  call p40_param ( 'GET', 'EPS', eps )

  t_stop = 1.0D+00

  y_stop(1) = - 2.0D+00 * sqrt ( eps ) &
    * exp ( ( 1.0D+00 - t_stop * t_stop ) / ( 2.0D+00 * eps ) ) / &
    ( &
      2 * sqrt ( eps ) + &
      ( &
        error_f ( 1.0D+00 / sqrt ( 2.0D+00 * eps ) )  + &
        error_f ( t_stop  / sqrt ( 2.0D+00 + eps ) ) &
      ) &
      * exp ( 0.5D+00 / eps ) * sqrt ( 2.0D+00 * pi ) &
    )

  return
end
subroutine p40_title ( title )

!*****************************************************************************80
!
!! P40_TITLE returns the title of problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Problem 40, the Knee problem'

  return
end
function r8_cube_root ( x )

!*****************************************************************************80
!
!! R8_CUBE_ROOT returns the cube root of a real number.
!
!  Discussion:
!
!    This routine is designed to avoid the possible problems that can occur
!    when formulas like 0.0**(1/3) or (-1.0)**(1/3) are to be evaluated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose cube root is desired.
!
!    Output, real ( kind = 8 ) R8_CUBE_ROOT, the cube root of X.
!
  implicit none

  real ( kind = 8 ) r8_cube_root
  real ( kind = 8 ) x

  if ( 0.0D+00 < x ) then
    r8_cube_root = x**(1.0D+00/3.0D+00)
  else if ( x == 0.0D+00 ) then
    r8_cube_root = 0.0D+00
  else
    r8_cube_root = - ( abs ( x ) )**(1.0D+00/3.0D+00)
  end if

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
