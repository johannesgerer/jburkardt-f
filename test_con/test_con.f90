subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
!    which guarantee the ASCII collating sequence.
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
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character ch
  integer ( kind = 4 ) itemp

  itemp = iachar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
  end if

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC increments a partially numeric filename.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected.
!
!    If the name is empty, then the routine stops.
!
!    If the name contains no digits, the empty string is returned.
!
!  Example:
!
!      Input            Output
!      -----            ------
!      'a7to11.txt'     'a7to12.txt'
!      'a7to99.txt'     'a8to00.txt'
!      'a9to99.txt'     'a0to00.txt'
!      'cat.txt'        ' '
!      ' '              STOP!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  integer  ( kind = 4 )change
  integer ( kind = 4 ) digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  if ( lens <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_NAME_INC - Fatal error!'
    write ( *, '(a)' ) '  The input string is empty.'
    stop
  end if

  change = 0

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      change = change + 1

      digit = ichar ( c ) - 48
      digit = digit + 1

      if ( digit == 10 ) then
        digit = 0
      end if

      c = char ( digit + 48 )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  if ( change == 0 ) then
    file_name = ' '
    return
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine p00_fun ( problem, option, nvar, x, fx )

!*****************************************************************************80
!
!! P00_FUN evaluates the function for any problem.
!
!  Discussion:
!
!    These problems were collected by Professor Werner Rheinboldt, of the
!    University of Pittsburgh, and were used in the development of the
!    PITCON program.
!
!  Index:
!
!     1  The Freudenstein-Roth function
!     2  The Boggs function
!     3  The Powell function
!     4  The Broyden function
!     5  The Wacker function
!     6  The Aircraft stability function
!     7  The Cell kinetic function
!     8  The Riks mechanical problem
!     9  The Oden mechanical problem
!    10  Torsion of a square rod, finite difference solution
!    11  Torsion of a square rod, finite element solution
!    12  The materially nonlinear problem
!    13  Simpson's mildly nonlinear boundary value problem
!    14  Keller's boundary value problem
!    15  The Trigger Circuit
!    16  The Moore-Spence Chemical Reaction Integral Equation
!    17  The Bremermann Propane Combustion System
!    18  The semiconductor problem
!    19  The nitric acid absorption flash
!    20  The buckling spring
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Rami Melhem, Werner Rheinboldt,
!    A Comparison of Methods for Determining Turning Points of Nonlinear Equations,
!    Computing,
!    Volume 29, Number 3, September 1982, pages 201-226.
!
!    Werner Rheinboldt,
!    Numerical Analysis of Parameterized Nonlinear Equations,
!    Wiley, 1986,
!    ISBN: 0-471-88814-1,
!    LC: QA372.R54.
!
!    Werner Rheinboldt,
!    Sample Problems for Continuation Processes,
!    Technical Report ICMA-80-?,
!    Institute for Computational Mathematics and Applications,
!    Department of Mathematics,
!    University of Pittsburgh, November 1980.
!
!    Werner Rheinboldt, John Burkardt,
!    A Locally Parameterized Continuation Process,
!    ACM Transactions on Mathematical Software,
!    Volume 9, Number 2, June 1983, pages 215-235.
!
!    Werner Rheinboldt, John Burkardt,
!    Algorithm 596:
!    A Program for a Locally Parameterized
!    Continuation Process,
!    ACM Transactions on Mathematical Software,
!    Volume 9, Number 2, June 1983, pages 236-241.
!
!    Werner Rheinboldt,
!    Computation of Critical Boundaries on Equilibrium Manifolds,
!    SIAM Journal on Numerical Analysis,
!    Volume 19, Number 3, June 1982, pages 653-669.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar-1)
  integer ( kind = 4 ) option
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(nvar)

  if ( problem == 1 ) then
    call p01_fun ( option, nvar, x, fx )
  else if ( problem == 2 ) then
    call p02_fun ( option, nvar, x, fx )
  else if ( problem == 3 ) then
    call p03_fun ( option, nvar, x, fx )
  else if ( problem == 4 ) then
    call p04_fun ( option, nvar, x, fx )
  else if ( problem == 5 ) then
    call p05_fun ( option, nvar, x, fx )
  else if ( problem == 6 ) then
    call p06_fun ( option, nvar, x, fx )
  else if ( problem == 7 ) then
    call p07_fun ( option, nvar, x, fx )
  else if ( problem == 8 ) then
    call p08_fun ( option, nvar, x, fx )
  else if ( problem == 9 ) then
    call p09_fun ( option, nvar, x, fx )
  else if ( problem == 10 ) then
    call p10_fun ( option, nvar, x, fx )
  else if ( problem == 11 ) then
    call p11_fun ( option, nvar, x, fx )
  else if ( problem == 12 ) then
    call p12_fun ( option, nvar, x, fx )
  else if ( problem == 13 ) then
    call p13_fun ( option, nvar, x, fx )
  else if ( problem == 14 ) then
    call p14_fun ( option, nvar, x, fx )
  else if ( problem == 15 ) then
    call p15_fun ( option, nvar, x, fx )
  else if ( problem == 16 ) then
    call p16_fun ( option, nvar, x, fx )
  else if ( problem == 17 ) then
    call p17_fun ( option, nvar, x, fx )
  else if ( problem == 18 ) then
    call p18_fun ( option, nvar, x, fx )
  else if ( problem == 19 ) then
    call p19_fun ( option, nvar, x, fx )
  else if ( problem == 20 ) then
    call p20_fun ( option, nvar, x, fx )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FUN - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized problem number = ', problem
    stop
  end if

  return
end
subroutine p00_jac ( problem, option, nvar, x, jac )

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
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) jac(nvar,nvar)
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(nvar)

  if ( problem == 1 ) then
    call p01_jac ( option, nvar, x, jac )
  else if ( problem == 2 ) then
    call p02_jac ( option, nvar, x, jac )
  else if ( problem == 3 ) then
    call p03_jac ( option, nvar, x, jac )
  else if ( problem == 4 ) then
    call p04_jac ( option, nvar, x, jac )
  else if ( problem == 5 ) then
    call p05_jac ( option, nvar, x, jac )
  else if ( problem == 6 ) then
    call p06_jac ( option, nvar, x, jac )
  else if ( problem == 7 ) then
    call p07_jac ( option, nvar, x, jac )
  else if ( problem == 8 ) then
    call p08_jac ( option, nvar, x, jac )
  else if ( problem == 9 ) then
    call p09_jac ( option, nvar, x, jac )
  else if ( problem == 10 ) then
    call p10_jac ( option, nvar, x, jac )
  else if ( problem == 11 ) then
    call p11_jac ( option, nvar, x, jac )
  else if ( problem == 12 ) then
    call p12_jac ( option, nvar, x, jac )
  else if ( problem == 13 ) then
    call p13_jac ( option, nvar, x, jac )
  else if ( problem == 14 ) then
    call p14_jac ( option, nvar, x, jac )
  else if ( problem == 15 ) then
    call p15_jac ( option, nvar, x, jac )
  else if ( problem == 16 ) then
    call p16_jac ( option, nvar, x, jac )
  else if ( problem == 17 ) then
    call p17_jac ( option, nvar, x, jac )
  else if ( problem == 18 ) then
    call p18_jac ( option, nvar, x, jac )
  else if ( problem == 19 ) then
    call p19_jac ( option, nvar, x, jac )
  else if ( problem == 20 ) then
    call p20_jac ( option, nvar, x, jac )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_JAC - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized problem number = ', problem
    stop
  end if
!
!  Guarantee that the last row is zeroed out.
!
  jac(nvar,1:nvar) = 0.0D+00

  return
end
subroutine p00_jac_check ( problem, option, nvar, x, max_adif, max_adif_i, &
  max_adif_j, max_rdif, max_rdif_i, max_rdif_j )

!*****************************************************************************80
!
!! P00_JAC_CHECK compares the jacobian with a finite difference estimate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) MAX_ADIF, the maximum absolute difference.
!
!    Output, integer ( kind = 4 ) MAX_ADIF_I, MAX_ADIF_J, the indices where
!    the maximmum absolute difference was found.
!
!    Output, real ( kind = 8 ) MAX_RDIF, the maximum relative difference.
!
!    Output, integer ( kind = 4 ) MAX_RDIF_I, MAX_RDIF_J, the indices where
!    the maximmum relative difference was found.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) dif
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) jac_dif(nvar,nvar)
  real ( kind = 8 ) max_adif
  integer ( kind = 4 ) max_adif_i
  integer ( kind = 4 ) max_adif_j
  real ( kind = 8 ) max_rdif
  integer ( kind = 4 ) max_rdif_i
  integer ( kind = 4 ) max_rdif_j
  integer ( kind = 4 ) problem
  real ( kind = 8 ), parameter :: rel = 0.0001D+00
  real ( kind = 8 ) x(nvar)
!
!  Compute the jacobian.
!
  call p00_jac ( problem, option, nvar, x, jac )
!
!  Estimate the jacobian via finite differences.
!
  call p00_jac_dif ( problem, option, nvar, x, jac_dif )
!
!  Compare the jacobians.
!
  max_rdif = 0.0D+00
  max_rdif_i = 0
  max_rdif_j = 0

  max_adif = 0.0D+00
  max_adif_i = 0
  max_adif_j = 0

  do i = 1, nvar - 1
    do j = 1, nvar

      dif = abs ( jac(i,j) - jac_dif(i,j) )

      if ( max_adif < dif ) then
        max_adif = dif
        max_adif_i = i
        max_adif_j = j
      end if

      if ( rel < abs ( jac(i,j) ) ) then
        if ( max_rdif * abs ( jac(i,j) ) < dif ) then
          max_rdif = dif / abs ( jac(i,j) )
          max_rdif_i = i
          max_rdif_j = j
        end if
      end if

    end do
  end do

  return
end
subroutine p00_jac_dif ( problem, option, nvar, x, jac_dif )

!*****************************************************************************80
!
!! P00_JAC_DIF estimates the jacobian via finite differences.
!
!  Discussion:
!
!    This is a relatively unsophisticated way of estimating the jacobian.
!    The value of the internal parameter REL, set below, can affect
!    the results in a strong way.  If the jacobian reported by this
!    routine seems unsatisfactory, check the results for values of
!    REL that are 10 times larger and smaller, and see if the trend
!    makes sense.  Values of REL that are too large for a given
!    problem will make crude estimates, but values that are too small
!    will result in roundoff, and in severe cases, the computation of
!    zeroes in the jacobian.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC_DIF(NVAR,NVAR), an estimate of the jacobian.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) dx
  real ( kind = 8 ) fxm(nvar)
  real ( kind = 8 ) fxp(nvar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac_dif(nvar,nvar)
  integer ( kind = 4 ) problem
  real ( kind = 8 ), parameter :: REL = 0.0001D+00
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xsave
!
!  Perturb each variable.
!
  do j = 1, nvar
!
!  Save X(J), and compute a perturbation.
!
    xsave = x(j)
    dx = REL * ( abs ( xsave ) + 1.0D+00 )
!
!  Compute the function value at X + dX.
!
    x(j) = xsave + dx

    call p00_fun ( problem, option, nvar, x, fxp )
!
!  Compute the function value at X - dX.
!
    x(j) = xsave - dx

    call p00_fun ( problem, option, nvar, x, fxm )
!
!  Restore X(J).
!
    x(j) = xsave
!
!  Compute column J of the finite difference jacobian.
!
    do i = 1, nvar - 1
      jac_dif(i,j) = 0.5D+00 * ( fxp(i) - fxm(i) ) / dx
    end do

  end do

  return
end
subroutine p00_limit ( problem, option, nvar, x1, tan1, x2, tan2, lim, &
  x, tan, status )

!*****************************************************************************80
!
!! P00_LIMIT seeks a limit point.
!
!  Discussion:
!
!    For a given index 1 <= LIM <= NVAR, a limit point X is a point which
!    satisfies F(X) = 0 and TAN(X)(LIM) = 0, that is, X is a point on the
!    solution curve, and the LIM-th component of the tangent vector at X
!    is zero.
!
!    This function may be called if a limit point has been bracketed,
!    that is, if X1 and X2 are points on the curve with the property that
!    there is a change in sign in the LIM-th component of the tangent
!    vector between X1 and X2.
!
!    The function carries out an iteration seeking a point X between
!    X1 and X2 for which the LIM-th tangent component is zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X1(NVAR), TAN1(NVAR), a point on the curve,
!    and its tangent vector.
!
!    Input, real ( kind = 8 ) X2(NVAR), TAN2(NVAR), a second point on the curve,
!    and its tangent vector.
!
!    Input, integer ( kind = 4 ) LIM, the index of the entry of TAN which
!    we are seeking to zero.
!
!    Output, real ( kind = 8 ) X(NVAR), TAN(NVAR), the computed limit point
!    and its tangent vector.
!
!    Output, integer ( kind = 4 ) STATUS.
!    nonnegative, the limit point was computed in STATUS steps.
!    negative, the limit point could not be computed.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  integer ( kind = 4 ) lim
  integer ( kind = 4 ) option
  integer ( kind = 4 ) par_index
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) status
  integer ( kind = 4 ) status_zero
  integer ( kind = 4 ) status_newton
  real ( kind = 8 ) tan(nvar)
  real ( kind = 8 ) tan1(nvar)
  real ( kind = 8 ) tan2(nvar)
  real ( kind = 8 ) tol
  real ( kind = 8 ) value
  logical, parameter :: VERBOSE = .false.
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) x1(nvar)
  real ( kind = 8 ) x2(nvar)
!
!  Use a fixed parameter index, but do NOT use LIM.
!
  x(1:nvar) = x2(1:nvar) - x1(1:nvar)
  x(lim) = 0.0D+00

  call r8vec_amax_index ( nvar, x, par_index )
!
!  Start the zero finding process.
!
  a = 0.0D+00
  b = 1.0D+00
  tol = sqrt ( epsilon ( tol ) )
  arg = 0.0D+00
  status_zero = 0
  value = 0.0D+00

  status = 0

  do

    call zero_rc ( a, b, tol, arg, status_zero, value )

    if ( status_zero < 0 ) then
      status = -1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P00_LIMIT - Fatal error!'
      write ( *, '(a)' ) '  ZERO_RC returned an error flag!'
      exit
    end if

    if ( arg == 0.0D+00 ) then

      x(1:nvar) = x1(1:nvar)
      tan(1:nvar) = tan1(1:nvar)

    else if ( arg == 1.0D+00 ) then

      x(1:nvar) = x2(1:nvar)
      tan(1:nvar) = tan2(1:nvar)

    else

      x(1:nvar) = ( 1.0D+00 - arg ) * x1(1:nvar)   &
                  +           arg   * x2(1:nvar)

      call p00_newton ( problem, option, nvar, x, par_index, status_newton )

      if ( status_newton < 0 ) then
        status = -2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_LIMIT - Fatal error!'
        write ( *, '(a)' ) '  ZERO_RC returned an error flag!'
        exit
      end if

      call p00_tan ( problem, option, nvar, x, tan )

    end if

    value = tan(lim)

    if ( VERBOSE ) then
      write ( *, '(2x,i8,2x,g14.8,2x,g14.8)' ) status_zero, arg, value
    end if

    status = status + 1

    if ( status_zero == 0 ) then
      exit
    end if

  end do

  return
end
subroutine p00_newton ( problem, option, nvar, x, par_index, status )

!*****************************************************************************80
!
!! P00_NEWTON applies Newton's method to an approximate root.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input/output, real ( kind = 8 ) X(NVAR).
!    On input, the starting point of Newton's method.
!    On output, an improved estimate of the root of F(X)=0, if the
!    algorithm converged.
!
!    Input, integer ( kind = 4 ) PAR_INDEX, the index of the parameter
!    to be fixed.  This variable should be between 1 and NVAR.  However, the user can
!    set it to 0, indicating that the program should make an intelligent
!    choice for the index.
!
!    Output, integer ( kind = 4 ) STATUS, the status of the iteration.
!    -3, the full number of steps was taken without convergence.
!        (however, the output X might be CLOSE to a good solution).
!    -2, the iteration seemed to be diverging, and was halted.
!    -1, the jacobian was singular, and the iteration was halted.
!     nonnegative, the convergence test was satisfied, and this is the
!        number of steps taken (possibly 0).
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar)
  real ( kind = 8 ), parameter :: FX_ABS_TOL = 0.000001D+00
  real ( kind = 8 ) fx_max
  real ( kind = 8 ) fx_max_init
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ipivot(nvar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) option
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: IT_MAX = 20
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) par_index
  real ( kind = 8 ) par_value
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) status
  logical, parameter :: VERBOSE = .false.
  real ( kind = 8 ) x(nvar)

  if ( par_index < 1 .or. nvar < par_index ) then
    call p00_par_index ( problem, option, nvar, x, par_index )

    if ( VERBOSE ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) '  Iteration will hold index ', par_index, ' fixed.'
    end if

  end if
  par_value = x(par_index)

  if ( VERBOSE ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_NEWTON'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Step   F(X)'
    write ( *, '(a)' ) ' '
  end if

  do it = 0, IT_MAX
!
!  Compute the function value.
!
    call p00_fun ( problem, option, nvar, x, fx )
    fx(nvar) = x(par_index) - par_value
!
!  Compute the norm of the function value.
!
    fx_max = maxval ( abs ( fx(1:nvar) ) )

    if ( VERBOSE ) then
      write ( *, '(2x,i8,2x,g14.6)' ) it, fx_max
    end if

    if ( it == 0 ) then
      fx_max_init = fx_max
    end if
!
!  If the function norm is small enough, return.
!
    if ( abs ( fx_max ) < FX_ABS_TOL ) then
      status = it
      exit
    end if
!
!  If the function norm seems to be exploding, halt.
!
    if ( 1000.0 * fx_max_init < abs ( fx_max ) ) then
      status = -2
      exit
    end if

    if ( it == IT_MAX ) then
      status = -3
      exit
    end if
!
!  Compute the jacobian.
!
    call p00_jac ( problem, option, nvar, x, jac )

    jac(nvar,1:nvar) = 0.0D+00
    jac(nvar,par_index) = 1.0D+00
!
!  Factor the jacobian.
!
    call sge_fa ( nvar, nvar, jac, ipivot, info )

    if ( info /= 0 ) then
      status = -1
      exit
    end if
!
!  Solve the system JAC * DX = FX
!
    job = 0
    call sge_sl ( nvar, nvar, jac, ipivot, fx, job )
!
!  Update X = X - DX.
!
    x(1:nvar) = x(1:nvar) - fx(1:nvar)

  end do

  return
end
subroutine p00_nvar ( problem, option, nvar )

!*****************************************************************************80
!
!! P00_NVAR sets the number of variables for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) nvar

  if (  problem == 1 ) then
    call p01_nvar ( option, nvar )
  else if ( problem == 2 ) then
    call p02_nvar ( option, nvar )
  else if ( problem == 3 ) then
    call p03_nvar ( option, nvar )
  else if ( problem == 4 ) then
    call p04_nvar ( option, nvar )
  else if ( problem == 5 ) then
    call p05_nvar ( option, nvar )
  else if ( problem == 6 ) then
    call p06_nvar ( option, nvar )
  else if ( problem == 7 ) then
    call p07_nvar ( option, nvar )
  else if ( problem == 8 ) then
    call p08_nvar ( option, nvar )
  else if ( problem == 9 ) then
    call p09_nvar ( option, nvar )
  else if ( problem == 10 ) then
    call p10_nvar ( option, nvar )
  else if ( problem == 11 ) then
    call p11_nvar ( option, nvar )
  else if ( problem == 12 ) then
    call p12_nvar ( option, nvar )
  else if ( problem == 13 ) then
    call p13_nvar ( option, nvar )
  else if ( problem == 14 ) then
    call p14_nvar ( option, nvar )
  else if ( problem == 15 ) then
    call p15_nvar ( option, nvar )
  else if ( problem == 16 ) then
    call p16_nvar ( option, nvar )
  else if ( problem == 17 ) then
    call p17_nvar ( option, nvar )
  else if ( problem == 18 ) then
    call p18_nvar ( option, nvar )
  else if ( problem == 19 ) then
    call p19_nvar ( option, nvar )
  else if ( problem == 20 ) then
    call p20_nvar ( option, nvar )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_NVAR - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized problem index  = ', problem
    stop
  end if

  return
end
subroutine p00_option_num ( problem, option_num )

!*****************************************************************************80
!
!! P00_OPTION_NUM returns the number of options available for a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options available for
!    this problem.  OPTION_NUM is always at least 1.
!
  implicit none

  integer ( kind = 4 ) option_num
  integer ( kind = 4 ) problem

  if (  problem == 1 ) then
    call p01_option_num ( option_num )
  else if ( problem == 2 ) then
    call p02_option_num ( option_num )
  else if ( problem == 3 ) then
    call p03_option_num ( option_num )
  else if ( problem == 4 ) then
    call p04_option_num ( option_num )
  else if ( problem == 5 ) then
    call p05_option_num ( option_num )
  else if ( problem == 6 ) then
    call p06_option_num ( option_num )
  else if ( problem == 7 ) then
    call p07_option_num ( option_num )
  else if ( problem == 8 ) then
    call p08_option_num ( option_num )
  else if ( problem == 9 ) then
    call p09_option_num ( option_num )
  else if ( problem == 10 ) then
    call p10_option_num ( option_num )
  else if ( problem == 11 ) then
    call p11_option_num ( option_num )
  else if ( problem == 12 ) then
    call p12_option_num ( option_num )
  else if ( problem == 13 ) then
    call p13_option_num ( option_num )
  else if ( problem == 14 ) then
    call p14_option_num ( option_num )
  else if ( problem == 15 ) then
    call p15_option_num ( option_num )
  else if ( problem == 16 ) then
    call p16_option_num ( option_num )
  else if ( problem == 17 ) then
    call p17_option_num ( option_num )
  else if ( problem == 18 ) then
    call p18_option_num ( option_num )
  else if ( problem == 19 ) then
    call p19_option_num ( option_num )
  else if ( problem == 20 ) then
    call p20_option_num ( option_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_OPTION_NUM - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized problem index  = ', problem
    stop
  end if

  return
end
subroutine p00_par_index ( problem, option, nvar, x, par_index )

!*****************************************************************************80
!
!! P00_PAR_INDEX chooses the index of the continuation parameter.
!
!  Discussion:
!
!    Given the NVAR-dimensional point X, the (NVAR-1)-dimensional function
!    F(X), and the NVAR-1 by NVAR jacobian matrix, let the NVAR-dimensional
!    vector TAN be any null vector of JAC.
!
!      JAC * TAN = 0
!
!    Choose PAR_INDEX to be the index of TAN of maximum absolute value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real X(NVAR), the starting point of Newton's method.
!
!    Output, integer ( kind = 4 ) PAR_INDEX, the index of the parameter
!    to be held fixed.  This variable will be between 1 and NVAR.  It is
!    the index of the variable which is currently changing most rapidly
!    along the curve F(X) = 0.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  integer ( kind = 4 ) par_index
  integer ( kind = 4 ) problem
  real ( kind = 8 ) tan(nvar)
  real ( kind = 8 ) x(nvar)

  call p00_tan ( problem, option, nvar, x, tan )

  call r8vec_amax_index ( nvar, tan, par_index )

  return
end
subroutine p00_problem_num ( problem_num )

!*****************************************************************************80
!
!! P00_PROBLEM_NUM returns the number of problems available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) PROBLEM_NUM, the number of problems.
!
  implicit none

  integer ( kind = 4 ) problem_num

  problem_num = 20

  return
end
subroutine p00_start ( problem, option, nvar, x )

!*****************************************************************************80
!
!! P00_START returns a starting point for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  integer ( kind = 4 ) problem
  real ( kind = 8 ) x(nvar)

  if ( problem == 1 ) then
    call p01_start ( option, nvar, x )
  else if ( problem == 2 ) then
    call p02_start ( option, nvar, x )
  else if ( problem == 3 ) then
    call p03_start ( option, nvar, x )
  else if ( problem == 4 ) then
    call p04_start ( option, nvar, x )
  else if ( problem == 5 ) then
    call p05_start ( option, nvar, x )
  else if ( problem == 6 ) then
    call p06_start ( option, nvar, x )
  else if ( problem == 7 ) then
    call p07_start ( option, nvar, x )
  else if ( problem == 8 ) then
    call p08_start ( option, nvar, x )
  else if ( problem == 9 ) then
    call p09_start ( option, nvar, x )
  else if ( problem == 10 ) then
    call p10_start ( option, nvar, x )
  else if ( problem == 11 ) then
    call p11_start ( option, nvar, x )
  else if ( problem == 12 ) then
    call p12_start ( option, nvar, x )
  else if ( problem == 13 ) then
    call p13_start ( option, nvar, x )
  else if ( problem == 14 ) then
    call p14_start ( option, nvar, x )
  else if ( problem == 15 ) then
    call p15_start ( option, nvar, x )
  else if ( problem == 16 ) then
    call p16_start ( option, nvar, x )
  else if ( problem == 17 ) then
    call p17_start ( option, nvar, x )
  else if ( problem == 18 ) then
    call p18_start ( option, nvar, x )
  else if ( problem == 19 ) then
    call p19_start ( option, nvar, x )
  else if ( problem == 20 ) then
    call p20_start ( option, nvar, x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_START - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized problem index = ', problem
    stop
  end if

  return
end
subroutine p00_step ( problem, option, nvar, x, par_index, h, hmin, hmax, &
  status )

!*****************************************************************************80
!
!! P00_STEP takes one continuation step.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real X(NVAR), the starting point.
!
!    Input, integer ( kind = 4 ) PAR_INDEX, the continuation parameter.
!    If the program is free to choose this value, set it to 0.
!
!    Input, real H, HMIN, HMAX, the suggested step, and the minimum
!    and maximum stepsizes.  H may be negative.
!
!    Output, real X(NVAR), the new point, if STATUS = 0.
!
!    Output, real H, the stepsize that was used.
!
!    Output, integer ( kind = 4 ) STATUS, the status of the calculation.
!    0, successful.
!    nonzero, the Newton iteration failed repeatedly even when the
!    minimum stepsize was used.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option
  integer ( kind = 4 ) par_index
  integer ( kind = 4 ) problem
  real ( kind = 8 ) r8_sign
  integer ( kind = 4 ) status
  real ( kind = 8 ) tan(nvar)
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xt(nvar)
!
!  Compute the tangent.
!
  call p00_tan ( problem, option, nvar, x, tan )
!
!  Estimate the next point.
!
  do

    xt(1:nvar) = x(1:nvar) + h * tan(1:nvar)
!
!  Use the Newton method.
!
    call p00_newton ( problem, option, nvar, xt, par_index, status )

    if ( status == 0 ) then
      exit
    end if

    if ( abs ( h ) <= hmin ) then
      exit
    end if

    if ( hmin < abs ( h ) ) then
      h = r8_sign ( h ) * max ( abs ( h ) / 2.0D+00, hmin )
    end if

  end do

  x(1:nvar) = xt(1:nvar)

  return
end
subroutine p00_stepsize ( problem, option, h, hmin, hmax )

!*****************************************************************************80
!
!! P00_STEPSIZE returns step sizes for any problem.
!
!  Discussion:
!
!    The routine returns a suggested initial stepsize, and suggestions for
!    the minimum and maximum stepsizes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option
  integer ( kind = 4 ) problem

  if ( problem == 1 ) then
    call p01_stepsize ( option, h, hmin, hmax )
  else if ( problem == 2 ) then
    call p02_stepsize ( option, h, hmin, hmax )
  else if ( problem == 3 ) then
    call p03_stepsize ( option, h, hmin, hmax )
  else if ( problem == 4 ) then
    call p04_stepsize ( option, h, hmin, hmax )
  else if ( problem == 5 ) then
    call p05_stepsize ( option, h, hmin, hmax )
  else if ( problem == 6 ) then
    call p06_stepsize ( option, h, hmin, hmax )
  else if ( problem == 7 ) then
    call p07_stepsize ( option, h, hmin, hmax )
  else if ( problem == 8 ) then
    call p08_stepsize ( option, h, hmin, hmax )
  else if ( problem == 9 ) then
    call p09_stepsize ( option, h, hmin, hmax )
  else if ( problem == 10 ) then
    call p10_stepsize ( option, h, hmin, hmax )
  else if ( problem == 11 ) then
    call p11_stepsize ( option, h, hmin, hmax )
  else if ( problem == 12 ) then
    call p12_stepsize ( option, h, hmin, hmax )
  else if ( problem == 13 ) then
    call p13_stepsize ( option, h, hmin, hmax )
  else if ( problem == 14 ) then
    call p14_stepsize ( option, h, hmin, hmax )
  else if ( problem == 15 ) then
    call p15_stepsize ( option, h, hmin, hmax )
  else if ( problem == 16 ) then
    call p16_stepsize ( option, h, hmin, hmax )
  else if ( problem == 17 ) then
    call p17_stepsize ( option, h, hmin, hmax )
  else if ( problem == 18 ) then
    call p18_stepsize ( option, h, hmin, hmax )
  else if ( problem == 19 ) then
    call p19_stepsize ( option, h, hmin, hmax )
  else if ( problem == 20 ) then
    call p20_stepsize ( option, h, hmin, hmax )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_STEPSIZE - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized problem number = ', problem
    stop
  end if

  return
end
subroutine p00_tan ( problem, option, nvar, x, tan )

!*****************************************************************************80
!
!! P00_TAN determines a tangent vector at X.
!
!  Discussion:
!
!    If X is a solution of F(Y) = 0, then the vector TAN
!    is tangent to the curve of solutions at X.
!
!    If X is not a solution of F(Y) = 0, then the vector TAN
!    is tangent to the curve F(Y) = F(X) at X.
!
!    The vector will have unit euclidean norm.
!
!    The sign of TAN will be chosen so that the determinant
!    of F'(X) augmented with a final row equal to TAN will be positive.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the evaluation point.
!
!    Output, real ( kind = 8 ) TAN(NVAR), a tangent vector at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) jac_det
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: nullspace
  integer ( kind = 4 ) nullspace_size
  integer ( kind = 4 ) option
  integer ( kind = 4 ) problem
  real ( kind = 8 ) tan(nvar)
  real ( kind = 8 ) tan_norm
  real ( kind = 8 ) x(nvar)
!
!  Compute the jacobian.
!
  call p00_jac ( problem, option, nvar, x, jac )
!
!  Compute the nullspace size.
!
  call r8mat_nullspace_size ( nvar, nvar, jac, nullspace_size )

  if ( nullspace_size < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TAN - Fatal error!'
    write ( *, '(a)' ) '  The matrix seems to have no nullspace.'
    write ( *, '(a)' ) '  The tangent vector could not be computed.'
    stop
  end if
!
!  Compute the nullspace.
!
  allocate ( nullspace(1:nvar,1:nullspace_size) )

  call r8mat_nullspace ( nvar, nvar, jac, nullspace_size, nullspace )

  tan(1:nvar) = nullspace(1:nvar,1)

  deallocate ( nullspace )
!
!  Choose the sign of TAN by the determinant condition.
!
  jac(nvar,1:nvar) = tan(1:nvar)

  call r8mat_det ( nvar, jac, jac_det )

  if ( jac_det < 0.0D+00 ) then
    tan(1:nvar) = - tan(1:nvar)
  end if

  tan_norm = sqrt ( sum ( tan(1:nvar)**2 ) )

  tan(1:nvar) = tan(1:nvar) / tan_norm

  return
end
subroutine p00_target ( problem, option, nvar, x1, x2, tar_index, tar_value, &
  x, status )

!*****************************************************************************80
!
!! P00_TARGET computes a solution with a given component value.
!
!  Discussion:
!
!    If we write G(X) = X(TAR_INDEX) - TAR_VALUE, then we are seeking a
!    solution of
!
!      ( F(X) )  =  ( 0 )
!      ( G(X) )     ( 0 )
!
!    We treat the index TAR_INDEX as the parameter to be held fixed.
!
!    Typically, this routine would be called when the user has computed
!    two successive solutions X1 and X2, with the property that the
!
!      X1(TAR_INDEX) < TAR_VALUE < X2(TAR_INDEX)
!
!    or vice-versa.
!
!    In that case, the appropriate estimate for the starting point X is
!
!      X = ( ( X2(TAR_INDEX) - TAR_VALUE                 ) * X1
!          + (                 TAR_VALUE - X1(TAR_INDEX) ) * X2 )
!          / ( X2(TAR_INDEX)             - X1(TAR_INDEX) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X1(NVAR), X2(NVAR), two points satisying F(X) = 0.
!    It is assumed that X1(TAR_INDEX) and X2(TAR_INDEX) bracket the
!    desired value TAR_VALUE.
!
!    Input, integer ( kind = 4 ) TAR_INDEX, the index of the entry of X whose
!    value is being specified.
!
!    Input, real ( kind = 8 ) TAR_VALUE, the desired value of X(TAR_INDEX).
!
!    Output, integer ( kind = 4 ) STATUS.
!    nonnegative, the target point was computed.
!    negative, the target point could not be computed.
!
!    Output, real ( kind = 8 ) X(NVAR), a point satisfying F(X) = 0 and
!    X(TAR_INDEX)=TAR_VALUE.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  integer ( kind = 4 ) par_index
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) status
  integer ( kind = 4 ) tar_index
  real ( kind = 8 ) tar_value
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) x1(nvar)
  real ( kind = 8 ) x2(nvar)

  x(1:nvar) = ( ( x2(tar_index) - tar_value                 ) * x1(1:nvar)   &
              + (                 tar_value - x1(tar_index) ) * x2(1:nvar) ) &
              / ( x2(tar_index)             - x1(tar_index) );

  par_index = tar_index

  call p00_newton ( problem, option, nvar, x, par_index, status )

  return
end
subroutine p00_title ( problem, option, title )

!*****************************************************************************80
!
!! P00_TITLE sets the title for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) problem
  character ( len = * ) title

  if (  problem == 1 ) then
    call p01_title ( option, title )
  else if ( problem == 2 ) then
    call p02_title ( option, title )
  else if ( problem == 3 ) then
    call p03_title ( option, title )
  else if ( problem == 4 ) then
    call p04_title ( option, title )
  else if ( problem == 5 ) then
    call p05_title ( option, title )
  else if ( problem == 6 ) then
    call p06_title ( option, title )
  else if ( problem == 7 ) then
    call p07_title ( option, title )
  else if ( problem == 8 ) then
    call p08_title ( option, title )
  else if ( problem == 9 ) then
    call p09_title ( option, title )
  else if ( problem == 10 ) then
    call p10_title ( option, title )
  else if ( problem == 11 ) then
    call p11_title ( option, title )
  else if ( problem == 12 ) then
    call p12_title ( option, title )
  else if ( problem == 13 ) then
    call p13_title ( option, title )
  else if ( problem == 14 ) then
    call p14_title ( option, title )
  else if ( problem == 15 ) then
    call p15_title ( option, title )
  else if ( problem == 16 ) then
    call p16_title ( option, title )
  else if ( problem == 17 ) then
    call p17_title ( option, title )
  else if ( problem == 18 ) then
    call p18_title ( option, title )
  else if ( problem == 19 ) then
    call p19_title ( option, title )
  else if ( problem == 20 ) then
    call p20_title ( option, title )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized problem index = ', problem
    stop
  end if

  return
end
subroutine p01_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P01_FUN evaluates the function for problem 1.
!
!  Title:
!
!    The Freudenstein-Roth function
!
!  Description:
!
!    One way to use a continuation code as a nonlinear root finder
!    is to start with a set of nonlinear equations G(X), and an
!    approximate root A, and create a "homotopy" function F(X,Y)
!    with the properties that F(A,0.0) = 0 and F(X,1.0) = G(X).
!    Thus, the homotopy function F has a known exact solution
!    from which we can start with no difficulty.  If the continuation
!    code can take us from Y = 0 to Y = 1, then we have found
!    an X so that F(X,1.0) = 0, so we have found a solution to G(X)=0.
!
!    The Freudenstein-Roth function F(X) is derived in this way
!    from a homotopy of G(X):
!
!      F ( X(1), X(2), X(3) ) =
!        G ( X(1), X(2) ) - ( 1 - X(3) ) * G ( Y1, Y2 )
!
!    where Y1 and Y2 are some fixed values, and
!
!      G(1) = X(1) - X(2)*X(2)*X(2) + 5*X(2)*X(2) -  2*X(2) - 13
!      G(2) = X(1) + X(2)*X(2)*X(2) +   X(2)*X(2) - 14*X(2) - 29
!
!  Options 1, 2, 3:
!
!    The starting point is X0 = ( 15, -2, 0 ).
!
!    A great deal of information is available about the homotopy curve
!    generated by this starting point:
!
!    The function F(X) has the form
!
!      F(1) = X(1) - X(2)**3 + 5*X(2)**2 -  2*X(2) - 13 + 34*(X(3)-1)
!      F(2) = X(1) + X(2)**3 +   X(2)**2 - 14*X(2) - 29 + 10*(X(3)-1)
!
!    There is a closed form representation of the curve in terms of the
!    second parameter:
!
!      X(1) = (-11*X(2)**3 + 4*X(2)**2 + 114*X(2) + 214) /  6
!      X(2) = X(2)
!      X(3) = (    X(2)**3 - 2*X(2)**2 -   6*X(2) +   4) / 12
!
!    The first option simply requests the production of solution points
!    along the curve until a point is reached whose third component is
!    exactly 1.
!
!    Options 2 and 3 use the same starting point, and also stop when the
!    third component is 1.  However, these options in addition search
!    for limit points in the first and third components of the solution,
!    respectively.
!
!    The target solution has X(3) = 1, and is ( 5, 4, 1 ).
!
!    Limit points for X1:
!
!      ( 14.28309, -1.741377,  0.2585779 )
!      ( 61.66936,  1.983801, -0.6638797 )
!
!    Limit points for X3:
!
!     (20.48586, -0.8968053, 0.5875873)
!     (61.02031,  2.230139, -0.6863528)
!
!    The curve has several dramatic bends.
!
!
!  Options 4, 5, and 6:
!
!    The starting point is (4, 3, 0).
!
!    The function F(X) has the form
!
!      F(1) = X(1) - X(2)**3 + 5*X(2)**2 -  2*X(2) - 13 +  3*(X(3)-1)
!      F(2) = X(1) + X(2)**3 +   X(2)**2 - 14*X(2) - 29 - 31*(X(3)-1)
!
!    There is a closed form representation of the curve in terms of the
!    second parameter:
!
!      X(1) = (14*X(2)**3 -79*X(2)**2 +52*X(2) + 245) / 17
!      X(2) = X(2)
!      X(3) = (   X(2)**3 - 2*X(2)**2 - 6*X(2) +   9) / 17
!
!    The correct value of the solution at X(3)=1 is:
!
!      (5, 4, 1)
!
!    In option 5, limit points in the first component are sought,
!    and in option 6, limit points in the third component are
!    sought.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ferdinand Freudenstein, Bernhard Roth,
!    Numerical Solutions of Nonlinear Equations,
!    Journal of the Association for Computing Machinery,
!    Volume 10, 1963, Pages 550-556.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar-1)
  real ( kind = 8 ) gx(2)
  real ( kind = 8 ) gy(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) y(3)
!
!  Get the starting point, Y.
!
  call p01_start ( option, nvar, y )
!
!  G is the function value at the starting point,
!  F the function value at the current point.
!
  call p01_gx ( y, gy )

  call p01_gx ( x, gx )
!
!  The parameter X3 generates the homotopy curve.
!
  fx(1:nvar-1) = gx(1:nvar-1) + ( x(3) - 1.0D+00 ) * gy(1:nvar-1)

  return
end
subroutine p01_gx ( x, g )

!*****************************************************************************80
!
!! P01_GX evaluates the underlying function for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(2), the point at which the function is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) G(2), the value of the function at X.
!
  implicit none

  real ( kind = 8 ) g(2)
  real ( kind = 8 ) x(2)

  g(1) = x(1) - ( ( x(2) - 5.0D+00 ) * x(2) + 2.0D+00 ) * x(2) - 13.0D+00
  g(2) = x(1) + ( ( x(2) + 1.0D+00 ) * x(2) - 14.0D+00 ) * x(2) - 29.0D+00

  return
end
subroutine p01_jac ( option, nvar, x, jac )

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
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) gy(3)
  integer ( kind = 4 ) option
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) y(3)

  jac(1:nvar,1:nvar) = 0.0D+00

  jac(1,1) = 1.0D+00
  jac(2,1) = 1.0D+00

  jac(1,2) = ( - 3.0D+00 * x(2) + 10.0D+00 ) * x(2) - 2.0D+00
  jac(2,2) = ( 3.0D+00 * x(2) + 2.0D+00 ) * x(2) - 14.0D+00
!
!  Get the starting point
!
  call p01_start ( option, nvar, y )
!
!  Get the function value at the starting point
!
  call p01_gx ( y, gy )

  jac(1,3) = gy(1)
  jac(2,3) = gy(2)

  return
end
subroutine p01_nvar ( option, nvar )

!*****************************************************************************80
!
!! P01_NVAR sets the number of variables for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option chosen for this problem.
!    For some problems, several options are available.  At least,
!    OPTION = 1 is always legal.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 3

  return
end
subroutine p01_option_num ( option_num )

!*****************************************************************************80
!
!! P01_OPTION_NUM returns the number of options for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 6

  return
end
subroutine p01_start ( option, nvar, x )

!*****************************************************************************80
!
!! P01_START returns a starting point for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  if ( option == 1 .or. option == 2 .or. option == 3 ) then
    x(1:3) = (/ 15.0D+00, -2.0D+00, 0.0D+00 /)
  else if ( option == 4 .or. option == 5 .or. option == 6 ) then
    x(1:3) = (/ 4.0D+00, 3.0D+00, 0.0D+00 /)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P01_START - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
subroutine p01_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P01_STEPSIZE returns step sizes for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    0.30000D+00
  hmin = 0.03125D+00
  hmax = 4.00000D+00

  return
end
subroutine p01_title ( option, title )

!*****************************************************************************80
!
!! P01_TITLE sets the title for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  if ( option == 1 ) then
    title = 'Freudenstein-Roth function, (15,-2,0).'
  else if ( option == 2 ) then
    title = 'Freudenstein-Roth function, (15,-2,0), x1 limits.'
  else if ( option == 3 ) then
    title = 'Freudenstein-Roth function, (15,-2,0), x3 limits.'
  else if ( option == 4 ) then
    title = 'Freudenstein-Roth function, (4,3,0).'
  else if ( option == 5 ) then
    title = 'Freudenstein-Roth function, (4,3,0), x1 limits.'
  else if ( option == 6 ) then
    title = 'Freudenstein-Roth function, (4,3,0), x3 limits.'
  else
    title = '???'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P01_TITLE - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized option number.'
    stop
  end if

  return
end
subroutine p02_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P02_FUN evaluates the function for problem 2.
!
!  Title:
!
!    The Boggs function
!
!  Description:
!
!    The function F is derived via homotopy from a simpler function:
!
!      F(X(1),X(2),X(3)) = G(X(1),X(2)) + (X(3)-1) * G(Y1,Y2)
!
!    with
!
!      (Y1, Y2) some starting value,
!
!    and
!
!      G(1) = X(1)*X(1) - X(2) + 1
!      G(2) = X(1) - COS(PI*X(2)/2)
!
!  Options:
!
!    OPTION = 1,
!      use starting point (  1,  0, 0 ).
!    OPTION = 2,
!      use starting point (  1, -1, 0 ).
!    OPTION = 3,
!      use starting point ( 10, 10, 0 ).
!
!  Target Points:
!
!    For the target value X(3) = 1.0, the solution is ( 0, 1, 1 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Boggs,
!    The Solution of Nonlinear Systems by A-stable Integration Techniques,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 4, December 1971, pages 767-785.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar-1)
  real ( kind = 8 ) gx(2)
  real ( kind = 8 ) gy(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) y(3)
!
!  Get the starting point
!
  call p02_start ( option, nvar, y )
!
!  Get the function value at the starting point and at the
!  current point.
!
  call p02_gx ( y, gy )
  call p02_gx ( x, gx )
!
!  Use X3 to compute a homotopy.
!
  do i = 1, nvar - 1
    fx(i) = gx(i) + ( x(3) - 1.0D+00 ) * gy(i)
  end do

  return
end
subroutine p02_gx ( x, g )

!*****************************************************************************80
!
!! P02_GX evaluates the underlying function for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(2), the point at which the function is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) G(2), the value of the function at X.
!
  implicit none

  real ( kind = 8 ) g(2)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(2)

  g(1) = x(1) * x(1) - x(2) + 1.0D+00
  g(2) = x(1) - cos ( pi * x(2) / 2.0D+00 )

  return
end
subroutine p02_jac ( option, nvar, x, jac )

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
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) gy(2)
  integer ( kind = 4 ) option
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) y(3)

  jac(1:nvar,1:nvar) = 0.0D+00

  jac(1,1) = 2.0D+00 * x(1)
  jac(2,1) = 1.0D+00
  jac(1,2) = - 1.0D+00
  jac(2,2) = 0.5D+00 * pi * sin ( 0.5D+00 * pi * x(2) )
!
!  Get the starting point
!
  call p02_start ( option, nvar, y )
!
!  Get the function value at the starting point
!
  call p02_gx ( y, gy )

  jac(1,3) = gy(1)
  jac(2,3) = gy(2)

  return
end
subroutine p02_nvar ( option, nvar )

!*****************************************************************************80
!
!! P02_NVAR sets the number of variables for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 3

  return
end
subroutine p02_option_num ( option_num )

!*****************************************************************************80
!
!! P02_OPTION_NUM returns the number of options for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 3

  return
end
subroutine p02_start ( option, nvar, x )

!*****************************************************************************80
!
!! P02_START returns a starting point for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  if ( option == 1 ) then
    x(1:3) = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
  else if ( option == 2 ) then
    x(1:3) = (/ 1.0D+00, -1.0D+00, 0.0D+00 /)
  else if ( option == 3 ) then
    x(1:3) = (/ 10.0D+00, 10.0D+00, 0.0D+00 /)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P02_START - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
subroutine p02_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P02_STEPSIZE returns step sizes for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    0.250D+00
  hmin = 0.001D+00
  hmax = 1.000D+00

  return
end
subroutine p02_title ( option, title )

!*****************************************************************************80
!
!! P02_TITLE sets the title for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  if ( option == 1 ) then
    title = 'Boggs function, (1,0,0).'
  else if ( option == 2 ) then
    title = 'Boggs function, (1,-1,0).'
  else if ( option == 3 ) then
    title = 'Boggs function, (10,10,0).'
  else
    title = '???'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P02_TITLE - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized option number.'
    stop
  end if

  return
end
subroutine p03_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P03_FUN evaluates the function for problem 3.
!
!  Title:
!
!    The Powell function
!
!  Description:
!
!    The function F is derived via homotopy from a simpler function G:
!
!      F(X(1),X(2),X(3)) = G(X(1),X(2)) + (X(3)-1)*G(Y1,Y2)
!
!    with
!
!      Y1, Y2 some starting point,
!
!    and
!
!      G(1) = 10000 * X(1) * X(2) - 1.0D+00
!      G(2) = exp ( -X(1) ) + exp ( -X(2) ) - 1.0001
!
!  Options:
!
!    OPTION = 1,
!      use starting point ( 3, 6, 0 );
!    OPTION = 2,
!      use starting point ( 4, 5, 0 );
!    OPTION = 3,
!      use starting point ( 6, 3, 0 );
!    OPTION = 4,
!      use starting point ( 1, 1, 0 ).
!
!  Special points:
!
!    For all options, there is a solution with last component 1, whose
!    value is either:
!
!      (1.098159E-5, 9.106146, 1.0)
!    or
!      (9.106146, 1.098159E-5, 1.0)
!
!  Comments:
!
!    Note that the function G is symmetric in X(1) and X(2).  Hence,
!    the run with starting point (1,1,0) should be interesting!
!
!    It would be worthwhile to seek limit points in X(3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Michael Powell,
!    A FORTRAN Subroutine for Solving Systems of Nonlinear Algebraic Equations,
!    in Numerical Methods for Nonlinear Algebraic Equations,
!    Edited by Philip Rabinowitz,
!    Gordon and Breach, 1970,
!    ISBN13: 978-0677142302,
!    LC: QA218.N85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar-1)
  real ( kind = 8 ) gx(2)
  real ( kind = 8 ) gy(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) y(3)
!
!  Get the starting point
!
  call p03_start ( option, nvar, y )
!
!  Get the (underlying) function value at the starting point and at the
!  current point.
!
  call p03_gx ( y, gy )
  call p03_gx ( x, gx )
!
!  Use X3 to compute a homotopy.
!
  do i = 1, nvar - 1
    fx(i) = gx(i) + ( x(3) - 1.0D+00 ) * gy(i)
  end do

  return
end
subroutine p03_gx ( x, g )

!*****************************************************************************80
!
!! P03_GX evaluates the underlying function for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(2), the point at which the function is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) G(2), the value of the function at X.
!
  implicit none

  real ( kind = 8 ) g(2)
  real ( kind = 8 ) x(2)

  g(1) = 10000.0D+00 * x(1) * x(2) - 1.0D+00
  g(2) = exp ( - x(1) ) + exp ( - x(2) ) - 1.0001D+00

  return
end
subroutine p03_jac ( option, nvar, x, jac )

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
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) gy(2)
  integer ( kind = 4 ) option
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) y(3)

  jac(1:nvar,1:nvar) = 0.0D+00
!
!  Get the starting point.
!
  call p03_start ( option, nvar, y )
!
!  Get the (underlying) function value at the starting point.
!
  call p03_gx ( y, gy )
!
!  The last column of the jacobian depends on the (underlying) function
!  value at the starting point.
!
  jac(1,1) = 10000.0D+00 * x(2)
  jac(1,2) = 10000.0D+00 * x(1)
  jac(1,3) = gy(1)

  jac(2,1) = - exp ( - x(1) )
  jac(2,2) = - exp ( - x(2) )
  jac(2,3) = gy(2)

  return
end
subroutine p03_nvar ( option, nvar )

!*****************************************************************************80
!
!! P03_NVAR sets the number of variables for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 3

  return
end
subroutine p03_option_num ( option_num )

!*****************************************************************************80
!
!! P03_OPTION_NUM returns the number of options for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 4

  return
end
subroutine p03_start ( option, nvar, x )

!*****************************************************************************80
!
!! P03_START returns a starting point for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  if ( option == 1 ) then
    x(1:3) = (/ 3.0D+00, 6.0D+00, 0.0D+00 /)
  else if ( option == 2 ) then
    x(1:3) = (/ 5.0D+00, 4.0D+00, 0.0D+00 /)
  else if ( option == 3 ) then
    x(1:3) = (/ 6.0D+00, 3.0D+00, 0.0D+00 /)
  else if ( option == 4 ) then
    x(1:3) = (/ 1.0D+00, 1.0D+00, 0.0D+00 /)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P03_START - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
subroutine p03_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P03_STEPSIZE returns step sizes for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    0.50000D+00
  hmin = 0.00025D+00
  hmax = 3.00000D+00

  return
end
subroutine p03_title ( option, title )

!*****************************************************************************80
!
!! P03_TITLE sets the title for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  if ( option == 1 ) then
    title = 'Powell function, (3,6,0).'
  else if ( option == 2 ) then
    title = 'Powell function, (4,5,0).'
  else if ( option == 3 ) then
    title = 'Powell function, (6,3,0).'
  else if ( option == 4 ) then
    title = 'Powell function, (1,1,0).'
  else
    title = '???'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P03_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
subroutine p04_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P04_FUN evaluates the function for problem 4.
!
!  Title:
!
!    The Broyden function
!
!  Description:
!
!    The function F is derived via homotopy from a simpler function G:
!
!      F(X(1),X(2),X(3)) = g(X(1),X(2)) + (X(3)-1) * G(Y1,Y2).
!
!    with
!
!      (Y1,Y2) some starting point,
!
!    and
!
!      G(1) = 0.5*sin(X(1)*X(2)) - X(2)/PI - X(1)
!      G(2) = (1-1/(4*PI))*(exp(2*X(1))-E) + E*X(2)/PI- 2*E*X(1)
!
!    where "E" = exp(1).
!
!  Options:
!
!    The only option starts with (0.4, 3, 0), and seeks the target
!    solution whose third component is 1.  The correct value of the
!    target solution is
!
!      ( -0.2207014, 0.8207467, 1.0D+00 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Broyden,
!    A New Method of Solving Nonlinear Simultaneous Equations,
!    The Computer Journal,
!    Volume 12, 1969, pages 94-99.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar-1)
  real ( kind = 8 ) gx(2)
  real ( kind = 8 ) gy(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) y(3)
!
!  Get the starting point.
!
  call p04_start ( option, nvar, y )
!
!  Get the function value at the starting point and at the
!  current point.
!
  call p04_gx ( y, gy )
  call p04_gx ( x, gx )
!
!  Use X3 to compute a homotopy.
!
  do i = 1, nvar - 1
    fx(i) = gx(i) + ( x(3) - 1.0D+00 ) * gy(i)
  end do

  return
end
subroutine p04_gx ( x, g )

!*****************************************************************************80
!
!! P04_GX evaluates the underlying function for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(2), the point at which the function is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) G(2), the value of the function at X.
!
  implicit none

  real ( kind = 8 ), parameter :: E = 2.71828182845904523536D+00
  real ( kind = 8 ) g(2)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(2)

  g(1) = 0.5D+00 * sin ( x(1) * x(2) ) - x(2) / pi - x(1)

  g(2) = ( 4.0D+00 * pi - 1.0D+00 ) * ( exp ( 2.0D+00 * x(1) ) - E ) &
    / ( 4.0D+00 * pi ) + E * x(2) / pi - 2.0D+00 * E * x(1)

  return
end
subroutine p04_jac ( option, nvar, x, jac )

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
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ), parameter :: E = 2.71828182845904523536D+00
  real ( kind = 8 ) gy(2)
  integer ( kind = 4 ) option
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) y(3)

  jac(1:nvar,1:nvar) = 0.0D+00

  jac(1,1) = 0.5D+00 * x(2) * cos ( x(1) * x(2) ) - 1.0D+00
  jac(2,1) = ( 4.0D+00 * pi - 1.0D+00 ) &
    * 2.0D+00 * exp ( 2.0D+00 * x(1) ) &
    / ( 4.0D+00 * pi ) - 2.0D+00 * E
  jac(1,2) = 0.5D+00 * x(1) * cos ( x(1) * x(2) ) - 1.0D+00 / pi
  jac(2,2) = E / pi
!
!  Get the starting point
!
  call p04_start ( option, nvar, y )
!
!  Get the function value at the starting point
!
  call p04_gx ( y, gy )

  jac(1,3) = gy(1)
  jac(2,3) = gy(2)

  return
end
subroutine p04_nvar ( option, nvar )

!*****************************************************************************80
!
!! P04_NVAR sets the number of variables for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 3

  return
end
subroutine p04_option_num ( option_num )

!*****************************************************************************80
!
!! P04_OPTION_NUM returns the number of options for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 1

  return
end
subroutine p04_start ( option, nvar, x )

!*****************************************************************************80
!
!! P04_START returns a starting point for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1:3) = (/ 0.4D+00, 3.0D+00, 0.0D+00 /)

  return
end
subroutine p04_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P04_STEPSIZE returns step sizes for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =     0.300D+00
  hmin =  0.001D+00
  hmax = 25.000D+00

  return
end
subroutine p04_title ( option, title )

!*****************************************************************************80
!
!! P04_TITLE sets the title for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  title = 'Broyden function'

  return
end
subroutine p05_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P05_FUN evaluates the function for problem 5.
!
!  Title:
!
!    The Wacker function
!
!  Description:
!
!    The function is of the form
!
!      F(1) = ( 1 - A * X(4) ) * X(1) + X(4) * exp ( - X(2) ) / 3.0D+00
!           + X(4) * ( A - 1 - 1/(3*E) )
!
!      F(2) = ( 1 - A * X(4) ) * X(2) - X(4) * log ( 1 + X(3) * X(3) ) / 5
!           + X(4) * ( A - 1 - log(2)/5 )
!
!      F(3) = ( 1 - A * X(3) ) * X(3) + X(4) * sin ( X(1) )
!           + X(4) * ( A - 1 - sin(1) )
!
!    with
!
!      A is a parameter, and
!      E is the base of the natural logarithm system, EXP(1.0).
!
!  Starting Point:
!
!    ( 0, 0, 0, 0 )
!
!  Options:
!
!    OPTION = 1,
!      A = 0.1;
!    OPTION = 2,
!      A = 0.5;
!    OPTION = 3,
!      A = 1.0.
!
!  Special points:
!
!
!    The value of the solution for which X(3) is 1 depends on the option
!    chosen:
!
!    Option    X(1)       X(2)       X(3)      X(4)
!
!      1    ( 1.147009,  1.431931,  1.000000, 1.084425 ).
!      2    ( 0.2412182, 0.4558247, 1.000000, 0.4534797 ).
!      3    ( 0.0000000, 0.0000000, 1.000000, 0.000000 ).
!
!    For option 3, there is a limit point in variable X(4):
!
!      ( -0.07109918, 0.06921115, 0.5009694, 0.2739685 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hans-Joerg Wacker, Erich Zarzer, Werner Zulehner,
!    Optimal Stepsize Control for the Globalized Newton Method,
!    in Continuation Methods,
!    edited by Hans-Joerg Wacker,
!    Academic Press, 1978,
!    ISBN: 0127292500,
!    LC: QA1.S899.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) aval
  real ( kind = 8 ), parameter :: E = 2.71828182845904523536D+00
  real ( kind = 8 ) fx(nvar-1)
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  if ( option == 1 ) then
    aval = 0.1D+00
  else if ( option == 2 ) then
    aval = 0.5D+00
  else
    aval = 1.0D+00
  end if

  fx(1) = ( 1.0D+00 - aval * x(4) ) * x(1) &
        + x(4) * exp ( - x(2) ) / 3.0D+00 &
        + x(4) * ( aval - 1.0D+00 - 1.0D+00 / ( 3.0D+00 * E ) )

  fx(2) = ( 1.0D+00 - aval * x(4) ) * x(2) &
        - x(4) * log ( 1.0D+00 + x(3) * x(3) ) / 5.0D+00 &
        + x(4) * ( aval - 1.0D+00 - log ( 2.0D+00 ) / 5.0D+00 )

  fx(3) = ( 1.0D+00 - aval * x(3) ) * x(3) &
        + x(4) * sin ( x(1) ) &
        + x(4) * ( aval - 1.0D+00 - sin ( 1.0D+00 ) )

  return
end
subroutine p05_jac ( option, nvar, x, jac )

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
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) aval
  real ( kind = 8 ), parameter :: E = 2.71828182845904523536D+00
  integer ( kind = 4 ) option
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) x(nvar)

  jac(1:nvar,1:nvar) = 0.0D+00

  if ( option == 1 ) then
    aval = 0.1D+00
  else if ( option == 2 ) then
    aval = 0.5D+00
  else
    aval = 1.0D+00
  end if

  jac(1,1) = 1.0D+00 - aval * x(4)
  jac(1,2) = - x(4) * exp ( - x(2) ) / 3.0D+00
  jac(1,3) = 0.0D+00
  jac(1,4) = - aval * x(1) + exp ( - x(2) ) / 3.0D+00 &
           + ( aval - 1.0D+00 - 1.0D+00 / ( 3.0D+00 * E ) )

  jac(2,1) = 0.0D+00
  jac(2,2) = 1.0D+00 - aval * x(4)
  jac(2,3) = - 2.0D+00 * x(3) * x(4) / ( 5.0D+00 * ( 1.0D+00 + x(3) * x(3) ) )
  jac(2,4) = - aval * x(2) - log ( 1.0D+00 + x(3) * x(3) ) / 5.0D+00 &
           + ( aval - 1.0D+00 - log ( 2.0D+00 ) / 5.0D+00 )

  jac(3,1) = x(4) * cos ( x(1) )
  jac(3,2) = 0.0D+00
  jac(3,3) = 1.0D+00 - 2.0D+00 * aval * x(3)
  jac(3,4) = sin ( x(1) ) + ( aval - 1.0D+00 - sin ( 1.0D+00 ) )

  return
end
subroutine p05_nvar ( option, nvar )

!*****************************************************************************80
!
!! P05_NVAR sets the number of variables for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 4

  return
end
subroutine p05_option_num ( option_num )

!*****************************************************************************80
!
!! P05_OPTION_NUM returns the number of options for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 3

  return
end
subroutine p05_start ( option, nvar, x )

!*****************************************************************************80
!
!! P05_START returns a starting point for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1:4) = 0.0D+00

  return
end
subroutine p05_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P05_STEPSIZE returns step sizes for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =     0.300D+00
  hmin =  0.001D+00
  hmax = 25.000D+00

  return
end
subroutine p05_title ( option, title )

!*****************************************************************************80
!
!! P05_TITLE sets the title for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  if ( option == 1 ) then
    title = 'Wacker function, A = 0.1.'
  else if ( option == 2 ) then
    title = 'Wacker function, A = 0.5.'
  else if ( option == 3 ) then
    title = 'Wacker function, A = 1.0.'
  else
    title = '???'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P05_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
subroutine p06_barray ( b )

!*****************************************************************************80
!
!! P06_BARRAY sets the B array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) B(5,8), the array of coefficients for the linear
!    part of the aircraft stability function.
!
  implicit none

  real ( kind = 8 ) b(5,8)
  real ( kind = 8 ), parameter, dimension ( 5, 8 ) :: b_save = reshape ( (/ &
    -3.933D+00,  0.0D+00,    0.002D+00,  0.0D+00,    0.0D+00, &
     0.107D+00, -0.987D+00,  0.0D+00,    1.0D+00,    0.0D+00, &
     0.126D+00,  0.0D+00,   -0.235D+00,  0.0D+00,   -1.0D+00, &
     0.0D+00,  -22.95D+00,   0.0D+00,   -1.0D+00,    0.0D+00, &
    -9.99D+00,   0.0D+00,    5.67D+00,   0.0D+00,   -0.196D+00, &
     0.0D+00,  -28.37D+00,   0.0D+00,   -0.168D+00,  0.0D+00, &
   -45.83D+00,   0.0D+00,   -0.921D+00,  0.0D+00,   -0.0071D+00, &
    -7.64D+00,   0.0D+00,   -6.51D+00,   0.0D+00,    0.0D+00 /), (/ 5, 8 /) )

  b(1:5,1:8) = b_save(1:5,1:8)

  return
end
subroutine p06_carray ( c )

!*****************************************************************************80
!
!! P06_CARRAY sets the C array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) C(5,8,8), the array of coefficients for the nonlinear
!    part of the aircraft stability function.
!
  implicit none

  real ( kind = 8 ) c(5,8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  c(1:5,1:8,1:8) = 0.0D+00

  c(1,2,3) = -   0.727D+00
  c(1,3,4) =     8.39D+00
  c(1,4,5) = - 684.4D+00
  c(1,4,7) =  + 63.5D+00

  c(2,1,3) =   + 0.949D+00
  c(2,1,5) =   + 0.173D+00

  c(3,1,2) =   - 0.716D+00
  c(3,1,4) =   - 1.578D+00
  c(3,4,7) =   + 1.132D+00

  c(4,1,5) =   - 1.0D+00

  c(5,1,4) =   + 1.0D+00

  return
end
subroutine p06_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P06_FUN evaluates the function for problem 6.
!
!  Title:
!
!    The aircraft stability problem.
!
!  Description:
!
!    The equations describe the behavior of an aircraft under the
!    control of a pilot.  The variables are
!
!      X(1) = roll
!      X(2) = pitch
!      X(3) = yaw
!      X(4) = angle of attack
!      X(5) = sideslip
!      X(6) = elevator
!      X(7) = aileron
!      X(8) = rudder
!
!    The function is of the following form
!
!    For indices I=1 through 5,
!
!      F(I) = SUM ( 1 <= J <= 8 ) B(I,J) * X(J)
!           + SUM ( 1 <= J <= 8, 1 <= K <= 8 ) C(I,J,K) * X(J) * X(K)
!
!    with the last two equations fixing the values of the elevator
!    and rudder:
!
!      F(6) = X(6) - value
!      F(7) = X(8)
!
!    Note that in the paper by Melhem and Rheinboldt, there are two
!    mistakes in the description of the function PHI(Y,U).  In both
!    cases, the factor "Y4*Y2" should be replaced by "Y4*U2".
!
!  Options:
!
!    There are five options, which vary in the value they fix the
!    elevator value in function 6:
!
!      Option   Elevator Value    Limit Points in X(7)
!
!       1        -0.050              1
!       2        -0.008              3
!       3         0.0D+00            2
!       4         0.05               1
!       5         0.1                1
!
!  Special points:
!
!    Melhem and Rheinboldt list the following limit points in X(7)
!    (note that Melhem has B(4,1)=1.0, B(4,2)=0.0)
!
!       X(1)    X(2)     X(3)    X(4)    X(5)    X(6)   X(7)     X(8)
!
!    -2.9691  0.8307  -0.0727  0.4102 -0.2688 -0.05   0.5091   0.0
!
!    -2.8158 -0.1748  -0.0894  0.0263  0.0709 -0.008  0.2044   0.0
!    -3.7571 -0.6491  -0.3935  0.0918  0.1968 -0.008 -0.0038   0.0
!    -4.1637  0.0922  -0.0926  0.0224 -0.0171 -0.008  0.3782   0.0
!
!    -2.5839 -0.2212  -0.0540  0.0135  0.0908  0.0    0.1860   0.0
!    -3.9007 -1.1421  -0.5786  0.1328  0.3268  0.0   -0.5070   0.0
!
!    -2.3610 -0.7236   0.0327 -0.0391  0.2934  0.05   0.2927   0.0
!
!    -2.2982  1.4033   0.0632 -0.0793  0.5833  0.10   0.5833   0.0
!
!    Rheinboldt lists the following limit points in X(7), with
!    B(4,1)=0.0, B(4,2)=1.0:
!
!       X(1)    X(2)     X(3)    X(4)    X(5)    X(6)   X(7)     X(8)
!
!     2.9648  0.8255   0.0736  0.0413  0.2673 -0.050 -0.0504   0.0
!
!     2.8173 -0.1762   0.0899  0.0264 -0.0714 -0.008 -0.2049   0.0
!     3.7579 -0.6554   0.3865  0.0925 -0.1986 -0.008  0.0062   0.0
!     4.1638  0.0891   0.0948  0.0228  0.1623 -0.008 -0.3776   0.0
!
!     2.5873 -0.2235   0.0546  0.0136 -0.0916  0.000 -0.1869   0.0
!     3.9005 -1.1481   0.5815  0.1335 -0.3285  0.000  0.5101   0.0
!
!     2.3639 -0.7297  -0.3160 -0.0387 -0.2958  0.050 -0.2957   0.0
!
!     2.2992 -1.4102  -0.0618 -0.0790 -0.5862  0.100 -0.6897   0.0
!
!    Rheinboldt lists the following bifurcation points:
!
!       X(1)    X(2)     X(3)    X(4)    X(5)    X(6)    X(7)    X(8)
!
!     4.482   0.1632   0.0237  0.0062  0.0352 -0.0006 -0.3986  0.0
!     3.319  -0.1869   0.1605  0.0437 -0.0688 -0.0125 -0.2374  0.0
!     4.466   0.1467   0.0404  0.0097  0.0308 -0.0061 -0.3995  0.0
!    -3.325   0.1880  -0.1614  0.0439  0.0691 -0.0124  0.2367  0.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Raman Mehra, William Kessel, James Carroll,
!    Global stability and contral analysis of aircraft at high angles of attack,
!    Technical Report CR-215-248-1, -2, -3,
!    Office of Naval Research, June 1977.
!
!    Rami Melhem, Werner Rheinboldt,
!    A Comparison of Methods for Determining Turning Points of Nonlinear Equations,
!    Computing,
!    Volume 29, Number 3, September 1982, pages 201-226.
!
!    Albert Schy, Margery Hannah,
!    Prediction of Jump Phenomena in Roll-coupled Maneuvers of Airplanes,
!    Journal of Aircraft,
!    Volume 14, Number 4, 1977,  pages 375-382.
!
!    John Young, Albert Schy, Katherine Johnson,,
!    Prediction of Jump Phenomena in Aircraft Maneuvers, Including
!    Nonlinear Aerodynamic Effects,
!    Journal of Guidance and Control,
!    Volume 1, Number 1, 1978, pages 26-31.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) b(5,8)
  real ( kind = 8 ) c(5,8,8)
  real ( kind = 8 ) fx(nvar-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) val
  real ( kind = 8 ) x(nvar)
!
!  Compute the linear term.
!
  call p06_barray ( b )

  fx(1:5) = matmul ( b(1:5,1:8), x(1:8) )
!
!  Compute the nonlinear terms.
!
  call p06_carray ( c )

  do i = 1, 5
    do j = 1, 8
      do k = 1, 8
        fx(i) = fx(i) + c(i,j,k) * x(j) * x(k)
      end do
    end do
  end do
!
!  Set function values for two fixed variables.
!
  if ( option == 1 ) then
    val = - 0.050D+00
  else if ( option == 2 ) then
    val = - 0.008D+00
  else if ( option == 3 ) then
    val =   0.000D+00
  else if ( option == 4 ) then
    val =   0.050D+00
  else if ( option == 5 ) then
    val =   0.100D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P06_FUN - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  fx(6) = x(6) - val
  fx(7) = x(8)

  return
end
subroutine p06_jac ( option, nvar, x, jac )

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
!    13 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) b(5,8)
  real ( kind = 8 ) c(5,8,8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(nvar)

  jac(1:nvar,1:nvar) = 0.0D+00
!
!  Set the jacobian to the linear coefficients.
!
  call p06_barray ( b )

  jac(1:5,1:8) = b(1:5,1:8)
!
!  Add the nonlinear terms.
!
  call p06_carray ( c )

  do i = 1, 5
    do j = 1, 8
      do k = 1, 8
        jac(i,j) = jac(i,j) + ( c(i,j,k) + c(i,k,j) ) * x(k)
      end do
    end do
  end do
!
!  Constraint equations.
!
  jac(6,6) = 1.0D+00
  jac(7,8) = 1.0D+00

  return
end
subroutine p06_nvar ( option, nvar )

!*****************************************************************************80
!
!! P06_NVAR sets the number of variables for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 8

  return
end
subroutine p06_option_num ( option_num )

!*****************************************************************************80
!
!! P06_OPTION_NUM returns the number of options for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 5

  return
end
subroutine p06_start ( option, nvar, x )

!*****************************************************************************80
!
!! P06_START returns a starting point for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  if ( option == 1 ) then

    x(1:nvar) = (/ &
       1.06001162985175758D-03, &
       5.12061216467178115D-02, &
       5.79953409787390485D-05, &
       5.96060845777059631D-02, &
       2.64683802731226678D-05, &
      -5.00000000000000000D-02, &
       0.00000000000000000D+00, &
       0.00000000000000000D+00 /)

  else if ( option == 2 ) then

    x(1:nvar) = (/ &
       0.000001548268247D+00, &
       0.008192973225663D+00, &
      -0.000000682134573D+00, &
       0.009536973221178D+00, &
       0.000002896734870D+00, &
      -0.008000000000000D+00, &
       0.000018188778989D+00, &
       0.000000000000000D+00 /)

  else if ( option == 3 ) then

    x(1:nvar) = (/ &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00, &
       0.0D+00 /)

  else if ( option == 4 ) then

    x(1:nvar) = (/ &
      -0.000010655314069D+00, &
      -0.051206082422980D+00, &
       0.000005600187501D+00, &
      -0.059606082643400D+00, &
      -0.000020891016199D+00, &
       0.050000000000000D+00, &
      -0.000122595323216D+00, &
       0.000000000000000D+00 /)

  else if ( option == 5 ) then

    x(1:nvar) = (/ &
      -0.000027083319493D+00, &
      -0.102412164106124D+00, &
       0.000014540858026D+00, &
      -0.119212165322433D+00, &
      -0.000048014067202D+00, &
       0.100000000000000D+00, &
      -0.000267808407544D+00, &
       0.000000000000000D+00 /)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P06_START - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
subroutine p06_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P06_STEPSIZE returns step sizes for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =   - 0.250D+00
  hmin =  0.001D+00
  hmax =  0.500D+00

  return
end
subroutine p06_title ( option, title )

!*****************************************************************************80
!
!! P06_TITLE sets the title for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  if ( option == 1 ) then
    title = 'Aircraft function, x(6) = - 0.050.'
  else if ( option == 2 ) then
    title = 'Aircraft function, x(6) = - 0.008.'
  else if ( option == 3 ) then
    title = 'Aircraft function, x(6) =   0.000.'
  else if ( option == 4 ) then
    title = 'Aircraft function, x(6) = + 0.050.'
  else if ( option == 5 ) then
    title = 'Aircraft function, x(6) = + 0.100.'
  else
    title = '???'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P06_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
subroutine p07_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P07_FUN evaluates the function for problem 7.
!
!  Title:
!
!    Cell kinetics problem.
!
!  Description:
!
!    The function is of the form
!
!      F(I) = Sum ( 1 <= J <= NVAR-1)
!        A(I,J) * X(J) + RHO ( X(I) ) - X(NVAR)
!
!    with tridiagonal matrix A.
!
!  Special points:
!
!    Limit points in the variable NVAR are sought.  There are two:
!
!       X(1)      X(2)      X(3)      X(4)      X(5)       X(6)
!
!    ( 1.048362, 1.048362, 1.048362, 1.048362, 1.048362, 34.35693 ).
!    ( 8.822219, 8.822219, 8.822219, 8.822219, 8.822218, 18.88707 ).
!
!    There are also four bifurcation points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Werner Rheinboldt,
!    Solution Fields of Nonlinear Equations and Continuation Methods,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 221-237.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)
!
!  RHO(X) = 100.0D+00 * X / ( 1 + X + X * X )
!
  do i = 1, nvar - 1
    fx(i) = 100.0D+00 * x(i) / ( 1.0D+00 + x(i) + x(i) * x(i) ) - x(nvar)
  end do
!
!  The tridiagonal matrix A = (  2 -1  0  0  0  0 ... )
!                             ( -1  2 -1  0  0  0 ... )
!                             (  0 -1  2 -1  0  0 ... )
!
  fx(1) = fx(1) + 2.0D+00 * x(1) - x(2)
  do i = 2, nvar - 2
    fx(i) = fx(i) - x(i-1) + 3.0D+00 * x(i) - x(i+1)
  end do
  fx(nvar-1) = fx(nvar-1) - x(nvar-2) + 2.0D+00 * x(nvar-1)

  return
end
subroutine p07_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P07_JAC evaluates the jacobian for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) x(nvar)

  jac(1:nvar,1:nvar) = 0.0D+00

  do i = 1, nvar - 1
    jac(i,i) = 100.0D+00 * ( 1.0D+00 - x(i) * x(i) ) &
      / ( 1.0D+00 + x(i) + x(i) * x(i) )**2
  end do

  jac(1,1) = jac(1,1) + 2.0D+00
  jac(1,2) = jac(1,2) - 1.0D+00
  jac(1,nvar) = jac(1,nvar) - 1.0D+00

  do i = 2, nvar - 2
    jac(i,i-1) = jac(i,i-1) - 1.0D+00
    jac(i,i) = jac(i,i) + 3.0D+00
    jac(i,i+1) = jac(i,i+1) - 1.0D+00
    jac(i,nvar) = jac(i,nvar) - 1.0D+00
  end do

  jac(nvar-1,nvar-2) = jac(nvar-1,nvar-2) - 1.0D+00
  jac(nvar-1,nvar-1) = jac(nvar-1,nvar-1) + 2.0D+00
  jac(nvar-1,nvar) = jac(nvar-1,nvar) - 1.0D+00

  return
end
subroutine p07_nvar ( option, nvar )

!*****************************************************************************80
!
!! P07_NVAR sets the number of variables for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 6

  return
end
subroutine p07_option_num ( option_num )

!*****************************************************************************80
!
!! P07_OPTION_NUM returns the number of options for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 1

  return
end
subroutine p07_start ( option, nvar, x )

!*****************************************************************************80
!
!! P07_START returns a starting point for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1:nvar) = 0.0D+00

  return
end
subroutine p07_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P07_STEPSIZE returns step sizes for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    1.000D+00
  hmin = 0.001D+00
  hmax = 1.000D+00

  return
end
subroutine p07_title ( option, title )

!*****************************************************************************80
!
!! P07_TITLE sets the title for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  title = 'Cell kinetics problem, seeking limit points.'

  return
end
subroutine p08_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P08_FUN evaluates the function for problem 8.
!
!  Title:
!
!    Riks's mechanical problem.
!
!  Description:
!
!    The equations describe the equilibrium state of a structure made of
!    three springs with a common movable endpoint and the other
!    endpoints fixed.  A load is applied to the common endpoint.
!
!      X(1), X(2), and X(3) are the x, y, and z coordinates of the
!        common point.
!      X(4) is the magnitude of the load which is applied in the
!        X direction.
!
!    If C(I) is the spring constant for the I-th spring, and A(I,J)
!    is the J-th coordinate of the I-th fixed endpoint, then the
!    equation is:
!
!      F(J) = SUM(I=1,3) COEF(I)*(A(I,J)-X(J)) + P(J)
!
!    where
!
!      COEF(I) = C(I) * (NORM(A(I,*)-NORM(X-A(I,*))) / NORM(X-A(I,*) )
!
!    and
!
!      P=(X(4),X(5),X(6)) is an applied load, and
!
!      NORM(X) is the euclidean norm, and
!
!      c(1) + c(2) + c(3) = 1.0D+00
!
!    Two augmenting equations control the load vector P:
!
!      F(4) = X(ival1) - val1.
!      F(5) = X(ival2) - val2.
!
!    For this example,
!
!      ival1=4, val1=0
!      ival2=5, val2=0
!
!    and hence the load is all in the Z direction.
!
!    We seek limit points in X(6).
!
!    In Riks's paper, there seem to be limit points in X(6) at 4.10 and
!    -3.84.  The current code does not confirm this.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    E Riks,
!    The Application of Newton's Method to the Problem of Elastic Stability,
!    Transactions of the ASME, Journal of Applied Mechanics,
!    December 1972, pages 1060-1065.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) anrm
  real ( kind = 8 ) aval(3,3)
  real ( kind = 8 ) cval(3)
  real ( kind = 8 ) fx(nvar-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) ival1
  integer ( kind = 4 ) ival2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) val1
  real ( kind = 8 ) val2
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xmanrm

  call p08_gx ( aval, cval )

  do i = 1, 3

    fx(i) = 0.0D+00

    do j = 1, 3
!
!  Compute norms.
!
      anrm = 0.0D+00
      xmanrm = 0.0D+00

      do k = 1, 3
        anrm = anrm + aval(j,k)**2
        xmanrm = xmanrm + ( x(k) - aval(j,k) )**2
      end do

      anrm = sqrt ( anrm )
      xmanrm = sqrt ( xmanrm )

      fx(i) = fx(i) + cval(j) * ( 1.0D+00 - anrm / xmanrm ) &
        * ( x(i) - aval(j,i) )

    end do

  end do
!
!  Add the load vector: ( X(4), X(5), X(6) ).
!
  do i = 1, 3
    fx(i) = fx(i) + x(i+3)
  end do
!
!  Get constraints.
!
  call p08_hx ( option, ival1, ival2, val1, val2 )

  fx(4) = x(ival1) - val1
  fx(5) = x(ival2) - val2

  return
end
subroutine p08_gx ( aval, cval )

!*****************************************************************************80
!
!! P08_GX sets data used for Rik's mechanical problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) AVAL(3,3); for each I, the values of AVAL(I,*)
!    record the (X,Y,Z) coordinates of the I-th support point.
!
!    Output, real ( kind = 8 ) CVAL(3), the values of the normalized spring
!    constants.
!
  implicit none

  real ( kind = 8 ) aval(3,3)
  real ( kind = 8 ) cval(3)

  aval(1,1) = 2.0D+00
  aval(1,2) = 0.0D+00
  aval(1,3) = 0.0D+00

  aval(2,1) = - 1.0D+00
  aval(2,2) =   1.0D+00
  aval(2,3) =   0.0D+00

  aval(3,1) = - 1.0D+00
  aval(3,2) = - 2.0D+00
  aval(3,3) =   1.0D+00

  cval(1) = 10.0D+00 / 21.0D+00
  cval(2) =  6.0D+00 / 21.0D+00
  cval(3) =  5.0D+00 / 21.0D+00

  return
end
subroutine p08_hx ( option, ival1, ival2,  val1, val2 )

!*****************************************************************************80
!
!! P08_HX reports the constraint equation data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) IVAL1, IVAL2, the indices of the two
!    constrained variables.
!
!    Output, real ( kind = 8 ) VAL1, VAL2, the values to which the two
!    constrained variables are to be set.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) ival1
  integer ( kind = 4 ) ival2
  real ( kind = 8 ) val1
  real ( kind = 8 ) val2

  ival1 = 4
  val1 = 0.0D+00
  ival2 = 5
  val2 = 0.0D+00

  return
end
subroutine p08_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P08_JAC evaluates the jacobian for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) anrm
  real ( kind = 8 ) aval(3,3)
  real ( kind = 8 ) cval(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) ival1
  integer ( kind = 4 ) ival2
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) val1
  real ( kind = 8 ) val2
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xmanrm

  jac(1:nvar,1:nvar) = 0.0D+00

  call p08_gx ( aval, cval )

  do i = 1, 3
    do j = 1, 3
      do k = 1, 3
!
!  Compute norms.
!
        anrm = 0.0D+00
        xmanrm = 0.0D+00

        do l = 1, 3
          anrm = anrm + aval(k,l)**2
          xmanrm = xmanrm + ( x(l) - aval(k,l) )**2
        end do

        anrm = sqrt ( anrm )
        xmanrm = sqrt ( xmanrm )

        jac(i,j) = jac(i,j) + cval(k) * anrm * ( x(i) - aval(k,i) ) &
                 * ( x(j) - aval(k,j) ) / xmanrm**3

        if ( i == j ) then
          jac(i,j) = jac(i,j) - cval(k) * anrm / xmanrm
        end if

      end do
    end do
  end do

  do i = 1, 3
    jac(i,i) = jac(i,i) + 1.0D+00
  end do
!
!  Add the loads.
!
  jac(1,4) = 1.0D+00
  jac(2,5) = 1.0D+00
  jac(3,6) = 1.0D+00
!
!  Apply the constraints.
!
  call p08_hx ( option, ival1, ival2, val1, val2 )

  jac(4,ival1) = 1.0D+00
  jac(5,ival2) = 1.0D+00

  return
end
subroutine p08_nvar ( option, nvar )

!*****************************************************************************80
!
!! P08_NVAR sets the number of variables for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 6

  return
end
subroutine p08_option_num ( option_num )

!*****************************************************************************80
!
!! P08_OPTION_NUM returns the number of options for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 1

  return
end
subroutine p08_start ( option, nvar, x )

!*****************************************************************************80
!
!! P08_START returns a starting point for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) ival1
  integer ( kind = 4 ) ival2
  real ( kind = 8 ) val1
  real ( kind = 8 ) val2
  real ( kind = 8 ) x(nvar)

  x(1:nvar) = 0.0D+00

  call p08_hx ( option, ival1, ival2, val1, val2 )

  x(ival1) = val1
  x(ival2) = val2

  return
end
subroutine p08_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P08_STEPSIZE returns step sizes for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    1.000D+00
  hmin = 0.001D+00
  hmax = 1.000D+00

  return
end
subroutine p08_title ( option, title )

!*****************************************************************************80
!
!! P08_TITLE sets the title for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  title = 'Riks mechanical problem, seeking limit points.'

  return
end
subroutine p09_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P09_FUN evaluates the function for problem 9.
!
!  Title:
!
!    Oden mechanical problem.
!
!  Description:
!
!    The equations describe the equilibrium of a simple two bar
!    framework, with one common endpoint, and the other endpoints
!    fixed.  A load is applied to the common endpoint.  The bars are
!    constructed of an isotropic hookean material.
!
!    The function is of the form
!
!      F(1) = X(1)**3 - 3*height*X(1)**2 + 2*height**2*X(1)
!            +X(1)*X(2)**2 - height*X(2)**2 - X(3)*cos(X(4))
!
!      F(2) = X(1)**2*X(2) - 2*height*X(1)*X(2) + X(2)**3 + 2*X(2)
!            -X(3)*sin(X(4))
!
!      F(3) = X(IVAL) - VAL
!
!    with
!
!      HEIGHT=2.0D+00
!      IVAL=4
!      VAL varying, depending on the option
!
!  Options:
!
!        VAL   IT XIT  LIM
!
!     1  0.00, 1, 4.0, 1
!     2  0.25, 1, 4.0, 1
!     3  0.50, 1, 4.0, 1
!     4  1.00, 1, 4.0, 1
!     5  0.00, 1, 4.0, 2
!     6  0.25, 1, 4.0, 2
!     7  0.50, 1, 4.0, 2
!     8  1.00, 1, 4.0, 2
!     9  0.00, 1, 4.0, 3
!    10  0.25, 1, 4.0, 3
!    11  0.50, 1, 4.0, 3
!    12  1.00, 1, 4.0, 3
!    13  0.00, 0,      0.
!
!    For options 1, 5, and 9, the target point is (4,0,0,0).
!
!    For option 9, there are the following limit points in X(3)
!
!      (2+-2/sqrt(3), 0, +-16/sqrt(27), 0)
!
!    For skew loads (X(4) nonzero) there are various limit points.
!
!    Melhem lists,
!
!      (0.5903206, 0.8391448, 0.9581753, 1.252346)
!      (2.705446,  0.6177675, 0.9581753, 1.252346)
!
!    with X(3),X(4) corresponding to a load vector of (.30,.91).
!
!    Computational results with this program are:
!
!    OPTION = 2  limit points in X(1)
!
!    2.816913  0.7396444  -2.348587  0.2500000
!    1.183087 -0.7396445   2.348587  0.2500000
!
!    OPTION=3  limit points in X(1)
!
!    2.520900  0.8598542  -1.774344  0.5000000
!    1.479100 -0.8598521   1.774346  0.5000000
!
!    OPTION=4  limit points in X(1)
!
!    2.210747  0.9241686  -1.209751  1.0000000
!    (limit point finder failed at second limit point)
!
!    OPTION=6  limit points in X(2)
!
!    1.831179  1.424861  0.3392428  0.2500000
!    (apparently did not reach second limit point)
!
!    OPTION=7  limit points in X(2)
!
!    1.697061  1.453503  0.6198216  0.2500000
!    2.302939 -1.453503 -0.6198219  0.2500000
!
!    OPTION=8  limit points in X(2)
!
!    1.534293  1.555364  1.175649  1.0000000
!    2.465706 -1.555364 -1.175648  1.0000000
!
!    OPTION=9  limit points in X(3)
!
!    0.8452995  0.0000000  3.079199  0.0000000
!    3.154701   0.0000000 -3.079197  0.0000000
!
!    OPTION=10  limit points in X(3)
!
!    0.5800046  0.7846684  2.004746  0.2500000
!    2.777765   0.5695726 -2.464886  0.2500000
!
!    OPTION=11  limit points in X(3)
!
!    0.6305253  0.9921379  1.779294  0.5000000
!    2.501894   0.7202593 -1.846869  0.5000000
!
!    OPTION=12  limit points in X(3)
!
!    0.7650624  1.292679   1.837450  1.000000
!    2.204188   0.8010838 -1.253382  1.000000
!
!    Bifurcation points occur at
!
!    (2+-sqrt(2), 0, +-sqrt(2), 0)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Rami Melhem, Werner Rheinboldt,
!    A Comparison of Methods for Determining Turning Points of Nonlinear Equations,
!    Computing,
!    Volume 29, Number 3, September 1982, pages 201-226.
!
!    John Oden,
!    Finite Elements of Nonlinear Continua,
!    Dover, 2006,
!    ISBN: 0486449734,
!    LC: QA808.2.O33.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar-1)
  real ( kind = 8 ) height
  integer ( kind = 4 ) option
  integer ( kind = 4 ) ival
  real ( kind = 8 ) val
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4

  call p09_gx ( option, height, ival, val )

  x1 = x(1)
  x2 = x(2)
  x3 = x(3)
  x4 = x(4)

  fx(1) = x1**3 - 3.0D+00 * height * x1 * x1 &
    + 2.0D+00 * height * height * x1 &
    + x1 * x2 * x2 - height * x2 * x2 - x3 * cos ( x4 )

  fx(2) = x1 * x1 * x2 - 2.0D+00 * height * x1 * x2 + x2**3 &
        + 2.0D+00 * x2 - x3 * sin ( x4 )

  fx(3) = x(ival) - val

  return
end
subroutine p09_gx ( option, height, ival, val )

!*****************************************************************************80
!
!! P09_GX is used by problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) HEIGHT, the height of the structure.
!
!    Output, integer ( kind = 4 ) IVAL, the index of the variable being fixed.
!
!    Output, real ( kind = 8 ) VAL, the value of the fixed variable.
!
  implicit none

  real ( kind = 8 ) height
  integer ( kind = 4 ) option
  integer ( kind = 4 ) ival
  real ( kind = 8 ) val

  height = 2.0D+00
  ival = 4

  if ( option == 1 .or. option == 5 .or. option == 9 ) then
    val = 0.00D+00
  else if ( option == 2 .or. option == 6 .or. option == 10 ) then
    val = 0.25D+00
  else if ( option == 3 .or. option == 7 .or. option == 11 ) then
    val = 0.50D+00
  else if ( option == 4 .or. option == 8 .or. option == 12 ) then
    val = 1.00D+00
  else if ( option == 13 ) then
    val = 0.00D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P09_GX - Fatal error'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
subroutine p09_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P09_JAC evaluates the jacobian for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) height
  integer ( kind = 4 ) option
  integer ( kind = 4 ) ival
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) val
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4

  jac(1:nvar,1:nvar) = 0.0D+00

  call p09_gx ( option, height, ival, val )

  x1 = x(1)
  x2 = x(2)
  x3 = x(3)
  x4 = x(4)

  jac(1,1) = 3.0D+00 * x1 * x1 - 6.0D+00 * height * x1 &
    + 2.0D+00 * height * height + x2 * x2
  jac(1,2) = 2.0D+00 * x1 * x2 - 2.0D+00 * height * x2
  jac(1,3) = - cos ( x4 )
  jac(1,4) = x3 * sin ( x4 )

  jac(2,1) = 2.0D+00 * x1 * x2 - 2.0D+00 * height * x2
  jac(2,2) = x1 * x1 - 2.0D+00 * height * x1 + 3.0D+00 * x2 * x2 + 2.0D+00
  jac(2,3) = - sin ( x4 )
  jac(2,4) = - x3 * cos ( x4 )

  jac(3,1) = 0.0D+00
  jac(3,2) = 0.0D+00
  jac(3,3) = 0.0D+00
  jac(3,4) = 0.0D+00

  jac(3,ival) = 1.0D+00

  return
end
subroutine p09_nvar ( option, nvar )

!*****************************************************************************80
!
!! P09_NVAR sets the number of variables for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 4

  return
end
subroutine p09_option_num ( option_num )

!*****************************************************************************80
!
!! P09_OPTION_NUM returns the number of options for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 13

  return
end
subroutine p09_start ( option, nvar, x )

!*****************************************************************************80
!
!! P09_START returns a starting point for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) height
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) ival
  real ( kind = 8 ) val
  real ( kind = 8 ) x(nvar)

  x(1:nvar) = 0.0D+00

  call p09_gx ( option, height, ival, val )

  x(ival) = val

  return
end
subroutine p09_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P09_STEPSIZE returns step sizes for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    0.300D+00
  hmin = 0.001D+00
  hmax = 0.600D+00

  return
end
subroutine p09_title ( option, title )

!*****************************************************************************80
!
!! P09_TITLE sets the title for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  if ( option == 1 ) then
    title = 'Oden problem, VAL=0.00, Target X(1)=4.0, Limits in X(1).'
  else if ( option == 2 ) then
    title = 'Oden problem, VAL=0.25, Target X(1)=4.0, Limits in X(1).'
  else if ( option == 3 ) then
    title = 'Oden problem, VAL=0.50, Target X(1)=4.0, Limits in X(1).'
  else if ( option == 4 ) then
    title = 'Oden problem, VAL=1.00, Target X(1)=4.0, Limits in X(1).'
  else if ( option == 5 ) then
    title = 'Oden problem, VAL=0.00, Target X(1)=4.0, Limits in X(2).'
  else if ( option == 6 ) then
    title = 'Oden problem, VAL=0.25, Target X(1)=4.0, Limits in X(2).'
  else if ( option == 7 ) then
    title = 'Oden problem, VAL=0.50, Target X(1)=4.0, Limits in X(2).'
  else if ( option == 8 ) then
    title = 'Oden problem, VAL=1.00, Target X(1)=4.0, Limits in X(2).'
  else if ( option == 9 ) then
    title = 'Oden problem, VAL=0.00, Target X(1)=4.0, Limits in X(3).'
  else if ( option == 10 ) then
    title = 'Oden problem, VAL=0.25, Target X(1)=4.0, Limits in X(3).'
  else if ( option == 11 ) then
    title = 'Oden problem, VAL=0.50, Target X(1)=4.0, Limits in X(3).'
  else if ( option == 12 ) then
    title = 'Oden problem, VAL=1.00, Target X(1)=4.0, Limits in X(3).'
  else if ( option == 13 ) then
    title = 'Oden problem, VAL=0.00, no targets, no limits.'
  else
    title = '???'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P09_TITLE - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized option number.'
    stop
  end if

  return
end
subroutine p10_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P10_FUN evaluates the function for problem 10.
!
!  Title:
!
!    Torsion of a square rod, finite difference solution
!
!  Description:
!
!    The problem is a boundary value problem on (0,1) x (0,1)
!    of the form:
!
!      - d/dx ( PHI ( dU/dx, dU/dy ) * dU/dx )
!      - d/dy ( PHI ( dU/dx, dU/dy ) * dU/dy ) = G ( U, LAMBDA )
!
!    A standard finite difference approximation on a uniform mesh is
!    applied to yield the equations, with X(1) through X(NVAR-1) storing
!    the value of U at the mesh points, and X(NVAR) holding the value
!    of LAMBDA.
!
!  Options:
!
!    Let S = dU/dX**2 + dU/dY**2.
!
!    OPTION=1
!
!      PHI(S) = exp ( 5 * S )
!
!    OPTION=2
!
!      Let SBAR = ( 40 * S - 13 ) / 7
!
!      if ( S <= 0.15 ) then
!        PHI = 1.0D+00
!      else if ( 0.15 <= S <= 0.50 ) then
!        PHI = 5.5 + 4.5 * ( 3 * SBAR - SBAR**3 )
!      else if ( 0.50 <= S ) then
!        PHI = 10.0D+00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Werner Rheinboldt,
!    On the Solution of Some Nonlinear Equations Arising in the
!    Application of Finite Element Methods,
!    in The Mathematics of Finite Elements and Applications II,
!    edited by John Whiteman
!    Academic Press, London, 1976, pages 465-482,
!    LC: TA347.F5.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) cx(2,4)
  real ( kind = 8 ) fx(nvar-1)
  real ( kind = 8 ) g
  real ( kind = 8 ) gp
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) rlk
  real ( kind = 8 ) sc
  real ( kind = 8 ) uc(2)
  real ( kind = 8 ) ux(4)
  real ( kind = 8 ) x(nvar)

  nrow = 6
  ncol = 6

  h = 1.0D+00 / real ( nrow + 1, kind = 8 )

  do i = 1, nrow

    do j = 1, ncol

      ij = ( j - 1 ) * nrow + i
!
!  UC contains the two cornerpoints,
!
      if ( i == 1 .or. j == 1 ) then
        uc(1) = 0.0D+00
      else
        jk = ij - nrow
        uc(1) = x(jk-1)
      end if

      if ( j == ncol .or. i == nrow ) then
        uc(2) = 0.0D+00
      else
        jk = ij + nrow
        uc(2) = x(jk+1)
      end if
!
!  UX contains the four side-points,
!
      if ( i == 1 ) then
        ux(1) = 0.0D+00
      else
        ux(1) = x(ij-1)
      end if

      if ( i < nrow ) then
        ux(2) = x(ij+1)
      else
        ux(2) = 0.0D+00
      end if

      if ( j == 1 ) then
        ux(3) = 0.0D+00
      else
        jk = ij - nrow
        ux(3) = x(jk)
      end if

      if ( j < ncol ) then
        jk = ij + nrow
        ux(4) = x(jk)
      else
        ux(4) = 0.0D+00
      end if
!
!  CX contains the elements connected to the side points.
!
!  k = 1, 2*qw calculated and stored in ( cx(1,k) + cx(2,k) )
!  k = 2, 2*qe calculated and stored in ( cx(1,k) + cx(2,k) )
!  k = 3, 2*qs calculated and stored in ( cx(1,k) + cx(2,k) )
!  k = 4, 2*qn calculated and stored in ( cx(1,k) + cx(2,k) )
!
      sc = 0.0D+00

      do k = 1, 4

        if ( k == 1 .or. k == 3 ) then
          k1 = 1
          k2 = 4
        else
          k1 = 2
          k2 = 3
        end if

        do l = 1, 2

          rlk = ( ux(k) - x(ij) )**2

          if ( l == 1 ) then
            rlk = ( rlk + ( ux(k) - uc(k1) )**2 ) / h / h
          else
            rlk = ( rlk + ( x(ij) - ux(k2) )**2 ) / h / h
          end if

          call p10_gx ( option, rlk, g, gp )
          cx(l,k) = g
          sc = sc + cx(l,k)

        end do

      end do
!
!  sc = qn + qs + qe + qw
!
      fx(ij) = 0.5D+00 * sc * x(ij) - x(nvar) * h * h
      do k = 1, 4
        fx(ij) = fx(ij) - 0.5D+00 * ux(k) * ( cx(1,k) + cx(2,k) )
      end do

    end do

  end do

  return
end
subroutine p10_gx ( option, s, g, gp )

!*****************************************************************************80
!
!! P10_GX is used by problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, real ( kind = 8 ) S, the value of the argument of G.
!    S = (d U/d X)**2 + (d U/d Y)**2.
!
!    Output, real ( kind = 8 ) G, GP, the value of G(S) and d G(S)/d S.
!
  implicit none

  real ( kind = 8 ) g
  real ( kind = 8 ) gp
  integer ( kind = 4 ) option
  real ( kind = 8 ) s
  real ( kind = 8 ) sbar

  if ( option == 1 ) then

    g = exp ( 5.0D+00 * s )
    gp = 5.0D+00 * exp ( 5.0D+00 * s )

  else if ( option == 2 ) then

    if ( s <= 0.15D+00 ) then

      g = 1.0D+00
      gp = 0.0D+00

    else if ( 0.15D+00 < s .and. s < 0.5D+00 ) then

      sbar = ( 40.0D+00 * s - 13.0D+00 ) / 7.0D+00

      g = 5.5D+00 + 2.25D+00 * sbar * ( 3.0D+00 - sbar * sbar )

      gp = 270.0D+00 * ( 1.0D+00 - sbar * sbar ) / 7.0D+00

    else

      g = 10.0D+00
      gp = 0.0D+00

    end if

  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P10_GX - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
subroutine p10_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P10_JAC evaluates the jacobian for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) cx(2,4)
  real ( kind = 8 ) dx(2,4)
  real ( kind = 8 ) g
  real ( kind = 8 ) gp
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) rlk
  real ( kind = 8 ) sc
  real ( kind = 8 ) uc(2)
  real ( kind = 8 ) ux(4)
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xjac

  jac(1:nvar,1:nvar) = 0.0D+00

  nrow = 6
  ncol = 6

  h = 1.0D+00 / real ( nrow + 1, kind = 8 )

  do i = 1, nrow

    do j = 1, ncol

      uc(1) = 0.0D+00
      uc(2) = 0.0D+00
      ux(1:4) = 0.0D+00

      ij = i + ( j - 1 ) * nrow

      if ( i /= 1 ) then
        ux(1) = x(ij-1)
      end if

      if ( i /= nrow ) then
        ux(2) = x(ij+1)
      end if

      if ( 1 < j ) then
        jk = ij - nrow
        ux(3) = x(jk)
        if ( i /= 1 ) then
          uc(1) = x(jk-1)
        end if
      end if

      if ( j < ncol ) then
        jk = ij + nrow
        ux(4) = x(jk)
        if ( i /= nrow ) then
          uc(2) = x(jk+1)
        end if
      end if

      sc = 0.0D+00

      do k = 1, 4

        if ( k == 1 .or. k == 3 ) then
          k1 = 1
        else
          k1 = 2
        end if
        k2 = 5 - k

        do l = 1, 2

          rlk = ( ux(k) - x(ij) )**2

          if ( l == 1 ) then
            rlk = ( rlk + ( ux(k) - uc(k1) )**2 ) / h / h
          else
            rlk = ( rlk + ( x(ij) - ux(k2) )**2 ) / h / h
          end if

          call p10_gx ( option, rlk, g, gp )
          cx(l,k) = g
          dx(l,k) = gp
          sc = sc + cx(l,k)
        end do

      end do
!
!  diagonal
!
      xjac = 0.5D+00 * sc

      do k = 1, 4
        k2 = 5 - k
        xjac = xjac + dx(2,k) * ( x(ij) - ux(k) ) &
              * ( 2.0D+00 * x(ij) - ux(k) - ux(k2) ) / h / h
        xjac = xjac + dx(1,k) * ( x(ij) - ux(k) )**2 / h / h
      end do

      jac(ij,ij) = xjac
!
!  off-diagonals
!
      do k = 1, 4

        if ( k == 1 ) then
          if ( i == 1 ) then
            continue
          end if
          jk = ij - 1
        else if ( k == 2 ) then
          if ( i == nrow ) then
            continue
          end if
          jk = ij + 1
        else if ( k == 3 ) then
          if ( j == 1 ) then
            continue
          end if
          jk = ij - nrow
        else if ( k == 4 ) then
          if ( j == ncol ) then
            continue
          end if
          jk = ij + nrow
        end if

        if ( k == 1 .or. k == 3 ) then
          k1 = 1
        else
          k1 = 2
        end if

        k2 = 5 - k
        xjac = ( x(ij) - ux(k) ) &
          * ( dx(1,k) * ( 2.0D+00 * ux(k) - x(ij) - uc(k1) ) &
          + dx(2,k) * ( ux(k) - x(ij) ) + dx(2,k2) * ( ux(k2) - x(ij) ) )

        xjac = xjac / h / h - 0.5D+00 * ( cx(1,k) + cx(2,k) )

        jac(ij,jk) = xjac

      end do

      if ( i /= 1 .and. j /= 1 ) then
        jk = ij - nrow - 1
        xjac =   ( x(ij) - ux(1) ) * dx(1,1) * ( uc(1) - ux(1) ) &
               + ( x(ij) - ux(3) ) * dx(1,3) * ( uc(1) - ux(3) )
        jac(ij,jk) = xjac / h / h
      end if

      if ( i /= nrow .and. j /= ncol ) then
        jk = ij + nrow + 1
        xjac =   ( x(ij) - ux(2) ) * dx(1,2) * ( uc(2) - ux(2) ) &
               + ( x(ij) - ux(4) ) * dx(1,4) * ( uc(2) - ux(4) )
        jac(ij,jk) = xjac / h / h
      end if

    end do
  end do

  do i = 1, nvar - 1
    jac(i,nvar) = - h * h
  end do

  return
end
subroutine p10_nvar ( option, nvar )

!*****************************************************************************80
!
!! P10_NVAR sets the number of variables for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nvar

  nrow = 6
  ncol = 6
  nvar = nrow * ncol + 1

  return
end
subroutine p10_option_num ( option_num )

!*****************************************************************************80
!
!! P10_OPTION_NUM returns the number of options for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 2

  return
end
subroutine p10_start ( option, nvar, x )

!*****************************************************************************80
!
!! P10_START returns a starting point for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1:nvar) = 0.0D+00

  return
end
subroutine p10_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P10_STEPSIZE returns step sizes for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =     2.000D+00
  hmin =  0.001D+00
  hmax = 10.000D+00

  return
end
subroutine p10_title ( option, title )

!*****************************************************************************80
!
!! P10_TITLE sets the title for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  if ( option == 1 ) then
    title = 'Torsion of a square rod, finite difference, PHI(S)=EXP(5*S).'
  else if ( option == 2 ) then
    title = 'Torsion of a square rod, finite difference, PHI(S)=two levels.'
  else
    title = '???'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P10_TITLE - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized option number.'
    stop
  end if

  return
end
subroutine p11_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P11_FUN evaluates the function for problem 11.
!
!  Title:
!
!    Torsion of a square rod, finite element solution.
!
!  Description:
!
!    The problem is a boundary value problem on (0,1) x (0,1)
!    of the form:
!
!      - d/dx ( PHI ( dU/dx, dU/dy ) * dU/dx )
!      - d/dy ( PHI ( dU/dx, dU/dy ) * dU/dy ) = G ( U, LAMBDA )
!
!    On the 2-dimensional region [0,1] x [0,1], a regular square grid
!    is used.  If there are NSIDE nodes on a side, then the spacing
!    is H=1/(NSIDE-1).  The nodes are ordered from left to right, and
!    from bottom to top, as are the resulting square elements:
!
!      21---22---23---24---25
!       !    !    !    !    !
!       ! 13 ! 14 ! 15 ! 16 !
!       !    !    !    !    !
!      16---17---18---19---20
!       !    !    !    !    !
!       ! 09 ! 10 ! 11 ! 12 !
!       !    !    !    !    !
!      11---12---13---14---15
!       !    !    !    !    !
!       ! 05 ! 06 ! 07 ! 08 !
!       !    !    !    !    !
!      06---07---08---09---10
!       !    !    !    !    !
!       ! 01 ! 02 ! 03 ! 04 !
!       !    !    !    !    !
!      01---02---03---04---05
!
!    On a single element, the local ordering of nodes and shape
!    functions is
!
!      3----4
!      !    !
!      !    !
!      !    !
!      1----2
!
!    Linear elements are used.  If H is the length of a side, the shape
!    function in a particular element associated with node 1 is:
!
!      PSI(X,Y) = ( X - XRIGHT ) * ( Y - YTOP ) / H**2
!
!    where
!
!      XRIGHT is the X coordinate of the right hand side of the element,
!      YTOP is the Y coordinate of the top side of the element.
!
!  Options:
!
!    OPTION = 1:
!
!      PHI = exp ( 5 * ( dUdX**2 + dUdY**2 ) )
!
!      G ( U, LAMBDA ) = - 5 * LAMBDA
!
!    OPTION = 2:
!
!      Let S = ( dUdX**2 + dUdY**2 ),
!      SBAR = ( 40 * S - 13 ) / 7
!
!      if ( S <= 0.15 ) then
!        PHI = 1.0D+00
!      else if ( 0.15 <= S <= 0.50 ) then
!        PHI = 5.5 + 4.5 * ( 3 * SBAR - SBAR**3 )
!      else if ( 0.50 <= S ) then
!        PHI = 10.0D+00
!
!      G ( U, LAMBDA ) = - 10 * LAMBDA
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Werner Rheinboldt,
!    On the Solution of Some Nonlinear Equations Arising in the
!    Application of Finite Element Methods,
!    in The Mathematics of Finite Elements and Applications II,
!    edited by John Whiteman,
!    Academic Press, London, 1976, pages 465-482,
!    LC: TA347.F5.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ), parameter :: ngauss = 4
  integer ( kind = 4 ), parameter :: nshape = 4
  integer ( kind = 4 ), parameter :: nside = 5

  real ( kind = 8 ) dpsidx(nshape)
  real ( kind = 8 ) dpsidy(nshape)
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) fx(nvar-1)
  real ( kind = 8 ) flam
  real ( kind = 8 ) flampl
  real ( kind = 8 ) flampu
  real ( kind = 8 ) hside
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) option
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) irowm1
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) jgauss
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) jrowm1
  integer ( kind = 4 ) kshape
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) nod
  integer ( kind = 4 ) node(nshape)
  real ( kind = 8 ) phi
  real ( kind = 8 ) phip
  real ( kind = 8 ) psi(nshape)
  real ( kind = 8 ) uval
  real ( kind = 8 ) wgauss(ngauss)
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xgauss(ngauss)
  real ( kind = 8 ) xmid
  real ( kind = 8 ) xval
  real ( kind = 8 ) ygauss(ngauss)
  real ( kind = 8 ) ymid
  real ( kind = 8 ) yval

  fx(1:nvar-1) = 0.0D+00

  nelem = ( nside - 1 ) * ( nside - 1 )
  hside = 1.0D+00 / real ( nside - 1, kind = 8 )

  do ielem = 1, nelem
!
!  From the element number IELEM, compute the indices of the four
!  corners, in the order SW, SE, NW, NE.
!
    irowm1 = ( ielem - 1 ) / ( nside - 1 )
    icol = ielem - irowm1 * ( nside - 1 )
    irow = irowm1 + 1

    xmid = hside * real ( 2 * icol - 1 ) / 2.0D+00
    ymid = hside * real ( 2 * irow - 1 ) / 2.0D+00

    node(1) = irowm1 * nside + icol
    node(2) = node(1) + 1
    node(3) = node(1) + nside
    node(4) = node(3) + 1
!
!  Get the Gauss points for this element.
!
    call p11_gauss ( hside, xmid, ymid, wgauss, xgauss, ygauss )
!
!  For each Gauss point in this element, evaluate the integrand.
!
    do jgauss = 1, ngauss

      xval = xgauss(jgauss)
      yval = ygauss(jgauss)
!
!  Evaluate the shape functions PSI and their derivatives.
!
      call p11_shape ( hside, xmid, xval, ymid, yval, psi, dpsidx, dpsidy )
!
!  Evaluate U and its derivatives.
!
      uval = 0.0D+00
      dudx = 0.0D+00
      dudy = 0.0D+00
      do i = 1, nshape
        uval = uval + x(node(i)) * psi(i)
        dudx = dudx + x(node(i)) * dpsidx(i)
        dudy = dudy + x(node(i)) * dpsidy(i)
      end do
!
!  Evaluate PHI ( DUDX, DUDY ).
!
      call p11_phi ( dudx, dudy, option, phi, phip )
!
!  Evaluate G ( U, LAMBDA ).
!
      call p11_gul ( option, x(nvar), flam, flampl, flampu )
!
!  Compute the inner product of the equation with each shape function
!  and add to the appropriate function.
!
      do kshape = 1, nshape

        nod = node(kshape)
        jrowm1 = ( nod - 1 ) / nside
        jcol = nod - jrowm1 * nside
        jrow = jrowm1 + 1

        if ( jrow == 1 .or. jrow == nside .or. &
             jcol == 1 .or. jcol == nside ) then

          fx(nod) = x(nod)

        else

          fx(nod) = fx(nod) + wgauss(jgauss) * hside * hside * &
                  ( phi * ( dudx * dpsidx(kshape) + dudy * dpsidy(kshape) ) &
                  + flam * psi(kshape) )

        end if

      end do

    end do

  end do

  return
end
subroutine p11_gauss ( hside, xmid, ymid, wgauss, xgauss, ygauss )

!*****************************************************************************80
!
!! P11_GAUSS returns the Gauss quadrature abscissas and weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) HSIDE, the length of a side of the square.
!
!    Input, real ( kind = 8 ) XMID, YMID, the coordinates of the center
!    of the square.
!
!    Output, real ( kind = 8 ) WGAUSS(4), the weights of the Gauss points.
!    The weights are normalized for a square of unit area.
!
!    Output, real ( kind = 8 ) XGAUSS(4), YGAUSS(4), the coordinates of the
!    Gauss points.
!
  implicit none

  real ( kind = 8 ) alfa
  real ( kind = 8 ) hside
  real ( kind = 8 ) wgauss(4)
  real ( kind = 8 ) xgauss(4)
  real ( kind = 8 ) xmid
  real ( kind = 8 ) ygauss(4)
  real ( kind = 8 ) ymid

  alfa = 1.0D+00 / ( 2.0D+00 * sqrt ( 3.0D+00 ) )

  wgauss(1:4) = 0.25D+00

  xgauss(1) = xmid - alfa * hside
  xgauss(2) = xmid + alfa * hside
  xgauss(3) = xmid - alfa * hside
  xgauss(4) = xmid + alfa * hside

  ygauss(1) = ymid - alfa * hside
  ygauss(2) = ymid - alfa * hside
  ygauss(3) = ymid + alfa * hside
  ygauss(4) = ymid + alfa * hside

  return
end
subroutine p11_gul ( option, lambda, flam, flampl, flampu )

!*****************************************************************************80
!
!! P11_GUL computes G(U,LAMBDA) and dG/dU and dG/dLAMBDA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, real ( kind = 8 ) LAMBDA, the value of LAMBDA.
!
!    Output, real ( kind = 8 ) FLAM, FLAMPL, FLAMPU, the values of F(U,LAMBDA),
!    d F(U,LAMBDA)/d LAMBDA, and d F(U,LAMBDA)/d U.
!
  implicit none

  real ( kind = 8 ) flam
  real ( kind = 8 ) flampl
  real ( kind = 8 ) flampu
  integer ( kind = 4 ) option
  real ( kind = 8 ) lambda

  if ( option == 1 ) then
    flam = - 5.0D+00 * lambda
    flampl = - 5.0D+00
    flampu = 0.0D+00
  else if ( option == 2 ) then
    flam = - 10.0D+00 * lambda
    flampl = - 10.0D+00
    flampu = 0.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P11_GUL - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal option = ', option
    stop
  end if

  return
end
subroutine p11_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P11_JAC evaluates the jacobian for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ), parameter :: ngauss = 4
  integer ( kind = 4 ), parameter :: nshape = 4
  integer ( kind = 4 ), parameter :: nside = 5

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) dpsidx(nshape)
  real ( kind = 8 ) dpsidy(nshape)
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) flam
  real ( kind = 8 ) flampl
  real ( kind = 8 ) flampu
  real ( kind = 8 ) hside
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) option
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) irowm1
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) jgauss
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) jrowm1
  integer ( kind = 4 ) kshape
  integer ( kind = 4 ) lshape
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) nod
  integer ( kind = 4 ) nod2
  integer ( kind = 4 ) node(nshape)
  real ( kind = 8 ) phi
  real ( kind = 8 ) phip
  real ( kind = 8 ) psi(nshape)
  real ( kind = 8 ) term1
  real ( kind = 8 ) term2
  real ( kind = 8 ) uval
  real ( kind = 8 ) wgauss(ngauss)
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xgauss(ngauss)
  real ( kind = 8 ) xmid
  real ( kind = 8 ) xval
  real ( kind = 8 ) ygauss(ngauss)
  real ( kind = 8 ) ymid
  real ( kind = 8 ) yval

  jac(1:nvar,1:nvar) = 0.0D+00

  nelem = ( nside - 1 ) * ( nside - 1 )
  hside = 1.0D+00 / real ( nside - 1, kind = 8 )

  do ielem = 1, nelem
!
!  From element number, compute 4 node numbers
!  in the order sw, se, nw, ne.
!
    irowm1 = ( ielem - 1 ) / ( nside - 1 )
    icol = ielem - irowm1 * ( nside - 1 )
    irow = irowm1 + 1
    xmid = hside * real ( 2 * icol - 1 ) / 2.0D+00
    ymid = hside * real ( 2 * irow - 1 ) / 2.0D+00

    node(1) = irowm1 * nside + icol
    node(2) = node(1) + 1
    node(3) = node(1) + nside
    node(4) = node(3) + 1
!
!  Get the Gauss quadrature points for this element.
!
    call p11_gauss ( hside, xmid, ymid, wgauss, xgauss, ygauss )
!
!  At each Gauss point in this element, evaluate the integrand.
!
    do jgauss = 1, ngauss

      xval = xgauss(jgauss)
      yval = ygauss(jgauss)
!
!  Evaluate the shape functions.
!
      call p11_shape ( hside, xmid, xval, ymid, yval, psi, dpsidx, dpsidy )
!
!  Evaluate U and its derivatives.
!
      uval = 0.0D+00
      dudx = 0.0D+00
      dudy = 0.0D+00
      do i = 1, nshape
        uval = uval + x(node(i)) * psi(i)
        dudx = dudx + x(node(i)) * dpsidx(i)
        dudy = dudy + x(node(i)) * dpsidy(i)
      end do
!
!  Evaluate PHI ( DUDX, DUDY ).
!
      call p11_phi ( dudx, dudy, option, phi, phip )
!
!  Evaluate G ( U, LAMBDA ).
!
      call p11_gul ( option, x(nvar), flam, flampl, flampu )
!
!  Compute inner product of equation with each shape function
!  and add to appropriate function.
!
      do kshape = 1, nshape

        nod = node(kshape)
        jrowm1 = ( nod - 1 ) / nside
        jcol = nod - jrowm1 * nside
        jrow = jrowm1 + 1

        if ( jrow == 1 .or. jrow == nside .or. &
             jcol == 1 .or. jcol == nside ) then

          jac(nod,nod) = 1.0D+00

        else

          do lshape = 1, nshape

            nod2 = node(lshape)

            term1 = phi * dpsidx(lshape) + 2.0D+00 * phip * &
              ( dudx * dudx * dpsidx(lshape) + dudx * dudy * dpsidy(lshape) )

            term2 = phi * dpsidy(lshape) + 2.0D+00 * phip * &
              ( dudy * dudx * dpsidx(lshape) + dudy * dudy * dpsidy(lshape) )

            jac(nod,nod2) = jac(nod,nod2) + wgauss(jgauss) * hside * hside * &
              ( term1 * dpsidx(kshape) + term2 * dpsidy(kshape) &
              + flampu * psi(lshape) * psi(kshape) )

          end do

          jac(nod,nvar) = jac(nod,nvar) + wgauss(jgauss) * hside * hside * &
            flampl * psi(kshape)

        end if

      end do

    end do

  end do

  return
end
subroutine p11_nvar ( option, nvar )

!*****************************************************************************80
!
!! P11_NVAR sets the number of variables for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 26

  return
end
subroutine p11_option_num ( option_num )

!*****************************************************************************80
!
!! P11_OPTION_NUM returns the number of options for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 2

  return
end
subroutine p11_phi ( dudx, dudy, option, phi, phip )

!*****************************************************************************80
!
!! P11_PHI is used by problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DUDX, DUDY, the values of dU/dX and dU/dY.
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) PHI, PHIP, the values of PHI(S) and d PHI(S)/d S,
!    where S = dU/dX**2 + dU/dY**2.
!
  implicit none

  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  integer ( kind = 4 ) option
  real ( kind = 8 ) phi
  real ( kind = 8 ) phip
  real ( kind = 8 ) s
  real ( kind = 8 ) sbar

  s = dudx * dudx + dudy * dudy

  if ( option == 1 ) then

    phi = exp ( 5.0D+00 * s )
    phip = 5.0D+00 * exp ( 5.0D+00 * s )

  else if ( option == 2 ) then

    sbar = ( 40.0D+00 * s - 13.0D+00 ) / 7.0D+00

    if ( s <= 0.15D+00 ) then
      phi = 1.0D+00
      phip = 0.0D+00
    else if ( 0.15D+00 <= s .and. s <= 0.5D+00 ) then
      phi = 5.5D+00 + 2.25D+00 * sbar * ( 3.0D+00 - sbar * sbar )
      phip = 2.25D+00 * ( 3.0D+00 - 3.0D+00 * sbar * sbar ) * 40.0D+00 / 7.0D+00
    else
      phi = 10.0D+00
      phip = 0.0D+00
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P11_PHI - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop

  end if

  return
end
subroutine p11_shape ( hside, xmid, xval, ymid, yval, psi, dpsidx, dpsidy )

!*****************************************************************************80
!
!! P11_SHAPE evaluates the shape functions for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) HSIDE, the length of a side of the square.
!
!    Input, real ( kind = 8 ) XMID, the X coordinate of the center of the square.
!
!    Input, real ( kind = 8 ) XVAL, the X coordinate of the point where the shape
!    functions are to be evaluated.
!
!    Input, real ( kind = 8 ) YMID, the Y coordinate of the center of the square.
!
!    Input, real ( kind = 8 ) YVAL, the Y coordinate of the point where the shape
!    functions are to be evaluated.
!
!    Output, real ( kind = 8 ) PSI(4), the value of PSI (the shape functions) at
!    (XVAL,YVAL).  The shape functions are stored in the order
!    SW, SE, NW, NE.
!
!    Output, real ( kind = 8 ) DPSIDX(4), DPSIDY(4), the values of dPSI/dX
!    and dPSI/dY at (XVAL,YVAL).
!
  implicit none

  real ( kind = 8 ) dpsidx(4)
  real ( kind = 8 ) dpsidy(4)
  real ( kind = 8 ) hside
  real ( kind = 8 ) psi(4)
  real ( kind = 8 ) xleft
  real ( kind = 8 ) xmid
  real ( kind = 8 ) xrite
  real ( kind = 8 ) xval
  real ( kind = 8 ) ybot
  real ( kind = 8 ) ymid
  real ( kind = 8 ) ytop
  real ( kind = 8 ) yval
!
!  Set coordinates.
!
  xleft = xmid - 0.5D+00 * hside
  xrite = xmid + 0.5D+00 * hside
  ybot = ymid - 0.5D+00 * hside
  ytop = ymid + 0.5D+00 * hside
!
!  Evaluate the shape functions.
!
  psi(1) =   ( xval - xrite ) * ( yval - ytop ) / hside / hside
  psi(2) = - ( xval - xleft ) * ( yval - ytop ) / hside / hside
  psi(3) = - ( xval - xrite ) * ( yval - ybot ) / hside / hside
  psi(4) =   ( xval - xleft ) * ( yval - ybot ) / hside / hside
!
!  Evaluate the partial derivatives.
!
  dpsidx(1) =   ( yval - ytop ) / hside / hside
  dpsidx(2) = - ( yval - ytop ) / hside / hside
  dpsidx(3) = - ( yval - ybot ) / hside / hside
  dpsidx(4) =   ( yval - ybot ) / hside / hside

  dpsidy(1) =   ( xval - xrite ) / hside / hside
  dpsidy(2) = - ( xval - xleft ) / hside / hside
  dpsidy(3) = - ( xval - xrite ) / hside / hside
  dpsidy(4) =   ( xval - xleft ) / hside / hside

  return
end
subroutine p11_start ( option, nvar, x )

!*****************************************************************************80
!
!! P11_START returns a starting point for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1:nvar) = 0.0D+00

  return
end
subroutine p11_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P11_STEPSIZE returns step sizes for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =     0.12500D+00
  hmin =  0.03125D+00
  hmax =  4.00000D+00

  return
end
subroutine p11_title ( option, title )

!*****************************************************************************80
!
!! P11_TITLE sets the title for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  title = 'Torsion of a square rod, finite element solution.'

  return
end
subroutine p12_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P12_FUN evaluates the function for problem 12.
!
!  Title:
!
!    Materially nonlinear problem.
!
!  Description:
!
!    The problem is the two point boundary value problem
!
!      U'' + LAMBDA * SIN ( U + U**2 + U**3 ) = 0
!
!    with boundary conditions
!
!      U(0) = 0.0D+00
!      U(1) = 0.0D+00
!
!    U is approximated by piecewise polynomials whose coefficients are
!    the unknowns U(1), ..., U(NVAR-1), and the value of LAMBDA is
!    stored as U(NVAR).
!
!  Options:
!
!    OPTION  Polynomials   Continuity
!      1     linear         1
!      2     cubic          1
!      3     cubic          2
!      4     quintic        1
!      5     quintic        2
!      6     quintic        3
!
!    All options use 8 intervals.
!
!  Comments:
!
!    The current program has zero as solution for all X(nvar).
!    Must find bifurcation branch and jump on to it.
!    Perhaps add X(nvar+1) a perturbation to right hand side.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ivo Babuska, Werner Rheinboldt,
!    Reliable Error Estimations and Mesh Adaptation for the Finite
!    Element Method,
!    in International Conference on Computational Methods
!    in Nonlinear Mechanics,
!    edited by John Oden,
!    Elsevier, 1980,
!    ISBN: 0444853820,
!    LC: QA808.I57.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ), parameter :: nbco = 1
  integer ( kind = 4 ), parameter :: nbcz = 1
  integer ( kind = 4 ), parameter :: nint = 8
  integer ( kind = 4 ), parameter :: maxpolys = 6

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) bcone(8)
  real ( kind = 8 ) bczero(8)
  real ( kind = 8 ) coef
  real ( kind = 8 ) dtdx
  real ( kind = 8 ) dtdxl
  real ( kind = 8 ) dtdxr
  real ( kind = 8 ) fx(nvar)
  real ( kind = 8 ) gcoef(8)
  real ( kind = 8 ) gpoint(8)
  real ( kind = 8 ) h2i
  real ( kind = 8 ) h2il
  real ( kind = 8 ) h2ir
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ieqn
  integer ( kind = 4 ) option
  integer ( kind = 4 ) iskip
  integer ( kind = 4 ) ivar
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lhil
  integer ( kind = 4 ) lhir
  integer ( kind = 4 ) lskip
  integer ( kind = 4 ) ncl
  integer ( kind = 4 ) ncr
  integer ( kind = 4 ) ndsum
  integer ( kind = 4 ) npsum
  integer ( kind = 4 ) nderiv
  integer ( kind = 4 ) npolys
  integer ( kind = 4 ) nvary
  real ( kind = 8 ) p12_theta
  real ( kind = 8 ) phi
  real ( kind = 8 ) pl(maxpolys)
  real ( kind = 8 ) pld(maxpolys)
  real ( kind = 8 ) psi
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) uprym
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xc
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr

  bcone(1) = 0.0D+00
  bczero(1) = 0.0D+00

  fx(1:nvar-1) = 0.0D+00

  if ( option == 1 ) then
    npolys = 2
    nderiv = 1
  else if ( option == 2 ) then
    npolys = 4
    nderiv = 1
  else if ( option == 3 ) then
    npolys = 4
    nderiv = 2
  else if ( option == 4 ) then
    npolys = 6
    nderiv = 1
  else if ( option == 5 ) then
    npolys = 6
    nderiv = 2
  else if ( option == 6 ) then
    npolys = 6
    nderiv = 3
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P12_FUN - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  nvary = nint * npolys
!
!  Get the Gauss quadrature rule.
!
  call p12_gauss8 ( gcoef, gpoint )
!
!  Set up the terms A * Y involving the bivariate form
!
!  For each interval I:
!
  do i = 1, nint

    iskip = ( i - 1 ) * npolys
    xl = real ( i - 1 ) / real ( nint, kind = 8 )
    xr = real ( i ) / real ( nint, kind = 8 )
    dtdx = 2.0D+00 / ( xr - xl )
!
!  For each Gauss point, J, evaluate the integrand.
!
    do j = 1, 8

      t = gpoint(j)
      coef = gcoef(j) * ( xr - xl ) / 2.0D+00
      call p12_legendre_val ( t, dtdx, npolys, pl, pld )

      u = 0.0D+00
      uprym = 0.0D+00
      do k = 1, npolys
        u = u + x(iskip+k) * pl(k)
        uprym = uprym + x(iskip+k) * pld(k)
      end do

      phi = - uprym
      psi = x(nvar) * sin ( u * ( 1.0D+00 + u * ( 1.0D+00 + u ) ) )
      lskip = iskip
!
!  Project onto each test function L.
!
      do l = 1, npolys
        ieqn = lskip + l
        fx(ieqn) = fx(ieqn) + coef * ( psi * pl(l) + phi * pld(l) )
      end do

      lskip = lskip + npolys

    end do

  end do
!
!  2. Add the terms B * Z for the continuity of the test functions.
!
!  For each interval I:
!
  do i = 1, nint

    if ( i == 1 ) then
      ncl = nvary
    else
      ncl = nvary + nbcz + ( i - 2 ) * nderiv
    end if

    ncr = nvary + nbcz + ( i - 1 ) * nderiv
    xl = real ( i - 1 ) / real ( nint, kind = 8 )
    xr = real ( i ) / real ( nint, kind = 8 )
    dtdx = 2.0D+00 / ( xr - xl )
!
!  Count conditions at the left endpoint, LHIL, and at right, LHIR.
!  If we are in the first or last interval, one of
!  these will be boundary conditions.
!
    if ( i == 1 ) then
      lhil = nbcz
    else
      lhil = nderiv
    end if

    if ( i == nint ) then
      lhir = nbco
    else
      lhir = nderiv
    end if
!
!  For each test function PL(K):
!
    do k = 1, npolys

      s = r8_mop ( k + 1 )
      ieqn = ( i - 1 ) * npolys + k
!
!  Apply the boundary conditions.
!
      h2i = 1.0D+00
      do l = 1, lhil
        s = - s
        ivar = ncl + l
        fx(ieqn) = fx(ieqn) + s * x(ivar) * h2i * p12_theta ( l, k )
        h2i = h2i * dtdx
      end do

      h2i = 1.0D+00
      do l = 1, lhir
        ivar = ncr + l
        fx(ieqn) = fx(ieqn) + x(ivar) * h2i * p12_theta ( l, k )
        h2i = h2i * dtdx
      end do

    end do

  end do
!
!  3. Create the C * Y terms for U and its derivatives.
!  One equation is generated for component and condition.
!
  npsum = 0
  dtdxr = 0.0D+00
  dtdxl = 0.0D+00
!
!  For each node:
!
  ndsum = nvary

  do i = 1, nint + 1

    if ( 1 < i ) then
      xl = real ( i - 2 ) / real ( nint, kind = 8 )
    end if

    xc = real ( i - 1 ) / real ( nint, kind = 8 )

    if ( i < nint + 1 ) then
      xr = real ( i ) / real ( nint, kind = 8 )
    end if

    if ( xc /= xl ) then
      dtdxl = 2.0D+00 / ( xc - xl )
    end if

    if ( xr /= xc ) then
      dtdxr = 2.0D+00 / ( xr - xc )
    end if

    h2il = 1.0D+00
    h2ir = 1.0D+00
!
!  Count the conditions:
!
    if ( i == 1 ) then
      khi = nbcz
    else if ( i < nint + 1 ) then
      khi = nderiv
    else if ( i == nint + 1 ) then
      khi = nbco
    end if

    do k = 1, khi

      s = r8_mop ( k + 1 )
!
!  Set up the term from the left hand interval.
!
      ieqn = ndsum + k

      if ( i == 1 ) then

        fx(ieqn) = fx(ieqn) + bczero(k)

      else

        do l = 1, npolys
          ivar = npsum + l - npolys
          fx(ieqn) = fx(ieqn) + x(ivar) * h2il * p12_theta ( k, l )
        end do

      end if
!
!  Set up the term from the right hand interval.
!
      if ( i == nint + 1 ) then

        fx(ieqn) = fx(ieqn) - bcone(k)

      else

        do l = 1, npolys
          ivar = npsum + l
          s = - s
          fx(ieqn) = fx(ieqn) + s * x(ivar) * h2ir * p12_theta(k,l)
        end do
      end if

      h2il = h2il * dtdxl
      h2ir = h2ir * dtdxr

    end do

    ndsum = ndsum + khi
    npsum = npsum + npolys

  end do

  return
end
subroutine p12_gauss8 ( gcoef, gpoint )

!*****************************************************************************80
!
!! P12_GAUSS8 returns an 8 point Gauss quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) GCOEF(8), the coefficients for the quadrature rule,
!    normalized for the interval [-1,1].
!
!    Output, real ( kind = 8 ) GPOINT(8), the abscissas for the quadrature rule,
!    normalized for the interval [-1,1].
!
  implicit none

  real ( kind = 8 ) gcoef(8)
  real ( kind = 8 ) gpoint(8)

  gcoef(1) = 0.1012285363D+00
  gcoef(2) = 0.2223810345D+00
  gcoef(3) = 0.3137066459D+00
  gcoef(4) = 0.3626837834D+00
  gcoef(5) = 0.3626837834D+00
  gcoef(6) = 0.3137066459D+00
  gcoef(7) = 0.2223810345D+00
  gcoef(8) = 0.1012285363D+00

  gpoint(1) = - 0.9602898565D+00
  gpoint(2) = - 0.7966664774D+00
  gpoint(3) = - 0.5255324099D+00
  gpoint(4) = - 0.1834346425D+00
  gpoint(5) =   0.1834346425D+00
  gpoint(6) =   0.5255324099D+00
  gpoint(7) =   0.7966664774D+00
  gpoint(8) =   0.9602898565D+00

  return
end
subroutine p12_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P12_JAC evaluates the jacobian for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ), parameter :: nbco = 1
  integer ( kind = 4 ), parameter :: nbcz = 1
  integer ( kind = 4 ), parameter :: nint = 8
  integer ( kind = 4 ), parameter :: maxpolys = 6

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) dbcodt(8)
  real ( kind = 8 ) dbczdt(8)
  real ( kind = 8 ) coef
  real ( kind = 8 ) dtdx
  real ( kind = 8 ) gcoef(8)
  real ( kind = 8 ) gpoint(8)
  real ( kind = 8 ) h2i
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ieqn
  integer ( kind = 4 ) option
  integer ( kind = 4 ) iskip
  integer ( kind = 4 ) ivar
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khil
  integer ( kind = 4 ) khir
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lhil
  integer ( kind = 4 ) lhir
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncl
  integer ( kind = 4 ) ncr
  integer ( kind = 4 ) nderiv
  integer ( kind = 4 ) npolys
  integer ( kind = 4 ) npsum
  integer ( kind = 4 ) nvary
  real ( kind = 8 ) p12_theta
  real ( kind = 8 ) phipt
  real ( kind = 8 ) phipu
  real ( kind = 8 ) phipup
  real ( kind = 8 ) pl(maxpolys)
  real ( kind = 8 ) pld(maxpolys)
  real ( kind = 8 ) psipt
  real ( kind = 8 ) psipu
  real ( kind = 8 ) psipup
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) uprym
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr

  jac(1:nvar,1:nvar) = 0.0D+00

  dbcodt(1) = 0.0D+00
  dbczdt(1) = 0.0D+00

  if ( option == 1 ) then
    npolys = 2
    nderiv = 1
  else if ( option == 2 ) then
    npolys = 4
    nderiv = 1
  else if ( option == 3 ) then
    npolys = 4
    nderiv = 2
  else if ( option == 4 ) then
    npolys = 6
    nderiv = 1
  else if ( option == 5 ) then
    npolys = 6
    nderiv = 2
  else if ( option == 6 ) then
    npolys = 6
    nderiv = 3
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P12_JAC - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  nvary = nint * npolys
!
!  Get the Gauss quadrature rule.
!
  call p12_gauss8 ( gcoef, gpoint )
!
!  1.  Set up the terms from the bivariate form A * Y:
!
  do i = 1, nint

    iskip = ( i - 1 ) * npolys
    xl = real ( i - 1 ) / real ( nint, kind = 8 )
    xr = real ( i ) / real ( nint, kind = 8 )
    dtdx = 2.0D+00 / ( xr - xl )
!
!  For each Gauss point in the interval:
!
    do j = 1, 8

      t = gpoint(j)
      coef = gcoef(j) * ( xr - xl ) / 2.0D+00
      call p12_legendre_val ( t, dtdx, npolys, pl, pld )

      u = 0.0D+00
      uprym = 0.0D+00
      do k = 1, npolys
        u = u + x(iskip+k) * pl(k)
        uprym = uprym + x(iskip+k) * pld(k)
      end do

      phipu = 0.0D+00
      phipup = - 1.0D+00
      phipt = 0.0D+00

      psipu = x(nvar) * ( 1.0D+00 + u * ( 2.0D+00 + 3.0D+00 * u ) ) &
        * cos ( u * ( 1.0D+00 + u * ( 1.0D+00 + u ) ) )
      psipup = 0.0D+00
      psipt = sin ( u * ( 1.0D+00 + u * ( 1.0D+00 + u ) ) )
!
!  For each Legendre polynomial coefficient:
!
      do l = 1, npolys

        ieqn = iskip + l
        jac(ieqn,nvar) = jac(ieqn,nvar) + coef * ( psipt * pl(l) + &
          phipt * pld(l) )
!
!  For each Y-coefficient of U:
!
        do n = 1, npolys
          ivar = npolys * ( i - 1 ) + n
          jac(ieqn,ivar) = jac(ieqn,ivar) + coef * ( &
                ( psipu * pl(n) + psipup * pld(n) ) * pl(l) &
              + ( phipu * pl(n) + phipup * pld(n) ) * pld(l) )

        end do

      end do

    end do

  end do
!
!  2. Add the terms involving the continuity of the test functions
!  which are the terms B * Z in F = A * Y + B * Z.
!
  do i = 1, nint

    if ( i == 1 ) then
      ncl = nvary
    else
      ncl = nvary + nbcz + ( i - 2 ) * nderiv
    end if

    ncr = nvary + nbcz + ( i - 1 ) * nderiv
    xl = real ( i - 1 ) / real ( nint, kind = 8 )
    xr = real ( i ) / real ( nint, kind = 8 )
    dtdx = 2.0D+00 / ( xr - xl )
!
!  For the polynomials used in approximating each U,
!  count conditions at left endpoint, LHIL, and at right, LHIR.
!
      if ( i == 1 ) then
        lhil = nbcz
      else
        lhil = nderiv
      end if

      if ( i == nint ) then
        lhir = nbco
      else
        lhir = nderiv
      end if
!
!  For each test function PL(K).
!
      do k = 1, npolys
        s = r8_mop ( k + 1 )
        ieqn = ( i - 1 ) * npolys + k
!
!  Consider the conditions:
!
        h2i = 1.0D+00
        do l = 1, lhil
          s = - s
          ivar = ncl + l
          jac(ieqn,ivar) = s * h2i * p12_theta ( l, k )
          h2i = h2i * dtdx
        end do
!
!  Evaluate contribution from right endpoint.
!
        h2i = 1.0D+00
        do l = 1, lhir
          ivar = ncr + l
          jac(ieqn,ivar) = h2i * p12_theta ( l, k )
          h2i = h2i * dtdx
        end do
      end do

    end do
!
!  3. Create the terms for the U functions and their derivatives
!  the matrix terms C * Y.
!
  do i = 1, nint

    if ( i == 1 ) then
      ncl = nvary
    else
      ncl = nvary + nbcz + ( i - 2 ) * nderiv
    end if

    ncr = nvary + nbcz + ( i - 1 ) * nderiv
    npsum = ( i - 1 ) * npolys
    xl = real ( i - 1 ) / real ( nint, kind = 8 )
    xr = real ( i ) / real ( nint, kind = 8 )
    dtdx = 2.0D+00 / ( xr - xl )

    h2i = 1.0D+00
!
!  Count the conditions:
!
    if ( i == 1 ) then
      khil = nbcz
    else
      khil = nderiv
    end if
!
!  Left hand term:
!
    do k = 1, khil

     ieqn = ncl + k

      if ( i == 1 ) then
        jac(ieqn,nvar) = dbczdt(k)
      end if

      s = r8_mop ( k + 1 )
      do l = 1, npolys
        ivar = npsum + l
        s = - s
        jac(ieqn,ivar) = s * h2i * p12_theta ( k, l )
      end do
      h2i = h2i * dtdx
    end do

    ncl = ncl + khil
!
!  Right hand term:
!
    h2i = 1.0D+00

    if ( i == nint ) then
      khir = nbco
    else
      khir = nderiv
    end if

    do k = 1, khir

      ieqn = ncr + k

      if ( i == nint ) then
        jac(ieqn,nvar) = - dbcodt(k)
      else
        jac(ieqn,nvar) = 0.0D+00
      end if

      do l = 1, npolys
        ivar = npsum + l
        jac(ieqn,ivar) = h2i * p12_theta ( k, l )
      end do

      h2i = h2i * dtdx

    end do

    ncr = ncr + khir
    npsum = npsum + npolys

  end do

  return
end
subroutine p12_legendre_val ( t, dtdx, npolys, pl, pld )

!*****************************************************************************80
!
!! P12_LEGENDRE_VAL evaluates the Legendre polynomials and derivatives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the argument of the Legendre polynomials, in
!    the normalized interval [-1,1].
!
!    Input, real ( kind = 8 ) DTDX, the value of the quantity dTdX at the point
!    X.  In the most common case, this is simply the relationship
!    between the width of the normalized T interval (2), and the
!    width of the X interval to which the Legendre polynomial
!    arguments have been mapped.  DTDX is needed so that the
!    computed values PLD can be converted from dPL/dT to dPL/dX.
!
!    Input, integer ( kind = 4 ) NPOLYS, the number of Legendre polynomials to
!    evaluate.  If NPOLYS is 1, then only the constant polynomial
!    is evaluated,  NPOLYS = 2 means the constant and linear, and so on.
!
!    Output, real ( kind = 8 ) PL(NPOLYS), PLD(NPOLYS), the values of PL(X)
!    and dPL(X)/dX at the point X which has normalized coordinate T.
!
  implicit none

  integer ( kind = 4 ) npolys

  real ( kind = 8 ) a
  real ( kind = 8 ) dtdx
  integer ( kind = 4 ) i
  real ( kind = 8 ) pl(npolys)
  real ( kind = 8 ) pld(npolys)
  real ( kind = 8 ) t

  if ( 1 <= npolys ) then
    pl(1) = 1.0D+00
    pld(1) = 0.0D+00
  end if

  if ( 2 <= npolys ) then
    pl(2) = t
    pld(2) = 1.0D+00
  end if

  a = 0.0D+00
  do i = 3, npolys

    a = a + 1.0D+00

    pl(i) = ( ( 2.0D+00 * a + 1.0D+00 ) * t * pl(i-1) - a * pl(i-2) ) &
      / ( a + 1.0D+00 )

    pld(i) = ( ( 2.0D+00 * a + 1.0D+00 ) * ( t * pld(i-1) + pl(i-1) ) &
           - a * pld(i-2) ) / ( a + 1.0D+00 )

  end do

  pld(1:npolys) = dtdx * pld(1:npolys)

  return
end
subroutine p12_nvar ( option, nvar )

!*****************************************************************************80
!
!! P12_NVAR sets the number of variables for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nbco
  integer ( kind = 4 ) nbcz
  integer ( kind = 4 ) nderiv
  integer ( kind = 4 ) nint
  integer ( kind = 4 ) npolys
  integer ( kind = 4 ) nvar
  integer ( kind = 4 ) nvary
  integer ( kind = 4 ) nvarz

  if ( option == 1 ) then
    npolys = 2
    nderiv = 1
  else if ( option == 2 ) then
    npolys = 4
    nderiv = 1
  else if ( option == 3 ) then
    npolys = 4
    nderiv = 2
  else if ( option == 4 ) then
    npolys = 6
    nderiv = 1
  else if ( option == 5 ) then
    npolys = 6
    nderiv = 2
  else if ( option == 6 ) then
    npolys = 6
    nderiv = 3
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P12_NVAR - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  nint = 8
  nvary = nint * npolys
  nbcz = 1
  nbco = 1
  nvarz = nbcz + ( nint - 1 ) * nderiv + nbco
  nvar = nvary + nvarz + 1

  return
end
subroutine p12_option_num ( option_num )

!*****************************************************************************80
!
!! P12_OPTION_NUM returns the number of options for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 6

  return
end
subroutine p12_start ( option, nvar, x )

!*****************************************************************************80
!
!! P12_START returns a starting point for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1:nvar) = 0.0D+00

  return
end
subroutine p12_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P12_STEPSIZE returns step sizes for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =     2.000D+00
  hmin =  0.001D+00
  hmax = 10.000D+00

  return
end
function p12_theta ( i, j )

!*****************************************************************************80
!
!! P12_THETA is a utility routine used in problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the indices of THETA.
!
!    Output, real ( kind = 8 ) P12_THETA, the value of THETA(I,J).
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p12_theta
  real ( kind = 8 ), save, dimension(nmax,nmax) :: theta

  data theta / &
             1.0,       0.0,       0.0,       0.0,       0.0, &
             0.0,       0.0,       0.0,       0.0,       0.0, &
             1.0,       1.0,       0.0,       0.0,       0.0, &
             0.0,       0.0,       0.0,       0.0,       0.0, &
             1.0,       3.0,       3.0,       0.0,       0.0, &
             0.0,       0.0,       0.0,       0.0,       0.0, &
             1.0,       6.0,      15.0,      15.0,       0.0, &
             0.0,       0.0,       0.0,       0.0,       0.0, &
             1.0,      10.0,      45.0,     105.0,     105.0, &
             0.0,       0.0,       0.0,       0.0,       0.0, &
             1.0,      15.0,     105.0,     420.0,     945.0, &
           945.0,       0.0,       0.0,       0.0,       0.0, &
             1.0,      21.0,     210.0,    1260.0,    4725.0, &
         10395.0,   10395.0,       0.0,       0.0,       0.0, &
             1.0,      28.0,     378.0,    3150.0,   17325.0, &
         62370.0,  135135.0,  135135.0,       0.0,       0.0, &
             1.0,      36.0,     630.0,    6930.0,   51975.0, &
        270270.0,  945945.0, 2027025.0, 2027025.0,       0.0, &
             1.0,      45.0,     990.0,   13860.0,  135135.0, &
       945945.0,  4729725.0, 16216200.0,34459425.0,34459425.0 /

  if ( i < 1 .or. nmax < i .or. j < 1 .or. nmax < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P12_THETA - Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    stop
  end if

  p12_theta = theta ( i, j )

  return
end
subroutine p12_title ( option, title )

!*****************************************************************************80
!
!! P12_TITLE sets the title for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  if ( option == 1 ) then
    title = 'Materially nonlinear problem, NPOLYS = 2, NDERIV = 1.'
  else if ( option == 2 ) then
    title = 'Materially nonlinear problem, NPOLYS = 4, NDERIV = 1.'
  else if ( option == 3 ) then
    title = 'Materially nonlinear problem, NPOLYS = 4, NDERIV = 2.'
  else if ( option == 4 ) then
    title = 'Materially nonlinear problem, NPOLYS = 6, NDERIV = 1.'
  else if ( option == 5 ) then
    title = 'Materially nonlinear problem, NPOLYS = 6, NDERIV = 2.'
  else if ( option == 6 ) then
    title = 'Materially nonlinear problem, NPOLYS = 6, NDERIV = 3.'
  else
    title = '???'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P12_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
subroutine p13_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P13_FUN evaluates the function for problem 13.
!
!  Discussion:
!
!    Simpson's mildly nonlinear boundary value problem.
!
!    The continuous problem is defined on the unit square,
!    and has the form:
!
!      - Laplacian ( U(X,Y) ) = LAMBDA * F ( U(X,Y) )
!
!    for points within the unit square, and boundary condition
!
!      U(X,Y) = 0.
!
!    The continuous problem is discretized with a uniform M by M
!    mesh of point in the interior.  Let DEL9 be the nine point
!    discrete Laplacian operator, and DEL5 the five point discrete
!    Laplacian operator.  Then the discrete problem has the form:
!
!      DEL9 U + lambda * ( F(U) + H**2 * DEL5 ( F(U) ) / 12 )  =  0.0D+00
!
!    where H is the mesh spacing.
!
!    The options allow a choice of M and the right hand side function F.
!
!    OPTION   M    NVAR  F(U)
!
!       1   8     65   exp ( U )
!       2   8     65   ( 100 + 100 * U + 51 * U**2 ) / ( 100 + U**2 )
!       3  12    145   exp ( U )
!       4  12    145   ( 100 + 100 * U + 51 * U**2 ) / ( 100 + U**2 )
!       5  16    257   exp ( U )
!       6  16    257   ( 100 + 100 * U + 51 * U**2 ) / ( 100 + U**2 )
!
!    Melhem lists a limit point in LAMBDA for each of the option cases
!    above.  Letting U* be the value of the computed solution at the
!    center point of the grid, we have:
!
!    OPTION   Lambda          U*
!
!     1     6.807504     1.391598
!     2     7.980356     2.272364
!     3     6.808004     1.391657
!     4     7.981426     2.273045
!     5     6.808087     1.391656
!     6     7.981605     2.273159
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Rami Melhem, Werner Rheinboldt,
!    A Comparison of Methods for Determining Turning Points of Nonlinear Equations,
!    Computing,
!    Volume 29, Number 3, September 1982, pages 201-226.
!
!    Bruce Simpson,
!    A Method for the Numerical Determination of Bifurcation
!    States of Nonlinear Systems of Equations,
!    SIAM Journal on Numerical Analysis,
!    Volume 12, Number 3, June 1975, pages 439-451.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar-1)
  integer ( kind = 4 ) option
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(nvar)
!
!  Compute M, the order of the square grid, such that M*M = NVAR-1.
!
  m = nint ( sqrt ( real ( nvar - 1, kind = 8 ) ) )

  lambda = x(nvar)
  call p13_fx2 ( option, m, x, lambda, fx )

  return
end
subroutine p13_fx2 ( option, m, u, lambda, fx )

!*****************************************************************************80
!
!! P13_FX2 computes the function by recasting it on a square grid.
!
!  Discussion:
!
!    For M = 4, there are M*M = 16 U variables plus LAMBDA.
!
!    The ordering of the U variables is suggested by the diagram, in which
!    "0" indicates a point where U is zero, and a nonzero value indicates
!    the index in the vector U of the corresponding value:
!
!      |
!     1.0      0   0   0   0   0   0
!      |       0  13  14  15  16   0
!      |       0   9  10  11  12   0
!      Y       0   5   6   7   8   0
!      |       0   1   2   3   4   0
!     0.0      0   0   0   0   0   0
!      |
!      +--0.0------X----------1.0--->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) M, the number of grid points along a side of
!    the square.
!
!    Input, real ( kind = 8 ) U(M,M), the value of the grid function at the
!    grid points.
!
!    Input, real ( kind = 8 ) LAMBDA, the value of the parameter.
!
!    Output, real ( kind = 8 ) FX(M,M), the value of the function.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) del5f
  real ( kind = 8 ) del9u
  real ( kind = 8 ) fc
  real ( kind = 8 ) fe
  real ( kind = 8 ) fn
  real ( kind = 8 ) fs
  real ( kind = 8 ) fw
  real ( kind = 8 ) fx(m,m)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) lambda
  real ( kind = 8 ) p13_gx
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) uc
  real ( kind = 8 ) ue
  real ( kind = 8 ) un
  real ( kind = 8 ) une
  real ( kind = 8 ) unw
  real ( kind = 8 ) us
  real ( kind = 8 ) use
  real ( kind = 8 ) usw
  real ( kind = 8 ) uw

  h = 1.0D+00 / real ( m + 1, kind = 8 )

  do i = 1, m
    do j = 1, m
!
!  Evaluate the solution on the grid:
!
!   UNW-UN--UNE
!    |   |   |
!   UW--UC--UE
!    |   |   |
!   USW-US--USE
!
      uc = u(i,j)

      if ( i < m ) then
        un = u(i+1,j)
      else
        un = 0.0D+00
      end if

      if ( 1 < i ) then
        us = u(i-1,j)
      else
        us = 0.0D+00
      end if

      if ( j < m ) then
        ue = u(i,j+1)
      else
        ue = 0.0D+00
      end if

      if ( 1 < j ) then
        uw = u(i,j-1)
      else
        uw = 0.0D+00
      end if

      if ( 1 < i .and. 1 < j ) then
        usw = u(i-1,j-1)
      else
        usw = 0.0D+00
      end if

      if ( 1 < i .and. j < m ) then
        use = u(i-1,j+1)
      else
        use = 0.0D+00
      end if

      if ( i < m .and. 1 < j ) then
        unw = u(i+1,j-1)
      else
        unw = 0.0D+00
      end if

      if ( i < m .and. j < m ) then
        une = u(i+1,j+1)
      else
        une = 0.0D+00
      end if
!
!  Evaluate the right hand side on the grid.
!
!      FN
!      |
!   FW-FC-FE
!      |
!      FS
!
      fc = p13_gx ( option, uc )
      fn = p13_gx ( option, un )
      fs = p13_gx ( option, us )
      fe = p13_gx ( option, ue )
      fw = p13_gx ( option, uw )
!
!  Compute the 9 point approximation to Laplacian U.
!
      del9u = ( - 20.0D+00 * uc + 4.0D+00 * ( un + us + ue + uw ) &
              + une + unw + use + usw ) / ( 6.0D+00 * h * h )

      del5f = fc + h * h * ( - 4.0D+00 * fc + fn + fs + fe + fw ) / 12.0D+00

      fx(i,j) = del9u + lambda * del5f

    end do
  end do

  return
end
function p13_gp ( option, u )

!*****************************************************************************80
!
!! P13_GP evaluates the derivative of the right hand side function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, real ( kind = 8 ) U, the value of the argument.
!
!    Output, real ( kind = 8 ) P13_GP, the derivative of the right hand side function
!    at U.
!
  implicit none

  integer ( kind = 4 ) option
  real ( kind = 8 ) p13_gp
  real ( kind = 8 ) u

  if ( option == 1 .or. option == 3 .or. option == 5 ) then
    p13_gp = exp ( u )
  else if ( option == 2 .or. option == 4 .or. option == 6 ) then
    p13_gp = ( 1.0D+00 + u - 0.01D+00 * u * u ) &
      / ( 1.0D+00 + 0.01D+00 * u * u )**2
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P13_GP - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
function p13_gx ( option, u )

!*****************************************************************************80
!
!! P13_GX evaluates the right hand side function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, real ( kind = 8 ) U, the value of the argument.
!
!    Output, real ( kind = 8 ) P13_GX, the right hand side function at U.
!
  implicit none

  integer ( kind = 4 ) option
  real ( kind = 8 ) p13_gx
  real ( kind = 8 ) u

  if ( option == 1 .or. option == 3 .or. option == 5 ) then
    p13_gx = exp ( u )
  else if ( option == 2 .or. option == 4 .or. option == 6 ) then
    p13_gx = ( 100.0D+00 + 100.0D+00 * u + 51.0D+00 * u * u ) &
      / ( 100.0D+00 + u * u )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P13_GX - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
subroutine p13_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P13_JAC evaluates the jacobian for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(nvar)

  jac(1:nvar,1:nvar) = 0.0D+00

  m = nint ( sqrt ( real ( nvar - 1, kind = 8 ) ) )
  lambda = x(nvar)

  call p13_jac2 ( option, m, nvar, lambda, x, jac )

  return
end
subroutine p13_jac2 ( option, m, nvar, lambda, u, jac )

!*****************************************************************************80
!
!! P13_JAC2 computes the jacobian by recasting it on a square grid.
!
!  Discussion:
!
!    Actually, to stave off insanity, we only "recast" the variables into
!    a 2D array that corresponds to the spatial ordering of the grid.
!    We leave the jacobian in its original arrangement, which assumes
!    a linear ordering of variables and equations, and we simply
!    compute the equation and variable indices of the jacobian when
!    we are ready to put entries into it.  This approach seems to produce
!    a smaller amount of cosmic grief than the alternatives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) M, the number of grid points on a side
!    of the square.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) LAMBDA, the value of the parameter.
!
!    Input, real ( kind = 8 ) U(M,M), the value of the grid function.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nvar

  real ( kind = 8 ) del5f
  real ( kind = 8 ) fc
  real ( kind = 8 ) fcp
  real ( kind = 8 ) fe
  real ( kind = 8 ) fep
  real ( kind = 8 ) fn
  real ( kind = 8 ) fnp
  real ( kind = 8 ) fs
  real ( kind = 8 ) fsp
  real ( kind = 8 ) fw
  real ( kind = 8 ) fwp
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ieqn
  integer ( kind = 4 ) option
  integer ( kind = 4 ) ivar
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) lambda
  real ( kind = 8 ) p13_gp
  real ( kind = 8 ) p13_gx
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) uc
  real ( kind = 8 ) ue
  real ( kind = 8 ) un
  real ( kind = 8 ) us
  real ( kind = 8 ) uw

  h = 1.0D+00 / real ( m + 1, kind = 8 )

  ieqn = 0

  do i = 1, m
    do j = 1, m

      ieqn = ( j - 1 ) * m + i

      uc = u(i,j)

      if ( i < m ) then
        un = u(i+1,j)
      else
        un = 0.0D+00
      end if

      if ( 1 < i ) then
        us = u(i-1,j)
      else
        us = 0.0D+00
      end if

      if ( j < m ) then
        ue = u(i,j+1)
      else
        ue = 0.0D+00
      end if

      if ( 1 < j ) then
        uw = u(i,j-1)
      else
        uw = 0.0D+00
      end if

      fc = p13_gx ( option, uc )
      fn = p13_gx ( option, un )
      fs = p13_gx ( option, us )
      fe = p13_gx ( option, ue )
      fw = p13_gx ( option, uw )

      del5f = fc + h * h * ( - 4.0D+00 * fc + fn + fs + fe + fw ) / 12.0D+00

      fcp = p13_gp ( option, uc )
      fnp = p13_gp ( option, un )
      fsp = p13_gp ( option, us )
      fep = p13_gp ( option, ue )
      fwp = p13_gp ( option, uw )

      ivar = ( j - 1 ) * m + i
      jac(ieqn,ivar) = - 20.0D+00 / ( 6.0D+00 * h * h ) &
        + lambda * ( fcp - 4.0D+00 * h * h * fcp / 12.0D+00 )

      if ( i < m ) then
        ivar = ( j - 1 ) * m + i + 1
        jac(ieqn,ivar) = 4.0D+00 / ( 6.0D+00 * h * h ) &
          + lambda * h * h * fnp / 12.0D+00
      end if

      if ( 1 < i ) then
        ivar = ( j - 1 ) * m + i - 1
        jac(ieqn,ivar) = 4.0D+00 / ( 6.0D+00 * h * h ) &
          + lambda * h * h * fsp / 12.0D+00
      end if

      if ( j < m ) then
        ivar = j * m + i
        jac(ieqn,ivar) = 4.0D+00 / ( 6.0D+00 * h * h ) &
          + lambda * h * h * fep / 12.0D+00
      end if

      if ( 1 < j ) then
        ivar = ( j - 2 ) * m + i
        jac(ieqn,ivar) = 4.0D+00 / ( 6.0D+00 * h * h ) &
          + lambda * h * h * fwp / 12.0D+00
      end if

      if ( 1 < i .and. 1 < j ) then
        ivar = ( j - 2 ) * m + i - 1
        jac(ieqn,ivar) = 1.0D+00 / ( 6.0D+00 * h * h )
      end if

      if ( 1 < i .and. j < m ) then
        ivar = j * m + i - 1
        jac(ieqn,ivar) = 1.0D+00 / ( 6.0D+00 * h * h )
      end if

      if ( i < m .and. 1 < j ) then
        ivar = ( j - 2 ) * m + i + 1
        jac(ieqn,ivar) = 1.0D+00 / ( 6.0D+00 * h * h )
      end if

      if ( i < m .and. j < m ) then
        ivar = j * m + i + 1
        jac(ieqn,ivar) = 1.0D+00 / ( 6.0D+00 * h * h )
      end if

      ivar = nvar
      jac(ieqn,nvar) = del5f

    end do
  end do

  return
end
subroutine p13_nvar ( option, nvar )

!*****************************************************************************80
!
!! P13_NVAR sets the number of variables for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nvar

  if ( option == 1 ) then
    m = 8
  else if ( option == 2 ) then
    m = 8
  else if ( option == 3 ) then
    m = 12
  else if ( option == 4 ) then
    m = 12
  else if ( option == 5 ) then
    m = 16
  else if ( option == 6 ) then
    m = 16
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P13_NVAR - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  nvar = m * m + 1

  return
end
subroutine p13_option_num ( option_num )

!*****************************************************************************80
!
!! P13_OPTION_NUM returns the number of options for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 6

  return
end
subroutine p13_start ( option, nvar, x )

!*****************************************************************************80
!
!! P13_START returns a starting point for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1:nvar) = 0.0D+00

  return
end
subroutine p13_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P13_STEPSIZE returns step sizes for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =     2.000D+00
  hmin =  0.001D+00
  hmax = 10.000D+00

  return
end
subroutine p13_title ( option, title )

!*****************************************************************************80
!
!! P13_TITLE sets the title for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  if ( option == 1 ) then
    title = 'Simpson''s BVP, F(U) = EXP(U), M = 8.'
  else if ( option == 2 ) then
    title = 'Simpson''s BVP, F(U) = function 2, M = 8.'
  else if ( option == 3 ) then
    title = 'Simpson''s BVP, F(U) = EXP(U), M = 12.'
  else if ( option == 4 ) then
    title = 'Simpson''s BVP, F(U) = function 2, M = 12.'
  else if ( option == 5 ) then
    title = 'Simpson''s BVP, F(U) = EXP(U), M = 16.'
  else if ( option == 6 ) then
    title = 'Simpson''s BVP, F(U) = function 2, M = 16.'
  else
    title = '???'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P13_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
    stop
  end if

  return
end
function p14_fu ( lambda, u )

!*****************************************************************************80
!
!! P14_FU computes the auxilliary function F(LAMBDA,U).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LAMBDA, U, the arguments of the function.
!
!    Output, real ( kind = 8 ) P14_FU, the value of the function.
!
  implicit none

  real ( kind = 8 ), parameter :: alpha = 0.1D+00
  real ( kind = 8 ) lambda
  real ( kind = 8 ) p14_fu
  real ( kind = 8 ) u

  p14_fu = 1.0D+00 + lambda / ( u + alpha )**2

  return
end
function p14_fudl ( u )

!*****************************************************************************80
!
!! P14_FUDL computes d F(LAMBDA,U) / d LAMBDA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U, the argument of the function.
!
!    Output, real ( kind = 8 ) P14_FUDL, the value of the derivative
!    of the function with respect to LAMBDA.
!
  implicit none

  real ( kind = 8 ), parameter :: alpha = 0.1D+00
  real ( kind = 8 ) p14_fudl
  real ( kind = 8 ) u

  p14_fudl = 1.0D+00 / ( u + alpha )**2

  return
end
function p14_fudu ( lambda, u )

!*****************************************************************************80
!
!! P14_FUDU computes d F(LAMBDA,U) / d U
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LAMBDA, U, the arguments of the function.
!
!    Output, real ( kind = 8 ) P14_FUDU, the value of the derivative
!    of the function with respect to U.
!
  implicit none

  real ( kind = 8 ), parameter :: alpha = 0.1D+00
  real ( kind = 8 ) lambda
  real ( kind = 8 ) p14_fudu
  real ( kind = 8 ) u

  p14_fudu = - 2.0D+00 * lambda / ( u + alpha )**3

  return
end
subroutine p14_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P14_FUN computes the function for problem 14.
!
!  Discussion:
!
!    Keller's boundary value problem.
!
!    The continuous problem is a two point boundary value problem
!    describing a diffusion-kinetics system, of the form:
!
!      -d/dt ( t * t * F(U) * dU/dt ) + t * t * G(U) = 0
!
!    where F(U) and G(U) are given functions,
!    with boundary conditions
!
!      dU/dt(0) = 0,
!      U(1) = 1.
!
!    A finite difference approximation to this continous problem
!    is used.
!
!    M points T(I) are used.  With a spacing of H=1/(M-2), the points
!    are set so that
!
!      T(1)=-H, T(2)=0, T(3)=H, ..., T(M)=1.0D+00
!
!    First equation:
!
!      U(3) - U(1) = 0.0D+00
!
!    Equations I = 2 through I = M-1
!
!      TL**2 * F(UL) * ( U(I) - U(I-1) ) +
!      TR**2 * F(UR) * ( U(I) - U(I+1) ) +
!      H**2 * T**2 * G(U) = 0.0D+00
!
!    with
!
!      T  = T(I) = ( I - 2 ) * H
!      U  = U(I)
!      TL = 0.5 * ( T(I-1) + T(I) )
!      UL = 0.5 * ( U(I-1) + U(I) )
!      TR = 0.5 * ( T(I) + T(I+1) )
!      UR = 0.5 * ( U(I) + U(I+1) )
!
!    and the diffusion function F(U)
!
!      F(U) = 1 + LAMBDA / ( U + ALPHA )**2
!
!    and
!
!      G(U) = U / ( BETA * ( U + GAMMA ) )
!
!    Equation M-1:
!
!      U(M) = 1.0D+00
!
!    For this version ALPHA = BETA = GAMMA = 0.1.
!
!    The only choice for options is
!
!    OPTION = 1:
!      IT = NVAR,
!      XIT = 1.0,
!      LIM = NVAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Herbert Keller,
!    Numerical Methods for Two-point Boundary Value Problems,
!    Dover, 1992,
!    ISBN: 0486669254,
!    LC: QA372.K42.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar-1)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) m
  real ( kind = 8 ) p14_fu
  real ( kind = 8 ) p14_gu
  real ( kind = 8 ) t
  real ( kind = 8 ) tl
  real ( kind = 8 ) tr
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) x(nvar)

  m = nvar - 1
  h = 1.0D+00 / real ( m - 2, kind = 8 )
  lambda = x(nvar)

  fx(1) = x(3) - x(1)

  do i = 2, m - 1

    t = ( i - 2 ) * h
    tl = ( real ( i, kind = 8 ) - 2.5D+00 ) * h
    tr = ( real ( i, kind = 8 ) - 1.5D+00 ) * h
    ul = 0.5D+00 * ( x(i-1) + x(i) )
    ur = 0.5D+00 * ( x(i) + x(i+1) )

    fx(i) =   ( x(i) - x(i-1) ) * tl * tl * p14_fu ( lambda, ul ) &
            + ( x(i) - x(i+1) ) * tr * tr * p14_fu ( lambda, ur ) &
            + h * h * t * t * p14_gu ( x(i) )

  end do

  fx(m) = x(m) - 1.0D+00

  return
end
function p14_gu ( u )

!*****************************************************************************80
!
!! P14_GU computes the auxilliary function G(U).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U, the argument of the function.
!
!    Output, real ( kind = 8 ) P14_GU, the value of the function.
!
  implicit none

  real ( kind = 8 ), parameter :: beta = 0.1D+00
  real ( kind = 8 ), parameter :: gamma = 0.1D+00
  real ( kind = 8 ) p14_gu
  real ( kind = 8 ) u

  p14_gu = u / ( beta * ( u + gamma ) )

  return
end
function p14_gudu ( u )

!*****************************************************************************80
!
!! P14_GUDU computes d G(U) / d U.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U, the argument of the function.
!
!    Output, real ( kind = 8 ) P14_GUDU, the value of the function.
!
  implicit none

  real ( kind = 8 ), parameter :: beta = 0.1D+00
  real ( kind = 8 ), parameter :: gamma = 0.1D+00
  real ( kind = 8 ) p14_gudu
  real ( kind = 8 ) u

  p14_gudu = gamma / ( beta * ( u + gamma )**2 )

  return
end
subroutine p14_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P14_JAC computes the jacobian of problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) m
  real ( kind = 8 ) p14_fu
  real ( kind = 8 ) p14_fudl
  real ( kind = 8 ) p14_fudu
  real ( kind = 8 ) p14_gudu
  real ( kind = 8 ) t
  real ( kind = 8 ) tl
  real ( kind = 8 ) tr
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) x(nvar)

  jac(1:nvar,1:nvar) = 0.0D+00

  m = nvar - 1
  h = 1.0D+00 / real ( m - 2, kind = 8 )
  lambda = x(nvar)
!
!  First equation.
!
  jac(1,1) = - 1.0D+00
  jac(1,3) =   1.0D+00
!
!  Intermediate equations.
!
  do i = 2, m - 1

    t = ( i - 2 ) * h
    tl = ( real ( i, kind = 8 ) - 2.5D+00 ) * h
    tr = ( real ( i, kind = 8 ) - 1.5D+00 ) * h
    ul = 0.5D+00 * ( x(i-1) + x(i) )
    ur = 0.5D+00 * ( x(i) + x(i+1) )

    jac(i,i) =   tl * tl * p14_fu ( lambda, ul ) &
      + ( x(i) - x(i-1) ) * tl * tl * p14_fudu ( lambda, ul ) * 0.5D+00 &
      + tr * tr * p14_fu ( lambda, ur ) &
      + ( x(i) - x(i+1) ) * tr * tr * p14_fudu ( lambda, ur ) * 0.5D+00 &
      + h * h * t * t * p14_gudu ( x(i) )

    jac(i,i-1) = - tl * tl * p14_fu ( lambda, ul ) &
      + ( x(i) - x(i-1) ) * tl * tl * p14_fudu ( lambda, ul ) * 0.5D+00

    jac(i,i+1) =  - tr * tr * p14_fu ( lambda, ur ) &
      + ( x(i) - x(i+1) ) * tr * tr * p14_fudu ( lambda, ur ) * 0.5D+00

    jac(i,nvar) = ( x(i) - x(i-1) ) * tl * tl * p14_fudl ( ul ) &
      + ( x(i) - x(i+1) ) * tr * tr * p14_fudl ( ur )

  end do
!
!  Last equation.
!
  jac(m,m) = 1.0D+00

  return
end
subroutine p14_nvar ( option, nvar )

!*****************************************************************************80
!
!! P14_NVAR sets the number of variables for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nvar

  m = 12
  nvar = m + 1

  return
end
subroutine p14_option_num ( option_num )

!*****************************************************************************80
!
!! P14_OPTION_NUM returns the number of options for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 1

  return
end
subroutine p14_start ( option, nvar, x )

!*****************************************************************************80
!
!! P14_START returns a starting point for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(nvar)

  if ( option == 1 ) then

    x(1)  = 0.029742673007439D+00
    x(2)  = 0.029742673007439D+00
    x(3)  = 0.029742673007439D+00
    x(4)  = 0.039933250735582D+00
    x(5)  = 0.061866539016825D+00
    x(6)  = 0.101137641789028D+00
    x(7)  = 0.164623875371221D+00
    x(8)  = 0.258536575943466D+00
    x(9)  = 0.387217701462343D+00
    x(10) = 0.553103336509555D+00
    x(11) = 0.757271228030916D+00
    x(12) = 1.000000000000000D+00
    x(13) = 0.000000000000000D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P14_START - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized option value = ', option
  end if

  return
end
subroutine p14_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P14_STEPSIZE returns step sizes for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =     2.000D+00
  hmin =  0.001D+00
  hmax = 10.000D+00

  return
end
subroutine p14_title ( option, title )

!*****************************************************************************80
!
!! P14_TITLE sets the title for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  title = 'Keller''s BVP.'

  return
end
subroutine p15_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P15_FUN evaluates the function for problem 15.
!
!  Title:
!
!    The Trigger Circuit.
!
!  Description:
!
!    The current flow of a trigger circuit with an operational amplifier
!    is modeled.  The variables are voltages, with X(6) the output
!    voltage and X(7) the input voltage.
!
!    The function has the form
!
!      F(X) = A * X + PHI ( X )
!
!    where A is a 6 by 7 matrix, and PHI is a nonlinear term, that is,
!
!      F(I) = SUM ( 1 <= J <= 7 ) A(I,J) * X(J) + PHI ( X )
!
!  Options:
!
!    Melhem lists the following limit points in X(7):
!
!    ( 0.04936  0.54735  0.04944  0.04944  0.12920  1.16602  0.60185 )
!    ( 0.23577  0.66296  0.23759  0.23760  0.62083  9.60913  0.32286 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 20008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Rami Melhem, Werner Rheinboldt,
!    A Comparison of Methods for Determining Turning Points of Nonlinear Equations,
!    Computing,
!    Volume 29, Number 3, September 1982, pages 201-226.
!
!    Gerd Poenisch, Hubert Schwetlick,
!    Computing Turning Points of Curves Implicitly Defined by Nonlinear
!    Equations Depending on a Parameter,
!    Computing,
!    Volume 26, Number 2, June 1981, pages 107-121.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the point of evaluation.
!
!    Input, real ( kind = 8 ) FX(NVAR-1), the function value.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) array(6,7)
  real ( kind = 8 ) fx(nvar-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(nvar)
!
!  Get the linear coefficients.
!
  call p15_gx ( array )
!
!  Compute the linear portion of the function.
!
  fx(1:nvar-1) = 0.0D+00

  do i = 1, nvar - 1
    do j = 1, nvar
      fx(i) = fx(i) + array(i,j) * x(j)
    end do
  end do
!
!  Add the nonlinear terms.
!
  fx(2) = fx(2) + 5.6D-08 * ( exp ( 25.0D+00 * x(2) ) - 1.0D+00 )
  fx(5) = fx(5) + 5.6D-08 * ( exp ( 25.0D+00 * x(5) ) - 1.0D+00 )
  fx(6) = fx(6) - 7.65D+00 * atan ( 1962.0D+00 * ( x(3) - x(1) ) ) / 0.201D+00

  return
end
subroutine p15_gx ( array )

!*****************************************************************************80
!
!! P15_GX returns the coefficients of the linear portion of the function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ARRAY(6,7), the coefficients of the linear portion
!    of the function, which are sums of the inverses of resistances.
!
  implicit none

  real ( kind = 8 ) array(6,7)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: s0 = 1.0D+00 / 10.0D+00
  real ( kind = 8 ), parameter :: s1 = 1.0D+00 / 39.0D+00
  real ( kind = 8 ), parameter :: s2 = 1.0D+00 / 51.0D+00
  real ( kind = 8 ), parameter :: s3 = 1.0D+00 / 10.0D+00
  real ( kind = 8 ), parameter :: s4 = 1.0D+00 / 25.5D+00
  real ( kind = 8 ), parameter :: s5 = 1.0D+00 / 1.0D+00
  real ( kind = 8 ), parameter :: s6 = 1.0D+00 / 0.62D+00
  real ( kind = 8 ), parameter :: s7 = 1.0D+00 / 13.0D+00
  real ( kind = 8 ), parameter :: s8 = 1.0D+00 / 0.201D+00

  array(1:6,1:7) = 0.0D+00

  array(1,1) = + s0 + s1 + s2
  array(1,2) =      - s1
  array(1,3) = - s0
  array(1,7) =           - s2

  array(2,1) =      - s1
  array(2,2) =      + s1 + s2
  array(2,6) =                - s3

  array(3,1) = - s0
  array(3,3) = + s0                + s4
  array(3,4) =                     - s4

  array(4,3) =                     - s4
  array(4,4) =                     + s4 + s5 + s6
  array(4,5) =                          - s5

  array(5,4) =                          - s5
  array(5,5) =                          + s5      + s7
  array(5,6) =                                    - s7

  array(6,2) =                - s3
  array(6,5) =                                    - s7
  array(6,6) =                + s3                + s7 + s8

  return
end
subroutine p15_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P15_JAC computes the jacobian for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) array(6,7)
  real ( kind = 8 ) jac(nvar,nvar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) u
  real ( kind = 8 ) x(nvar)

  jac(1:nvar,1:nvar) = 0.0D+00
!
!  Get the coefficients of the linear part of the function.
!
  call p15_gx ( array )

  jac(1:nvar-1,1:nvar) = array(1:nvar-1,1:nvar)
!
!  Add the derivatives of the nonlinear part.
!
  jac(2,2) = jac(2,2) + 5.6D-08 * 25.0D+00 * exp ( 25.0D+00 * x(2) )

  jac(5,5) = jac(5,5) + 5.6D-08 * 25.0D+00 * exp ( 25.0D+00 * x(5) )

  u = 1962.0D+00 * ( x(3) - x(1) )

  jac(6,1) = jac(6,1) + 7.65D+00 * 1962.0D+00 / ( 1.0D+00 + u * u ) / 0.201D+00

  jac(6,3) = jac(6,3) - 7.65D+00 * 1962.0D+00 / ( 1.0D+00 + u * u ) / 0.201D+00

  return
end
subroutine p15_nvar ( option, nvar )

!*****************************************************************************80
!
!! P15_NVAR sets the number of variables for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 7

  return
end
subroutine p15_option_num ( option_num )

!*****************************************************************************80
!
!! P15_OPTION_NUM returns the number of options for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 1

  return
end
subroutine p15_start ( option, nvar, x )

!*****************************************************************************80
!
!! P15_START returns a starting point for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1:nvar) = 0.0D+00

  return
end
subroutine p15_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P15_STEPSIZE returns step sizes for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    0.300D+00
  hmin = 0.001D+00
  hmax = 0.600D+00

  return
end
subroutine p15_title ( option, title )

!*****************************************************************************80
!
!! P15_TITLE sets the title for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  title = 'The Trigger Circuit.'

  return
end
subroutine p16_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P16_FUN evaluates the function for problem 16.
!
!  Title:
!
!    The Moore-Spence Chemical Reaction Integral Equation
!
!  Description:
!
!    The continuous equation describes the heat and mass transfer in a
!    plate-shaped porous catalyst, and is of the form
!
!      Y(S) = 1 + integral ( 0 <= T <= 1) K(S,T) * G(Y(T)) dT
!
!    with
!
!      K(S,T) = MAX ( S, T ) - 1
!
!      G(Y) = Y * EXP ( BETA * GAMMA * ( 1 - Y )
!        / ( 1 + BETA * ( 1 - Y ) ) )
!
!    with
!
!      BETA = 0.4,
!      GAMMA = 20.0.
!
!    The integral equation is discretized using M equally spaced
!    abscissas T(I) from 0 to 1, and applying the trapezoidal rule to
!    compute the integrand.  Finally, the integral is multiplied
!    by a homotopy parameter LAMBDA so that the problem is easy
!    to solve for LAMBDA = 0, while the solution for LAMBDA = 1
!    is the solution of the original problem.  Thus:
!
!      F(I) = Y(I) - 1 - LAMBDA * Trapezoid ( K(S(I),T())*G(Y(T())), T() ).
!
!  Options:
!
!    The solution for LAMBDA = 1 is desired.
!
!    With NVAR = 17, Melhem lists the limit points in LAMBDA:
!
!      LAMBDA = 0.1375390,  x(16) = 0.8524311,
!      LAMBDA = 0.07791579, x(16) = 0.4657826.
!
!    Computational results with this program are:
!
!      LAMBDA = 0.1286312,  x(16) = 0.8977113,
!      LAMBDA = 0.0926850,  x(16) = 0.2956740.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Rami Melhem, Werner Rheinboldt,
!    A Comparison of Methods for Determining Turning Points of Nonlinear Equations,
!    Computing,
!    Volume 29, Number 3, September 1982, pages 201-226.
!
!    Gerald Moore, Alastair Spence,
!    The Calculation of Turning Points of Nonlinear Equations,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 4, August 1980, pages 567-576.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at U.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) arg
  real ( kind = 8 ), parameter :: beta = 0.4D+00
  real ( kind = 8 ) factor
  real ( kind = 8 ) fx(nvar-1)
  real ( kind = 8 ), parameter :: gamma = 20.0D+00
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) m
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) trapezoid
  real ( kind = 8 ) x(nvar)

  lambda = x(nvar)
  m = nvar - 1
  h = 1.0D+00 / real ( m - 1, kind = 8 )

  do i = 1, m

    s = h * ( i - 1 )
    trapezoid = 0.0D+00

    do j = 1, m

      t = h * ( j - 1 )

      arg = beta * gamma * ( 1.0D+00 - x(j) ) &
        / ( 1.0D+00 + beta * ( 1.0D+00 - x(j) ) )

      if ( j == 1 ) then
        factor = 0.5D+00
      else if ( j < m - 1 ) then
        factor = 1.0D+00
      else if ( j == m ) then
        factor = 0.5D+00
      end if

      trapezoid = trapezoid + h * factor * x(j) * exp ( arg ) * &
        ( max ( s, t ) - 1.0D+00 )

    end do

    fx(i) = x(i) - 1.0D+00 - lambda * trapezoid

  end do

  return
end
subroutine p16_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P16_JAC computes the jacobian for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) U(NVAR), the point where the jacobian is evaluated.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) arg
  real ( kind = 8 ), parameter :: beta = 0.4D+00
  real ( kind = 8 ) dg
  real ( kind = 8 ) factor
  real ( kind = 8 ), parameter :: gamma = 20.0D+00
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) m
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) trapezoid
  real ( kind = 8 ) x(nvar)

  jac(1:nvar,1:nvar) = 0.0D+00

  lambda = x(nvar)
  m = nvar - 1
  h = 1.0D+00 / real ( m - 1, kind = 8 )

  do i = 1, m

    s = h * ( i - 1 )
    trapezoid = 0.0D+00

    do j = 1, m

      t = h * ( j - 1 )

      arg = beta * gamma * ( 1.0D+00 - x(j) ) &
        / ( 1.0D+00 + beta * ( 1.0D+00 - x(j) ) )

      if ( j == 1 ) then
        factor = 0.5D+00
      else if ( j < m - 1 ) then
        factor = 1.0D+00
      else if ( j == m ) then
        factor = 0.5D+00
      end if

      trapezoid = trapezoid + h * factor * x(j) * exp ( arg ) * &
        ( max ( s, t ) - 1.0D+00 )

      dg = - beta * gamma / ( 1.0D+00 + beta * ( 1.0D+00 - x(j) ) )**2

      jac(i,j) = jac(i,j) - lambda * h * factor * exp ( arg ) * &
        ( 1.0D+00 + x(j) * dg ) * ( max ( s, t ) - 1.0D+00 )

    end do

    jac(i,i) = jac(i,i) + 1.0D+00

    jac(i,nvar) = - trapezoid

  end do

  return
end
subroutine p16_nvar ( option, nvar )

!*****************************************************************************80
!
!! P16_NVAR sets the number of variables for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 17

  return
end
subroutine p16_option_num ( option_num )

!*****************************************************************************80
!
!! P16_OPTION_NUM returns the number of options for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 1

  return
end
subroutine p16_start ( option, nvar, x )

!*****************************************************************************80
!
!! P16_START returns a starting point for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1:nvar-1) = 1.0D+00
  x(nvar) = 0.0D+00

  return
end
subroutine p16_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P16_STEPSIZE returns step sizes for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    0.200D+00
  hmin = 0.001D+00
  hmax = 2.000D+00

  return
end
subroutine p16_title ( option, title )

!*****************************************************************************80
!
!! P16_TITLE sets the title for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  title = 'The Moore Spence Chemical Reaction Integral Equation.'

  return
end
subroutine p17_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P17_FUN evaluates the function for problem 17.
!
!  Title:
!
!    The Bremermann Propane Combustion System
!
!  Description:
!
!    The equations describe the combustion of propane (C3H4) in air
!    (O2 and N2) at 2200 degrees Fahrenheit.  The chemical substances
!    monitored include:
!
!      X(1)  CO2  carbon dioxide
!      X(2)  H2O  water
!      X(3)  N2
!      X(4)  CO   carbon monoxide
!      X(5)  H2
!      X(6)  H
!      X(7)  OH
!      X(8)  O
!      X(9)  NO
!      X(10) O2
!
!    with auxilliary variables
!
!      X(11)  amount of air: 0.5*X(11) moles of O2, 2*X(11) moles of N2.
!      X(12)  air pressure in atmospheres.
!
!    The mass balance and reaction equations become, once square
!    roots are eliminated:
!
!      F(1) = X(1) + X(4) - 3.0D+00
!      F(2) = 2 * X(1) + X(2) + X(4) + X(7) + X(8) + X(9)
!             + 2 * X(10) - X(12)
!      F(3) = 2 * X(2) + 2 * X(5) + X(6) + X(7) - 8.0D+00
!      F(4) = 2 * X(3) + X(9) - 4 * X(12)
!      F(5) = X(1) * X(5) - 0.193 * X(2) * X(4)
!      F(6) = X(11) * X(1) * X(6)**2 - 0.002597**2 * X(2) * X(4) * XSUM
!      F(7) = X(11) * X(4) * X(7)**2 - 0.003448**2 * X(1) * X(2) * XSUM
!      F(8) = X(11) * X(4) * X(8) - 0.0001799 * X(1) * XSUM
!      F(9) = X(11) * X(4)**2 * X(9)**2
!             - 0.0002155**2 * X(1)**2 * X(3) * XSUM
!      F(10)= X(11) * X(4)**2 * X(10) - 0.00003846 * X(1)**2 * XSUM
!
!    where
!
!      XSUM = Sum ( 1 <= I <= 10 ) X(I)
!
!  Options:
!
!    OPTION = 1:
!
!      FX(11) = X(11) - 1.0D+00 (fixed concentration)
!
!    OPTION = 2:
!
!      FX(11) = X(12) - 5.0D+00 (fixed pressure)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hans Bremermann,
!    Calculation of Equilibrium Points for Models of Ecological and
!    Chemical Systems,
!    in Proceedings of a Conference on the Applications of Undergraduate
!    Mathematics in the Engineering, Life, Managerial and Social Sciences,
!    Georgia Institute of Technology, June 1973, pages 198-217.
!
!    K L Hiebert,
!    A Comparison of Software Which Solves Systems of Nonlinear Equations,
!    Technical Report SAND-80-0181, 1980,
!    Sandia National Laboratory, Albuquerque, New Mexico,
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) U(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at U.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xsum

  xsum = sum ( x(1:10) )

  fx(1) = x(1) + x(4) - 3.0D+00

  fx(2) = 2.0D+00 * x(1) + x(2) + x(4) + x(7) + x(8) + x(9) &
    + 2.0D+00 * x(10) - x(12)

  fx(3) = 2.0D+00 * x(2) + 2.0D+00 * x(5) + x(6) + x(7) - 8.0D+00

  fx(4) = 2.0D+00 * x(3) + x(9) - 4.0D+00 * x(12)

  fx(5) = x(1) * x(5) - 0.193D+00 * x(2) * x(4)

  fx(6) = x(11) * x(1) * x(6) * x(6) - 0.002597D+00**2 * x(2) * x(4) * xsum

  fx(7) = x(11) * x(4) * x(7) * x(7) - 0.003448D+00**2 * x(1) * x(2) * xsum

  fx(8) = x(11) * x(4) * x(8) - 0.0001799D+00 * x(1) * xsum

  fx(9) = x(11) * x(4) * x(4) * x(9) * x(9) &
    - 0.0002155D+00**2 * x(1) * x(1) * x(3) * xsum

  fx(10) = x(11) * x(4) * x(4) * x(10) - 0.00003846D+00 * x(1) * x(1) * xsum

  if ( option == 1 ) then
    fx(11) = x(11) - 1.0D+00
  else if ( option == 2 ) then
    fx(11) = x(12) - 5.0D+00
  end if

  return
end
subroutine p17_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P17_JAC evaluates the jacobian for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) U(NVAR), the point where the jacobian is evaluated.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) term
  real ( kind = 8 ) x(nvar)
  real ( kind = 8 ) xsum

  jac(1:nvar,1:nvar) = 0.0D+00

  xsum = sum ( x(1:10) )

  jac(1,1) =    1.0D+00
  jac(1,4) =    1.0D+00

  jac(2,1) =    2.0D+00
  jac(2,2) =    1.0D+00
  jac(2,4) =    1.0D+00
  jac(2,7) =    1.0D+00
  jac(2,8) =    1.0D+00
  jac(2,9) =    1.0D+00
  jac(2,10) =   2.0D+00
  jac(2,12) = - 1.0D+00

  jac(3,2) =    2.0D+00
  jac(3,5) =    2.0D+00
  jac(3,6) =    1.0D+00
  jac(3,7) =    1.0D+00

  jac(4,3) =    2.0D+00
  jac(4,9) =    1.0D+00
  jac(4,12) = - 4.0D+00

  jac(5,1) = x(5)
  jac(5,2) = - 0.193D+00 * x(4)
  jac(5,4) = - 0.193D+00 * x(2)
  jac(5,5) = x(1)

  term = - 0.002597D+00**2 * x(2) * x(4)

  jac(6,1) = x(6) * x(6) * x(11) + term
  jac(6,2) = - 0.002597D+00**2 * x(4) * ( xsum + x(2) )
  jac(6,3) = term
  jac(6,4) = - 0.002597D+00**2 * x(2) * ( xsum + x(4) )
  jac(6,5) = term
  jac(6,6) = term + 2.0D+00 * x(1) * x(6) * x(11)
  jac(6,7) = term
  jac(6,8) = term
  jac(6,9) = term
  jac(6,10) = term
  jac(6,11) = x(1) * x(6) * x(6)

  term = - 0.003448D+00**2 * x(1) * x(2)

  jac(7,1) = - 0.003448D+00**2 * x(2) * ( xsum + x(1) )
  jac(7,2) = - 0.003448D+00**2 * x(1) * ( xsum + x(2) )
  jac(7,3) = term
  jac(7,4) = x(7) * x(7) * x(11) + term
  jac(7,5) = term
  jac(7,6) = term
  jac(7,7) = 2.0D+00 * x(4) * x(7) * x(11) + term
  jac(7,8) = term
  jac(7,9) = term
  jac(7,10) = term
  jac(7,11) = x(4) * x(7) * x(7)

  term = - 0.0001799D+00 * x(1)

  jac(8,1) = - 0.0001799D+00 * ( xsum + x(1) )
  jac(8,2) = term
  jac(8,3) = term
  jac(8,4) = x(8) * x(11) + term
  jac(8,5) = term
  jac(8,6) = term
  jac(8,7) = term
  jac(8,8) = x(4) * x(11) + term
  jac(8,9) = term
  jac(8,10) = term
  jac(8,11) = x(4) * x(8)

  term = - 0.0002155D+00**2 * x(1) * x(1) * x(3)

  jac(9,1) = - 0.0002155D+00**2 * x(1) * x(3) * ( x(1) + 2.0D+00 * xsum )
  jac(9,2) = term
  jac(9,3) = - 0.0002155D+00**2 * x(1) * x(1) * ( x(3) + xsum )
  jac(9,4) = 2.0D+00 * x(4) * x(9) * x(9) * x(11) + term
  jac(9,5) = term
  jac(9,6) = term
  jac(9,7) = term
  jac(9,8) = term
  jac(9,9) = 2.0D+00 * x(4) * x(4) * x(9) * x(11) + term
  jac(9,10) = term
  jac(9,11) = x(4) * x(4) * x(9) * x(9)

  term = - 0.00003846D+00 * x(1) * x(1)

  jac(10,1) = - 0.00003846D+00 * x(1) * ( x(1) + 2.0D+00 * xsum )
  jac(10,2) = term
  jac(10,3) = term
  jac(10,4) = 2.0D+00 * x(4) * x(10) * x(11) + term
  jac(10,5) = term
  jac(10,6) = term
  jac(10,7) = term
  jac(10,8) = term
  jac(10,9) = term
  jac(10,10) = x(4) * x(4) * x(11) + term
  jac(10,11) = x(4) * x(4) * x(10)

  if ( option == 1 ) then
    jac(11,11) = 1.0D+00
  else if ( option == 2 ) then
    jac(11,12) = 1.0D+00
  end if

  return
end
subroutine p17_nvar ( option, nvar )

!*****************************************************************************80
!
!! P17_NVAR sets the number of variables for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 12

  return
end
subroutine p17_option_num ( option_num )

!*****************************************************************************80
!
!! P17_OPTION_NUM returns the number of options for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 2

  return
end
subroutine p17_start ( option, nvar, x )

!*****************************************************************************80
!
!! P17_START returns a starting point for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1) = 0.3564320D+00
  x(2) = 1.636071D+00
  x(3) = 9.999810D+00
  x(4) = 2.643568D+00
  x(5) = 2.341926D+00
  x(6) = 0.3732447D-01
  x(7) = 0.6681509D-02
  x(8) = 0.4128999D-03
  x(9) = 0.3790901D-03
  x(10) = 0.1190167D-04

  if ( option == 1 ) then
    x(11) = 1.0D+00
    x(12) = 5.0D+00
  else if ( option == 2 ) then
    x(11) = 1.0D+00
    x(12) = 5.0D+00
  end if

  return
end
subroutine p17_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P17_STEPSIZE returns step sizes for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    1.000D+00
  hmin = 0.001D+00
  hmax = 2.000D+00

  return
end
subroutine p17_title ( option, title )

!*****************************************************************************80
!
!! P17_TITLE sets the title for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  if ( option == 1 ) then
    title = 'Bremermann Propane Combustion System, fixed pressure.'
  else if ( option == 2 ) then
    title = 'Bremermann Propane Combustion System, fixed concentration.'
  else
    title = '???'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P17_TITLE - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized option number.'
    stop
  end if

  return
end
subroutine p18_fun ( option, nvar, u, fx )

!*****************************************************************************80
!
!! P18_FUN evaluates the function for problem 18.
!
!  Title:
!
!    The Semiconductor Problem.
!
!  Description:
!
!    The continuous problem is a two point boundary value problem
!    of the form
!
!      - U'' = G ( T, U, LAMBDA )
!
!    for A < T < B, with boundary conditions
!
!      U(A) = LAMBDA * UA,
!      U(B) = LAMBDA * UB.
!
!    and with right hand side:
!
!      G ( T, U, LAMBDA ) = LAMBDA *
!        ( CA * EXP ( LAMBDA * BETA * ( LAMBDA * UA - U ) )
!        - CB * EXP ( LAMBDA * BETA * ( U - LAMBDA * UB ) ) + H(T) )
!
!    and
!
!      H(T) = - CA for T <= 0,
!           =   CB for 0 < T.
!
!    The discrete version of the problem uses a mesh of M points
!    T(1) = A, T(2) = A + H, T(3) = A * 2 * H, ..., T(M) + B,
!    and corresponding solution values U(I).  The system of
!    M equations is:
!
!      U(1) = LAMBDA * UA
!
!      - U(I-1) + 2 * U(I) - U(I+1) = 2 * H * LAMBDA * G ( T(I), U(I), LAMBDA )
!
!      U(M) = LAMBDA * UB
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SJ Polak, A Wachten, H Vaes, A deBeer, Cor denHeijer,
!    A Continuation Method for the Calculation of Electrostatic
!    Potentials in Semiconductors,<br>
!    Technical Report ISA-TIS/CARD,<br>
!    NV Philips Gloeilampen-Fabrieken, 1979.
!
!    Cor denHeijer, Werner Rheinboldt,
!    On Steplength Algorithms for a Class of Continuation Methods,
!    SIAM Journal on Numerical Analysis,
!    Volume 18, Number 5, October 1981, pages 925-947.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) U(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at U.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ), parameter :: b = 0.010D+00
  real ( kind = 8 ) del2x
  real ( kind = 8 ) fx(nvar-1)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) m
  real ( kind = 8 ) p18_gx
  real ( kind = 8 ) t
  real ( kind = 8 ) u(nvar)
  real ( kind = 8 ), parameter :: ua = 0.0D+00
  real ( kind = 8 ), parameter :: ub = 25.0D+00

  lambda = u(nvar)
  m = nvar - 1
  h = 1.0D+00 / real ( m - 1, kind = 8 )
  fx(1) = u(1) - lambda * ua

  do i = 2, m - 1

    t = ( real ( m - i,     kind = 8 ) * a   &
        + real (     i - 1, kind = 8 ) * b ) &
        / real ( m     - 1, kind = 8 )

    del2x = ( u(i-1) - 2.0D+00 * u(i) + u(i+1) ) / ( 2.0D+00 * h )

    fx(i) = del2x - p18_gx ( t, u(i), lambda )

  end do

  fx(m) = u(m) - lambda * ub

  return
end
function p18_gpl ( t, u, lambda )

!*****************************************************************************80
!
!! P18_GPL evaluates d G ( T, U, LAMBDA ) / d LAMBDA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, U, LAMBDA, the arguments of the function.
!
!    Output, real ( kind = 8 ) P18_GPL, the derivative of the function with
!    respect to LAMBDA.
!
  implicit none

  real ( kind = 8 ), parameter :: beta = 20.0D+00
  real ( kind = 8 ), parameter :: ca = 1.0D+06
  real ( kind = 8 ), parameter :: cb = 1.0D+07
  real ( kind = 8 ) e1
  real ( kind = 8 ) e2
  real ( kind = 8 ) ht
  real ( kind = 8 ) lambda
  real ( kind = 8 ) p18_gpl
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ), parameter :: ua = 0.0D+00
  real ( kind = 8 ), parameter :: ub = 25.0D+00

  if ( t <= 0.0D+00 ) then
    ht = - ca
  else
    ht = cb
  end if

  e1 = exp ( lambda * beta * ( lambda * ua - u ) )
  e2 = exp ( lambda * beta * ( u - lambda * ub ) )

  p18_gpl = ht + ca * e1 - cb * e2 + lambda * &
       ( ca * beta * ( 2.0D+00 * lambda * ua - u ) * e1 &
       - cb * beta * ( u - 2.0D+00 * lambda * ub ) * e2 )

  return
end
function p18_gpu ( u, lambda )

!*****************************************************************************80
!
!! P18_GPU evaluates d G ( T, U, LAMBDA ) / dU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U, LAMBDA, the arguments of the function.
!
!    Output, real ( kind = 8 ) P18_GPU, the derivative of the function with
!    respect to U.
!
  implicit none

  real ( kind = 8 ), parameter :: beta = 20.0D+00
  real ( kind = 8 ), parameter :: ca = 1.0D+06
  real ( kind = 8 ), parameter :: cb = 1.0D+07
  real ( kind = 8 ) lambda
  real ( kind = 8 ) p18_gpu
  real ( kind = 8 ) u
  real ( kind = 8 ), parameter :: ua = 0.0D+00
  real ( kind = 8 ), parameter :: ub = 25.0D+00

  p18_gpu = - lambda * lambda * beta * ( &
       ca * exp ( lambda * beta * ( lambda * ua - u ) ) &
     + cb * exp ( lambda * beta * ( u - lambda * ub ) ) )

  return
end
function p18_gx ( t, u, lambda )

!*****************************************************************************80
!
!! P18_GX evaluates the auxilliary function G ( T, U, LAMBDA ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, U, LAMBDA, the arguments of the function.
!
!    Output, real ( kind = 8 ) P18_GX, the value of the function.
!
  implicit none

  real ( kind = 8 ), parameter :: beta = 20.0D+00
  real ( kind = 8 ), parameter :: ca = 1.0D+06
  real ( kind = 8 ), parameter :: cb = 1.0D+07
  real ( kind = 8 ) ht
  real ( kind = 8 ) lambda
  real ( kind = 8 ) p18_gx
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ), parameter :: ua = 0.0D+00
  real ( kind = 8 ), parameter :: ub = 25.0D+00

  if ( t <= 0.0D+00 ) then
    ht = - ca
  else
    ht = cb
  end if

  p18_gx = lambda * ( ht + ca * exp ( lambda * beta * ( lambda * ua - u ) ) &
                         - cb * exp ( lambda * beta * ( u - lambda * ub ) ) )

  return
end
subroutine p18_jac ( option, nvar, u, jac )

!*****************************************************************************80
!
!! P18_JAC evaluates the jacobian for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) U(NVAR), the point where the jacobian is evaluated.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ), parameter :: b = 0.01D+00
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ) p18_gpl
  real ( kind = 8 ) p18_gpu
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) m
  real ( kind = 8 ) t
  real ( kind = 8 ) u(nvar)
  real ( kind = 8 ), parameter :: ua = 0.0D+00
  real ( kind = 8 ), parameter :: ub = 25.0D+00

  jac(1:nvar,1:nvar) = 0.0D+00

  lambda = u(nvar)
  m = nvar - 1
  h = 1.0D+00 / real ( m - 1, kind = 8 )

  jac(1,1) = 1.0D+00
  jac(1,nvar) = - ua

  do i = 2, m - 1

    t = ( real ( m - i,     kind = 8 ) * a   &
        + real (     i - 1, kind = 8 ) * b ) &
        / real ( m     - 1, kind = 8 )

    jac(i,i-1) = 0.5D+00 / h
    jac(i,i) = - 1.0D+00 / h - p18_gpu ( u(i), lambda )
    jac(i,i+1) = 0.5D+00 / h
    jac(i,nvar) = - p18_gpl ( t, u(i), lambda )

  end do

  jac(m,m) = 1.0D+00
  jac(m,nvar) = - ub

  return
end
subroutine p18_nvar ( option, nvar )

!*****************************************************************************80
!
!! P18_NVAR sets the number of variables for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 12

  return
end
subroutine p18_option_num ( option_num )

!*****************************************************************************80
!
!! P18_OPTION_NUM returns the number of options for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 1

  return
end
subroutine p18_start ( option, nvar, x )

!*****************************************************************************80
!
!! P18_START returns a starting point for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1:nvar) = 0.0D+00

  return
end
subroutine p18_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P18_STEPSIZE returns step sizes for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    2.500D+00
  hmin = 0.001D+00
  hmax = 5.000D+00

  return
end
subroutine p18_title ( option, title )

!*****************************************************************************80
!
!! P18_TITLE sets the title for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  title = 'The Semiconductor Problem.'

  return
end
subroutine p19_con ( nvar, press, temper, x, con, flow )

!*****************************************************************************80
!
!! P19_CON returns physical constants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) PRESS, the pressure in atmospheres.
!
!    Input, real ( kind = 8 ) TEMPER, the temperature in degrees Kelvin.
!
!    Input, real ( kind = 8 ) X(NVAR), the point where the jacobian is evaluated.
!
!    Output, real ( kind = 8 ) CON(5), the equilibrium constants for the reagants.
!
!    Output, real ( kind = 8 ) FLOW(5), the flow rates for the reagants.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) arg
  real ( kind = 8 ) con(5)
  real ( kind = 8 ) flow(5)
  real ( kind = 8 ) press
  real ( kind = 8 ) temper
  real ( kind = 8 ) x(nvar)
!
!  Set flow rates.
!
  flow(1) = 10.0D+00
  flow(2) = 10.0D+00
  flow(3) = 10.0D+00
  flow(4) = 100.0D+00
  flow(5) = 100.0D+00
!
!  Set the equilibrium constants.
!
  con(1) = 1333.0D+00 / press

  con(2) = 33.0D+00 / press

  con(3) = 28780.0D+00 / press

  arg = 11.99D+00 - 4004.0D+00 / ( temper - 39.06 ) &
    - 8546.0D+00 * x(5) * x(5) / temper &
    + 4.0D+00 * x(5) * x(5) + 6754.0D+00 * x(5) * x(5) * x(4) / temper &
    - 8.0D+00 * x(5) * x(5) * x(4) - log ( press )

  con(4) = exp ( arg )

  arg = 10.98D+00 - 3362.0D+00 / ( temper - 50.79D+00 ) &
    - 2872.0D+00 * x(4) * x(4) / temper &
    - 6754.0D+00 * x(5) * x(4) * x(4) / temper &
    + 8.0D+00 * x(5) * x(4) * x(4) - log ( press )

  con(5) = exp ( arg )

  return
end
subroutine p19_conp ( con, nvar, temper, x, dc4dx4, dc4dx5, dc5dx4, dc5dx5 )

!*****************************************************************************80
!
!! P19_CONP returns physical constant derivatives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CON(5), the equilibrium constants for the reagants.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) TEMPER, the temperature in degrees Kelvin.
!
!    Input, real ( kind = 8 ) X(NVAR), the point where the jacobian is evaluated.
!
!    Output, real ( kind = 8 ) DC4DX4, DC4DX5, DC5DX4, DC5DX5, the values of
!    d CON(4)/d X(4), d CON(4)/d X(5), d CON(5)/d X(4) and d CON(5)/d X(5).
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) con(5)
  real ( kind = 8 ) dc4dx4
  real ( kind = 8 ) dc4dx5
  real ( kind = 8 ) dc5dx4
  real ( kind = 8 ) dc5dx5
  real ( kind = 8 ) temper
  real ( kind = 8 ) x(nvar)

  dc4dx4 = con(4) * ( 6754.0D+00 * x(5) * x(5) / temper - 8.0D+00 * x(5) * x(5) )

  dc4dx5 = con(4) * ( -8546.0D+00 * 2.0D+00 * x(5) / temper &
     + 8.0D+00 * x(5) + 6754.0D+00 * 2.0D+00 * x(4) * x(5) / temper &
     - 16.0D+00 * x(5) * x(4) )

  dc5dx4 = con(5) * ( -2872.0D+00 * 2.0D+00 * x(4) / temper &
     - 6754.0D+00 * 2.0D+00 * x(4) * x(5) / temper + 16.0D+00 * x(4) * x(5) )

  dc5dx5 = con(5) * ( -6754.0D+00 * x(4) * x(4) / temper + 8.0D+00 * x(4) * x(4) )

  return
end
subroutine p19_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P19_FUN evaluates the function for problem 19.
!
!  Title:
!
!    Nitric acid absorption problem
!
!  Description:
!
!    Physical Constants:
!
!      CON    - physical equilibrium constants for the five reagents.
!      FLOW   - flow rates for the five reagants in moles/hour.
!      PRESS  - pressure in atmospheres
!      TEMPER - temperature in Kelvin
!
!    Variables:
!
!      Entries 1 to 5 are the relative concentrations of liquid:
!
!        X(1)  = liquid NO2.
!        X(2)  = liquid N2O4.
!        X(3)  = liquid NO.
!        X(4)  = liquid H2O.
!        X(5)  = liquid HNO3.
!
!      Entries 6 through 10 are the relative concentrations of vapor:
!
!        X(6)  = vapor NO2.
!        X(7)  = vapor N2O4.
!        X(8)  = vapor NO.
!        X(9)  = vapor H2O.
!        X(10) = vapor HNO3.
!
!      Entries 11 and 12 are the number of moles:
!
!        X(11) = moles of liquid.
!        X(12) = moles of vapor.
!
!      Entry 13 is LAMBDA:
!
!        X(13) = LAMBDA, flowrate multiplier.
!
!    Equations:
!
!      Mole Balance equations, I = 1 to 5:
!
!        FX(I) = X(11) * X(I) + X(12) * X(I+5) - X(13) * FLOW(I)
!
!      Liquid-Vapor Transfer equations, I = 6 to 10:
!
!        FX(I) =  X(I) - CON(I-5) * X(I-5)
!
!      Liquid and Vapor proportions add to 1:
!
!        FX(11) =  X(1) + X(2) + X(3) + X(4) + X(5) - 1.0D+00
!        FX(12) =  X(6) + X(7) + X(8) + X(9) + X(10) - 1.0D+00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Tama Copeman,
!    Air Products and Chemicals, Inc.
!    Box 538,
!    Allentown, Pennsylvania, 18105.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) con(5)
  real ( kind = 8 ) flow(5)
  real ( kind = 8 ) fx(nvar-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: press = 7.0D+00
  real ( kind = 8 ), parameter :: temper = 323.0D+00
  real ( kind = 8 ) x(nvar)
!
!  Get chemical constants.
!
  call p19_con ( nvar, press, temper, x, con, flow )
!
!  Evaluate the Mole Balance equations:
!
  do i = 1, 5
    fx(i) = x(11) * x(i) + x(12) * x(i+5) - x(13) * flow(i)
  end do
!
!  Evaluate the Liquid-Vapor Transfer equations:
!
  do i = 6, 10
    fx(i) = x(i) - con(i-5) * x(i-5)
  end do
!
!  Evaluate the Liquid and Vapor Proportion equations:
!
  fx(11) = x(1) + x(2) + x(3) + x(4) + x(5) - 1.0D+00
  fx(12) = x(6) + x(7) + x(8) + x(9) + x(10) - 1.0D+00

  return
end
subroutine p19_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P19_JAC evaluates the jacobian for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the point where the jacobian is evaluated.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) con(5)
  real ( kind = 8 ) dc4dx4
  real ( kind = 8 ) dc4dx5
  real ( kind = 8 ) dc5dx4
  real ( kind = 8 ) dc5dx5
  real ( kind = 8 ) flow(5)
  real ( kind = 8 ) jac(nvar,nvar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: press = 7.0D+00
  real ( kind = 8 ), parameter :: temper = 323.0D+00
  real ( kind = 8 ) x(nvar)

  jac(1:nvar,1:nvar) = 0.0D+00
!
!  Get chemical constants.
!
  call p19_con ( nvar, press, temper, x, con, flow )
!
!  Get derivatives of chemical constants.
!
  call p19_conp ( con, nvar, temper, x, dc4dx4, dc4dx5, dc5dx4, dc5dx5 )
!
!  Differentiate the Mole Balance equations:
!
  do i = 1, 5
    jac(i,i) = x(11)
    jac(i,i+5) = x(12)
    jac(i,11) = x(i)
    jac(i,12) = x(i+5)
    jac(i,13) = - flow(i)
  end do
!
!  Differentiate the Liquid-Vapor Transfer equations:
!
  do i = 6, 10
    jac(i,i) = 1.0D+00
    jac(i,i-5) = - con(i-5)
  end do

  jac(9,4) = jac(9,4) - dc4dx4 * x(4)
  jac(9,5) = jac(9,5) - dc4dx5 * x(5)

  jac(10,4) = jac(10,4) - dc5dx4 * x(4)
  jac(10,5) = jac(10,5) - dc5dx5 * x(5)
!
!  Differentiate the Liquid and Vapor Proportion equations:
!
  jac(11,1:5) = 1.0D+00
  jac(12,6:10) = 1.0D+00

  return
end
subroutine p19_nvar ( option, nvar )

!*****************************************************************************80
!
!! P19_NVAR sets the number of variables for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 13

  return
end
subroutine p19_option_num ( option_num )

!*****************************************************************************80
!
!! P19_OPTION_NUM returns the number of options for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 1

  return
end
subroutine p19_start ( option, nvar, x )

!*****************************************************************************80
!
!! P19_START returns a starting point for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) i
  integer ( kind = 4 ) option
  real ( kind = 8 ) x(nvar)

  x(1:13) =   (/ &
      0.00218216D+00, &
      0.03171126D+00, &
      0.00010562D+00, &
      0.48301846D+00, &
      0.48298250D+00, &
      0.41554567D+00, &
      0.14949595D+00, &
      0.43425476D+00, &
      0.00018983D+00, &
      0.00051379D+00, &
    207.02239583D+00, &
     22.97760417D+00, &
      1.00000000D+00 /)

  return
end
subroutine p19_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P19_STEPSIZE returns step sizes for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    0.125000D+00
  hmin = 0.015625D+00
  hmax = 4.000000D+00

  return
end
subroutine p19_title ( option, title )

!*****************************************************************************80
!
!! P19_TITLE sets the title for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  title = 'Nitric Acid Absorption Flash.'

  return
end
subroutine p20_fun ( option, nvar, x, fx )

!*****************************************************************************80
!
!! P20_FUN evaluates the function for problem 20.
!
!  Title:
!
!    The Buckling Spring
!
!  Description:
!
!    We are given three points A, B, and C.
!    A is at the origin (0,0).
!    B has coordinates (X,Y) with Y nonnegative, and the ray from A to B
!    makes an angle of THETA with the horizontal axis.
!    C is at the point (2*X,0).
!
!    A spring extends from A to B, and is normally of length 1,
!    and is currently of length L.
!    A spring extends from B to C, and is normally of length 1,
!    and is currently of length L.
!    A spring force is also exerted, which tends to draw the two
!    springs together, proportional to the angle between the two springs.
!
!    A vertical load MU is applied at point B (downward is positive).
!    A horizontal load LAMBDA is applied at point C (leftware is positive).
!    The spring force is applied perpendicularly to the axes of the two springs.
!
!    If we compute F(1), the force along the axis of one spring, and
!    F(2), the force perpendicular to the axis of one spring, we have that
!    F(L,THETA,LAMBDA,MU) is given by:
!
!    F(1) =  - 2 ( 1 - L ) + 2 * LAMBDA * cos ( THETA ) + MU * sin ( THETA )
!    F(2) =  0.5 * THETA - 2 * LAMBDA * L * sin ( THETA ) + MU * L * cos ( THETA )
!
!    The user must specify a third, augmenting equation, of the form
!
!    F(3) = X(HOLD_INDEX) - HOLD_VALUE.
!
!    Typically, HOLD_INDEX is 2, for the varlable THETA, and HOLD_VALUE is
!    an angle measured in radians and strictly between 0 and PI/2.
!
!    Another choice for HOLD_INDEX would be 1, for the variable L, with
!    HOLD_VALUE greater than 0.  Values of L less than 1 represent compressed
!    springs; values greater than 1 indicate extended springs.
!
!    L represents the current length of the springs.
!    THETA represents the angle that the springs make with the horizontal axis.
!    MU is a vertical load applied at the midpoint B.
!    LAMBDA is a horizontal load applied at the right endpoint C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Tim Poston, Ian Stewart,
!    Catastrophe Theory and its Applications,
!    Dover, 1996,
!    ISBN13: 978-0486692715,
!    LC: QA614.58.P66.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the function.
!
!    Output, real ( kind = 8 ) FX(NVAR-1), the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) nvar

  real ( kind = 8 ) fx(nvar-1)
  integer ( kind = 4 ) hold_index
  real ( kind = 8 ) hold_value
  real ( kind = 8 ) l
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu
  integer ( kind = 4 ) option
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(nvar)

  l      = x(1)
  theta  = x(2)
  lambda = x(3)
  mu     = x(4)

  fx(1) = - 2.0D+00 * ( 1.0D+00 - l ) &
          + 2.0D+00 * lambda * cos ( theta ) &
          + mu * sin ( theta )

  fx(2) = 0.5D+00 * theta &
          - 2.0D+00 * lambda * l * sin ( theta ) &
          + mu * l * cos ( theta )

  call p20_i4_get ( 'hold_index', hold_index )
  call p20_r8_get ( 'hold_value', hold_value )

  fx(3) = x(hold_index) - hold_value

  return
end
subroutine p20_i4_get ( name, value )

!*****************************************************************************80
!
!! P20_I4_GET returns the value of an integer parameter for problem 20.
!
!  Discussion:
!
!    The only legal input name is 'hold_index', which indicates the
!    index of the variable to be held fixed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NAME, the name of the variable.
!
!    Output, integer ( kind = 4 ) VALUE, the value of the variable.
!
  implicit none

  character ( len = *  ) name
  integer   ( kind = 4 ) value

  call p20_i4_store ( 'get', name, value )

  return
end
subroutine p20_i4_set ( name, value )

!*****************************************************************************80
!
!! P20_I4_SET sets the value of an integer parameter for problem 20.
!
!  Discussion:
!
!    The only legal input name is 'hold_index', which indicates the
!    index of the variable to be held fixed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NAME, the name of the variable.
!
!    Input, integer ( kind = 4 ) VALUE, the value of the variable.
!
  implicit none

  character ( len = *  ) name
  integer   ( kind = 4 ) value

  call p20_i4_store ( 'set', name, value )

  return
end
subroutine p20_i4_store ( action, name, value )

!*****************************************************************************80
!
!! P20_I4_STORE sets or gets the value of an integer parameter for problem 20.
!
!  Discussion:
!
!    The only legal input name is 'hold_index', which indicates the
!    index of the variable to be held fixed.  This variable is given
!    the default value of 2.
!
!    The only legal values of 'hold_index' are 1 and 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, either 'get' or 'set'.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!
!    Input/output, integer ( kind = 4 ) VALUE, the value of the variable.
!
  implicit none

  character ( len = *  ) action
  integer   ( kind = 4 ), save :: hold_index = 2
  character ( len = *  ) name
  logical                s_eqi
  integer   ( kind = 4 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'hold_index' ) ) then
      value = hold_index
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_I4_STORE - Fatal error!'
      write ( *, '(a)' ) '  Unexpected "get" of NAME = "' // name // '".'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'hold_index' ) ) then
      if ( value == 1 .or. value == 2 ) then
        hold_index = value
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P20_I4_STORE - Fatal error!'
        write ( *, '(a,i8)' ) '  Unacceptable value for HOLD_INDEX = ', value
        write ( *, '(a)' ) '  Acceptable values are 1 and 2.'
        stop
      end if
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_I4_STORE - Fatal error!'
      write ( *, '(a)' ) '  Unexpected "set" of NAME = "' // name // '".'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P20_I4_STORE - Fatal error!'
    write ( *, '(a)' ) '  Unexpected ACTION = "' // action // '".'
    stop
  end if

  return
end
subroutine p20_jac ( option, nvar, x, jac )

!*****************************************************************************80
!
!! P20_JAC evaluates the jacobian for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Input, real ( kind = 8 ) X(NVAR), the argument of the jacobian.
!
!    Output, real ( kind = 8 ) JAC(NVAR,NVAR), the jacobian matrix evaluated
!    at X.  The NVAR-th row is not set by this routine.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) hold_index
  real ( kind = 8 ) jac(nvar,nvar)
  real ( kind = 8 ) l
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu
  integer ( kind = 4 ) option
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(nvar)

  jac(1:nvar,1:nvar) = 0.0D+00

  l      = x(1)
  theta  = x(2)
  lambda = x(3)
  mu     = x(4)

  jac(1,1) = 2.0D+00
  jac(1,2) = - 2.0D+00 * lambda * sin ( theta ) + mu * cos ( theta )
  jac(1,3) = 2.0D+00 * cos ( theta )
  jac(1,4) = sin ( theta )

  jac(2,1) = - 2.0D+00 * lambda * sin ( theta ) + mu * cos ( theta )
  jac(2,2) = 0.5D+00 - 2.0D+00 * lambda * l * cos ( theta ) &
           - mu * l * sin ( theta )
  jac(2,3) = - 2.0D+00 * l * sin ( theta )
  jac(2,4) = l * cos ( theta )

  call p20_i4_get ( 'hold_index', hold_index )

  jac(3,hold_index) = 1.0D+00

  return
end
subroutine p20_nvar ( option, nvar )

!*****************************************************************************80
!
!! P20_NVAR sets the number of variables for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option chosen for this problem.
!    For some problems, several options are available.  At least,
!    OPTION = 1 is always legal.
!
!    Output, integer ( kind = 4 ) NVAR, the number of variables.
!
  implicit none

  integer ( kind = 4 ) option
  integer ( kind = 4 ) nvar

  nvar = 4

  return
end
subroutine p20_option_num ( option_num )

!*****************************************************************************80
!
!! P20_OPTION_NUM returns the number of options for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OPTION_NUM, the number of options.
!
  implicit none

  integer ( kind = 4 ) option_num

  option_num = 1

  return
end
subroutine p20_r8_get ( name, value )

!*****************************************************************************80
!
!! P20_R8_GET returns the value of a real parameter for problem 20.
!
!  Discussion:
!
!    The only legal input name is 'hold_value', which indicates the
!    value of the variable to be held fixed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NAME, the name of the variable.
!
!    Output, real ( kind = 8 ) VALUE, the value of the variable.
!
  implicit none

  character ( len = *  ) name
  real      ( kind = 8 ) value

  call p20_r8_store ( 'get', name, value )

  return
end
subroutine p20_r8_set ( name, value )

!*****************************************************************************80
!
!! P20_R8_SET sets the value of a real parameter for problem 20.
!
!  Discussion:
!
!    The only legal input name is 'hold_value', which indicates the
!    value of the variable to be held fixed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NAME, the name of the variable.
!
!    Input, real ( kind = 8 ) VALUE, the value of the variable.
!
  implicit none

  character ( len = *  ) name
  real      ( kind = 8 ) value

  call p20_r8_store ( 'set', name, value )

  return
end
subroutine p20_r8_store ( action, name, value )

!*****************************************************************************80
!
!! P20_R8_STORE sets or gets the value of a real parameter for problem 20.
!
!  Discussion:
!
!    The only legal input name is 'hold_value', which indicates the
!    value of the variable to be held fixed.  This variable is given
!    the default value of pi/8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, either 'get' or 'set'.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!
!    Input/output, real ( kind = 8 ) VALUE, the value of the variable.
!
  implicit none

  character ( len = *  ) action
  real      ( kind = 8 ), save :: hold_value = 0.39269908169872415481
  character ( len = *  ) name
  logical                s_eqi
  real      ( kind = 8 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'hold_value' ) ) then
      value = hold_value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_R8_STORE - Fatal error!'
      write ( *, '(a)' ) '  Unexpected "get" of NAME = "' // name // '".'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'hold_value' ) ) then
      hold_value = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P20_R8_STORE - Fatal error!'
      write ( *, '(a)' ) '  Unexpected "set" of NAME = "' // name // '".'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P20_R8_STORE - Fatal error!'
    write ( *, '(a)' ) '  Unexpected ACTION = "' // action // '".'
    stop
  end if

  return
end
subroutine p20_setup ( l, theta, lambda, mu )

!*****************************************************************************80
!
!! P20_SETUP finds a solution (L,THETA,LAMBDA,MU) given L and THETA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) L, THETA, the values of L and THETA.
!
!    Output, real ( kind = 8 ) LAMBDA, MU, values for LAMBDA and MU.
!
  implicit none

  real ( kind = 8 ) l
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu
  real ( kind = 8 ) theta

  mu = 2.0D+00 * ( 1.0D+00 - l ) * sin ( theta ) &
    - 0.5D+00 * cos ( theta ) * theta / l

  lambda = ( ( 1.0D+00 - l ) - 0.5D+00 * mu * sin ( theta ) ) / cos ( theta )

  return
end
subroutine p20_start ( option, nvar, x )

!*****************************************************************************80
!
!! P20_START returns a starting point for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Input, integer ( kind = 4 ) NVAR, the number of variables.
!
!    Output, real ( kind = 8 ) X(NVAR), the starting point.
!
  implicit none

  integer ( kind = 4 ) nvar

  integer ( kind = 4 ) hold_index
  real ( kind = 8 ) hold_value
  real ( kind = 8 ) l
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(nvar)

  call p20_i4_get ( 'hold_index', hold_index )
  call p20_r8_get ( 'hold_value', hold_value )

  if ( hold_index == 1 ) then
    l = hold_value
    theta = pi / 8.0D+00
  else if ( hold_index == 2 ) then
    l = 0.25D+00
    theta = hold_value
  else
    l = 0.25D+00
    theta = pi / 8.0D+00
  end if

  call p20_setup ( l, theta, lambda, mu )

  x(1:nvar) = (/ l, theta, lambda, mu /)

  return
end
subroutine p20_stepsize ( option, h, hmin, hmax )

!*****************************************************************************80
!
!! P20_STEPSIZE returns step sizes for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, real ( kind = 8 ) H, HMIN, HMAX, suggested values for the
!    initial step, the minimum step, and the maximum step.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) hmax
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) option

  h =    0.0025000D+00
  hmin = 0.01000D+00
  hmax = 0.08000D+00

  return
end
subroutine p20_title ( option, title )

!*****************************************************************************80
!
!! P20_TITLE sets the title for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the option index.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!    TITLE will never be longer than 80 characters.
!
  implicit none

  integer ( kind = 4 ) option
  character ( len = * ) title

  if ( option == 1 ) then
    title = 'The Buckling Spring, F(L,Theta,Lambda,Mu).'
  else
    title = '???'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P20_TITLE - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized option number.'
    stop
  end if

  return
end
function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8 value.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, real ( kind = 8 ) R8_MOP, the I-th power of -1.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_mop

  if ( mod ( i, 2 ) == 0 ) then
    r8_mop = + 1.0D+00
  else
    r8_mop = - 1.0D+00
  end if

  return
end
function r8_sign ( x )

!*****************************************************************************80
!
!! R8_SIGN returns the sign of an R8.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value =  0 if X => 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose sign is desired.
!
!    Output, real ( kind = 8 ) R8_SIGN, the sign of X:
!
  implicit none

  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    r8_sign = -1.0D+00
  else
    r8_sign = +1.0D+00
  end if

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8mat_det ( n, a, det )

!*****************************************************************************80
!
!! R8MAT_DET computes the determinant of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Helmut Spaeth
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 125-127.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) piv(1)
  real ( kind = 8 ) t

  b(1:n,1:n) = a(1:n,1:n)

  det = 1.0D+00

  do k = 1, n

    piv = maxloc ( abs ( b(k:n,k) ) )

    m = piv(1) + k - 1

    if ( m /= k ) then
      det = - det
      t      = b(m,k)
      b(m,k) = b(k,k)
      b(k,k) = t
    end if

    det = det * b(k,k)

    if ( b(k,k) /= 0.0D+00 ) then

      b(k+1:n,k) = -b(k+1:n,k) / b(k,k)

      do j = k + 1, n
        if ( m /= k ) then
          t      = b(m,j)
          b(m,j) = b(k,j)
          b(k,j) = t
        end if
        b(k+1:n,j) = b(k+1:n,j) + b(k+1:n,k) * b(k,j)
      end do

    end if

  end do

  return
end
subroutine r8mat_nullspace ( m, n, a, nullspace_size, nullspace )

!*****************************************************************************80
!
!! R8MAT_NULLSPACE computes the nullspace of a matrix.
!
!  Discussion:
!
!    Let A be an MxN matrix.
!
!    If X is an N-vector, and A*X = 0, then X is a null vector of A.
!
!    The set of all null vectors of A is called the nullspace of A.
!
!    The 0 vector is always in the null space.
!
!    If the 0 vector is the only vector in the nullspace of A, then A
!    is said to have maximum column rank.  (Because A*X=0 can be regarded
!    as a linear combination of the columns of A).  In particular, if A
!    is square, and has maximum column rank, it is nonsingular.
!
!    The dimension of the nullspace is the number of linearly independent
!    vectors that span the nullspace.  If A has maximum column rank,
!    its nullspace has dimension 0.
!
!    This routine uses the reduced row echelon form of A to determine
!    a set of NULLSPACE_SIZE independent null vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix to be analyzed.
!
!    Input, integer ( kind = 4 ) NULLSPACE_SIZE, the size of the nullspace.
!
!    Output, real ( kind = 8 ) NULLSPACE(N,NULLSPACE_SIZE), vectors that
!    span the nullspace.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nullspace_size

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) col(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  real ( kind = 8 ) nullspace(n,nullspace_size)
  integer ( kind = 4 ) row(m)
  real ( kind = 8 ) rref(m,n)
!
!  Make a copy of A.
!
  rref(1:m,1:n) = a(1:m,1:n)
!
!  Get the reduced row echelon form of A.
!
  call r8mat_rref ( m, n, rref )
!
!  Note in ROW the columns of the leading nonzeros.
!  COL(J) = +J if there is a leading 1 in that column, and -J otherwise.
!
  row(1:m) = 0

  do j = 1, n
    col(j) = - j
  end do

  do i = 1, m
    do j = 1, n
      if ( rref(i,j) == 1.0D+00 ) then
        row(i) = j
        col(j) = j
        exit
      end if
    end do
  end do

  nullspace(1:n,1:nullspace_size) = 0.0D+00

  j2 = 0
!
!  If column J does not contain a leading 1, then it contains
!  information about a null vector.
!
  do j = 1, n

    if ( col(j) < 0 ) then

      j2 = j2 + 1

      do i = 1, m
        if ( rref(i,j) /= 0.0D+00 ) then
          i2 = row(i)
          nullspace(i2,j2) = - rref(i,j)
        end if
      end do

      nullspace(j,j2) = 1.0D+00

    end if

  end do

  return
end
subroutine r8mat_nullspace_size ( m, n, a, nullspace_size )

!*****************************************************************************80
!
!! R8MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.
!
!  Discussion:
!
!    Let A be an MxN matrix.
!
!    If X is an N-vector, and A*X = 0, then X is a null vector of A.
!
!    The set of all null vectors of A is called the nullspace of A.
!
!    The 0 vector is always in the null space.
!
!    If the 0 vector is the only vector in the nullspace of A, then A
!    is said to have maximum column rank.  (Because A*X=0 can be regarded
!    as a linear combination of the columns of A).  In particular, if A
!    is square, and has maximum column rank, it is nonsingular.
!
!    The dimension of the nullspace is the number of linearly independent
!    vectors that span the nullspace.  If A has maximum column rank,
!    its nullspace has dimension 0.
!
!    This routine ESTIMATES the dimension of the nullspace.  Cases of
!    singularity that depend on exact arithmetic will probably be missed.
!
!    The nullspace will be estimated by counting the leading 1's in the
!    reduced row echelon form of A, and subtracting this from N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix to be analyzed.
!
!    Output, integer ( kind = 4 ) NULLSPACE_SIZE, the estimated size
!    of the nullspace.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leading
  integer ( kind = 4 ) nullspace_size
  real ( kind = 8 ) rref(m,n)
!
!  Get the reduced row echelon form of A.
!
  rref(1:m,1:n) = a(1:m,1:n)

  call r8mat_rref ( m, n, rref )
!
!  Count the leading 1's in A.
!
  leading = 0
  do i = 1, m
    do j = 1, n
      if ( rref(i,j) == 1.0D+00 ) then
        leading = leading + 1
        exit
      end if
    end do
  end do

  nullspace_size = n - leading

  return
end
subroutine r8mat_rref ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_RREF computes the reduced row echelon form of a matrix.
!
!  Discussion:
!
!    A matrix is in row echelon form if:
!
!    * The first nonzero entry in each row is 1.
!
!    * The leading 1 in a given row occurs in a column to
!      the right of the leading 1 in the previous row.
!
!    * Rows which are entirely zero must occur last.
!
!    The matrix is in reduced row echelon form if, in addition to
!    the first three conditions, it also satisfies:
!
!    * Each column containing a leading 1 has no other nonzero entries.
!
!  Example:
!
!    Input matrix:
!
!     1.0  3.0  0.0  2.0  6.0  3.0  1.0
!    -2.0 -6.0  0.0 -2.0 -8.0  3.0  1.0
!     3.0  9.0  0.0  0.0  6.0  6.0  2.0
!    -1.0 -3.0  0.0  1.0  0.0  9.0  3.0
!
!    Output matrix:
!
!     1.0  3.0  0.0  0.0  2.0  0.0  0.0
!     0.0  0.0  0.0  1.0  2.0  0.0  0.0
!     0.0  0.0  0.0  0.0  0.0  1.0  0.3
!     0.0  0.0  0.0  0.0  0.0  0.0  0.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix A.
!
!    Input/output, real ( kind = 8 ) A(M,N).  On input, the matrix to be
!    analyzed.  On output, the RREF form of the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lead
  integer ( kind = 4 ) r
  real ( kind = 8 ) temp

  lead = 1

  do r = 1, m

    if ( n < lead ) then
      exit
    end if

    i = r

    do while ( a(i,lead) == 0.0 )

      i = i + 1

      if ( m < i ) then
        i = r
        lead = lead + 1
        if ( n < lead ) then
          lead = -1
          exit
        end if
      end if

    end do

    if ( lead < 0 ) then
      exit
    end if

    do j = 1, n
      temp   = a(i,j)
      a(i,j) = a(r,j)
      a(r,j) = temp
    end do

    a(r,1:n) = a(r,1:n) / a(r,lead)

    do i = 1, m
      if ( i /= r ) then
        a(i,1:n) = a(i,1:n) - a(i,lead) * a(r,1:n)
      end if
    end do

    lead = lead + 1

  end do

  return
end
subroutine r8vec_amax_index ( n, a, amax_index )

!*****************************************************************************80
!
!! R8VEC_AMAX_INDEX returns the index of the maximum absolute value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMAX_INDEX, the index of the entry of
!    largest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amax
  integer ( kind = 4 ) amax_index
  integer ( kind = 4 ) i

  if ( n <= 0 ) then

    amax_index = -1

  else

    amax_index = 1
    amax = abs ( a(1) )

    do i = 2, n
      if ( amax < abs ( a(i) ) ) then
        amax_index = i
        amax = abs ( a(i) )
      end if
    end do

  end if

  return
end
function r8vec_norm_l2 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_L2, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm_l2

  r8vec_norm_l2 = sqrt ( sum ( a(1:n)**2 ) )

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Discussion:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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

  character              c1
  character              c2
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lenc
  logical                s_eqi
  character ( len = *  ) s1
  integer   ( kind = 4 ) s1_length
  character ( len = *  ) s2
  integer   ( kind = 4 ) s2_length

  s1_length = len ( s1 )
  s2_length = len ( s2 )
  lenc = min ( s1_length, s2_length )

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

  do i = lenc + 1, s1_length
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, s2_length
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine sge_check ( lda, m, n, ierror )

!*****************************************************************************80
!
!! SGE_CHECK checks the dimensions of a general matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if LDA is illegal;
!    IERROR = IERROR + 2 if M is illegal;
!    IERROR = IERROR + 4 if N is illegal.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  ierror = 0

  if ( lda < m ) then
    ierror = ierror + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_CHECK - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal LDA = ', lda
    stop
  end if

  if ( m < 1 ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_CHECK - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal M = ', m
    stop
  end if

  if ( n < 1 ) then
    ierror = ierror + 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_CHECK - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal N = ', n
    stop
  end if

  return
end
subroutine sge_fa ( lda, n, a, pivot, info )

!*****************************************************************************80
!
!! SGE_FA factors a general matrix.
!
!  Discussion:
!
!    SGE_FA is a simplified version of the LINPACK routine SGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the matrix to be factored.
!    On output, A contains an upper triangular matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = L * U, where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Output, integer ( kind = 4 ) PIVOT(N), a vector of pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_FA - Fatal error!'
    write ( *, '(a)' ) '  Illegal dimensions.'
    stop
  end if

  info = 0

  do k = 1, n - 1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k + 1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGE_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      call r8_swap ( a(l,k), a(k,k) )
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = -a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k + 1, n

      if ( l /= k ) then
        call r8_swap ( a(l,j), a(k,j) )
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  pivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
subroutine sge_sl ( lda, n, a, pivot, b, job )

!*****************************************************************************80
!
!! SGE_SL solves a system factored by SGE_FA.
!
!  Discussion:
!
!    SGE_SL is a simplified version of the LINPACK routine SGESL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(LDA,N), the LU factors from SGE_FA.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from SGE_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB, specifies the operation.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_SL - Fatal error!'
    write ( *, '(a)' ) '  Illegal dimensions.'
    stop
  end if
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve PL * Y = B.
!
    do k = 1, n - 1

      l = pivot(k)

      if ( l /= k ) then
        call r8_swap ( b(l), b(k) )
      end if

      b(k+1:n) = b(k+1:n) + a(k+1:n,k) * b(k)

    end do
!
!  Solve U * X = Y.
!
    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      b(1:k-1) = b(1:k-1) - a(1:k-1,k) * b(k)
    end do
!
!  Solve transpose ( A ) * X = B.
!
  else
!
!  Solve transpose ( U ) * Y = B.
!
    do k = 1, n
      b(k) = ( b(k) - dot_product ( b(1:k-1), a(1:k-1,k) ) ) / a(k,k)
    end do
!
!  Solve transpose ( PL ) * X = Y.
!
    do k = n - 1, 1, -1

      b(k) = b(k) + dot_product ( b(k+1:n), a(k+1:n,k) )

      l = pivot(k)

      if ( l /= k ) then
        call r8_swap ( b(l), b(k) )
      end if

    end do

  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine zero_rc ( a, b, t, arg, status, value )

!*****************************************************************************80
!
!! ZERO_RC seeks the root of a function F(X) using reverse communication.
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    The routine is a revised version of the Brent zero finder
!    algorithm, using reverse communication.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the change of sign interval.
!
!    Input, real ( kind = 8 ) T, a positive error tolerance.
!
!    Output, real ( kind = 8 ) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS zero, ARG is the routine's
!    estimate for the function's zero.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between
!    the user and the routine.  The user only sets STATUS to zero on the first
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated zero
!
!    Input, real ( kind = 8 ) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ), save :: c
  real ( kind = 8 ), save :: d
  real ( kind = 8 ), save :: e
  real ( kind = 8 ), save :: fa
  real ( kind = 8 ), save :: fb
  real ( kind = 8 ), save :: fc
  real ( kind = 8 ) m
  real ( kind = 8 ), save :: machep
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ), save :: sa
  real ( kind = 8 ), save :: sb
  integer ( kind = 4 ) status
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) value
!
!  Input STATUS = 0.
!  Initialize, request F(A).
!
  if ( status == 0 ) then

    machep = epsilon ( a )

    sa = a
    sb = b
    e = sb - sa
    d = e

    status = 1
    arg = a
    return
!
!  Input STATUS = 1.
!  Receive F(A), request F(B).
!
  else if ( status == 1 ) then

    fa = value

    status = 2
    arg = sb
    return
!
!  Input STATUS = 2
!  Receive F(B).
!
  else if ( status == 2 ) then

    fb = value

    if ( 0.0D+00 < fa * fb ) then
      status = -1
      return
    end if

    c = sa
    fc = fa

  else

    fb = value

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

  end if
!
!  Compute the next point at which a function value is requested.
!
  if ( abs ( fc ) < abs ( fb ) ) then

    sa = sb
    sb = c
    c = sa
    fa = fb
    fb = fc
    fc = fa

  end if

  tol = 2.0D+00 * machep * abs ( sb ) + t
  m = 0.5D+00 * ( c - sb )

  if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
    status = 0
    arg = sb
    return
  end if

  if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

    e = m
    d = e

  else

    s = fb / fa

    if ( sa == c ) then

      p = 2.0D+00 * m * s
      q = 1.0D+00 - s

    else

      q = fa / fc
      r = fb / fc
      p = s * ( 2.0D+00 * m * a * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
      q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

    end if

    if ( 0.0D+00 < p ) then
      q = - q
    else
      p = - p
    end if

    s = e
    e = d

    if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
      p < abs ( 0.5D+00 * s * q ) ) then
      d = p / q
    else
      e = m
      d = e
    end if

  end if

  sa = sb
  fa = fb

  if ( tol < abs ( d ) ) then
    sb = sb + d
  else if ( 0.0D+00 < m ) then
    sb = sb + tol
  else
    sb = sb - tol
  end if

  arg = sb
  status = status + 1

  return
end
