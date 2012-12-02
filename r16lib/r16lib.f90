function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I  I4_LOG_10
!    -----  --------
!        0    0
!        1    0
!        2    0
!        9    0
!       10    1
!       11    1
!       99    1
!      100    2
!      101    2
!      999    2
!     1000    3
!     1001    3
!     9999    3
!    10000    4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 10
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_10, the integer part of the
!    logarithm base 10 of the absolute value of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) ten_pow

  if ( i == 0 ) then

    i4_log_10 = 0

  else

    i4_log_10 = 0
    ten_pow = 10

    i_abs = abs ( i )

    do while ( ten_pow <= i_abs )
      i4_log_10 = i4_log_10 + 1
      ten_pow = ten_pow * 10
    end do

  end if

  return
end
function r16_abs ( x )

!*****************************************************************************80
!
!! R16_ABS returns the absolute value of an R16.
!
!  Discussion:
!
!    An R16 is a real ( kind = 16 ) value.
!
!    FORTRAN90 supplies the ABS function, which should be used instead
!    of this function!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number whose absolute value is desired.
!
!    Output, real ( kind = 16 ) R16_ABS, the absolute value of X.
!
  implicit none

  real ( kind = 16 ) r16_abs
  real ( kind = 16 ) x

  if ( 0.0_16 <= x ) then
    r16_abs = x
  else
    r16_abs = - x
  end if

  return
end
function r16_add ( x, y )

!*****************************************************************************80
!
!! R16_ADD returns the sum of two R16's.
!
!  Discussion:
!
!    An R16 is a real ( kind = 16 ) value.
!
!    FORTRAN90 supplies the + operator, which should generally be used instead
!    of this function!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, Y, the numbers to be added.
!
!    Output, real ( kind = 16 ) R16_ADD, the sum.
!
  implicit none

  real ( kind = 16 ) r16_add
  real ( kind = 16 ) x
  real ( kind = 16 ) y

  r16_add = x + y

  return
end
function r16_atan ( y, x )

!*****************************************************************************80
!
!! R16_ATAN computes the inverse tangent of the ratio Y / X.
!
!  Discussion:
!
!    R16_ATAN returns an angle whose tangent is ( Y / X ), a job which
!    the built in functions ATAN and ATAN2 already do.
!
!    However:
!
!    * R16_ATAN always returns a positive angle, between 0 and 2 PI,
!      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
!      and [-PI,+PI] respectively;
!
!    * R16_ATAN accounts for the signs of X and Y, (as does ATAN2).  The ATAN
!     function by contrast always returns an angle in the first or fourth
!     quadrants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) Y, X, two quantities which represent the
!    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, real ( kind = 16 ) R16_ATAN, an angle between 0 and 2 * PI, whose
!    tangent is (Y/X), and which lies in the appropriate quadrant so that
!    the signs of its cosine and sine match those of X and Y.
!
  implicit none

  real ( kind = 16 ) abs_x
  real ( kind = 16 ) abs_y
  real ( kind = 16 ) r16_atan
  real ( kind = 16 ), parameter :: r16_pi &
    = 3.1415926535897932384626433832795028_16
  real ( kind = 16 ) theta
  real ( kind = 16 ) theta_0
  real ( kind = 16 ) x
  real ( kind = 16 ) y
!
!  Special cases:
!
  if ( x == 0.0_16 ) then

    if ( 0.0_16 < y ) then
      theta = r16_pi / 2.0_16
    else if ( y < 0.0_16 ) then
      theta = 3.0_16 * r16_pi / 2.0_16
    else if ( y == 0.0_16 ) then
      theta = 0.0_16
    end if

  else if ( y == 0.0_16 ) then

    if ( 0.0_16 < x ) then
      theta = 0.0_16
    else if ( x < 0.0_16 ) then
      theta = r16_pi
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( 0.0_16 < x .and. 0.0_16 < y ) then
      theta = theta_0
    else if ( x < 0.0_16 .and. 0.0_16 < y ) then
      theta = r16_pi - theta_0
    else if ( x < 0.0_16 .and. y < 0.0_16 ) then
      theta = r16_pi + theta_0
    else if ( 0.0_16 < x .and. y < 0.0_16 ) then
      theta = 2.0_16 * r16_pi - theta_0
    end if

  end if

  r16_atan = theta

  return
end
function r16_cas ( x )

!*****************************************************************************80
!
!! R16_CAS returns the "casine" of an R16.
!
!  Discussion:
!
!    The "casine", used in the discrete Hartley transform, is abbreviated
!    CAS(X), and defined by:
!
!      CAS(X) = cos ( X ) + sin( X )
!             = sqrt ( 2 ) * sin ( X + pi/4 )
!             = sqrt ( 2 ) * cos ( X - pi/4 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ralph Hartley,
!    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
!    Proceedings of the Institute of Radio Engineers,
!    Volume 30, pages 144-150, 1942.
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number whose casine is desired.
!
!    Output, real ( kind = 16 ) R16_CAS, the casine of X, which will be between
!    plus or minus the square root of 2.
!
  implicit none

  real ( kind = 16 ) r16_cas
  real ( kind = 16 ) x

  r16_cas = cos ( x ) + sin ( x )

  return
end
function r16_ceiling ( r )

!*****************************************************************************80
!
!! R16_CEILING rounds an R16 "up" (towards +oo) to the next I4.
!
!  Example:
!
!    R     Value
!
!   -1.1  -1
!   -1.0  -1
!   -0.9   0
!    0.0   0
!    5.0   5
!    5.1   6
!    5.9   6
!    6.0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) R, the value to be rounded up.
!
!    Output, integer ( kind = 4 ) R16_CEILING, the rounded value.
!
  implicit none

  real    ( kind = 16 ) r
  integer ( kind = 4 ) r16_ceiling
  integer ( kind = 4 ) value

  value = int ( r, kind = 4 )

  if ( real ( value, kind = 16 ) < r ) then
    value = value + 1
  end if

  r16_ceiling = value

  return
end
function r16_choose ( n, k )

!*****************************************************************************80
!
!! R16_CHOOSE computes the binomial coefficient C(N,K) as an R16.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in R16 arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, real ( kind = 16 ) R16_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  real ( kind = 16 ) r16_choose
  real ( kind = 16 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0.0_16

  else if ( mn == 0 ) then

    value = 1.0_16

  else

    mx = max ( k, n - k )
    value = real ( mx + 1, kind = 16 )

    do i = 2, mn
      value = ( value * real ( mx + i, kind = 16 ) ) / real ( i, kind = 16 )
    end do

  end if

  r16_choose = value

  return
end
function r16_chop ( place, x )

!*****************************************************************************80
!
!! R16_CHOP chops an R16 to a given number of binary places.
!
!  Example:
!
!    3.875 = 2 + 1 + 1/2 + 1/4 + 1/8.
!
!    The following values would be returned for the 'chopped' value of
!    3.875:
!
!    PLACE  Value
!
!       1      2
!       2      3     = 2 + 1
!       3      3.5   = 2 + 1 + 1/2
!       4      3.75  = 2 + 1 + 1/2 + 1/4
!       5+     3.875 = 2 + 1 + 1/2 + 1/4 + 1/8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PLACE, the number of binary places to preserve.
!    PLACE = 0 means return the integer part of X.
!    PLACE = 1 means return the value of X, correct to 1/2.
!    PLACE = 2 means return the value of X, correct to 1/4.
!    PLACE = -1 means return the value of X, correct to 2.
!
!    Input, real ( kind = 16 ) X, the number to be chopped.
!
!    Output, real ( kind = 16 ) R16_CHOP, the chopped number.
!
  implicit none

  real ( kind = 16 ) fac
  integer ( kind = 4 ) place
  real ( kind = 16 ) r16_chop
  real ( kind = 16 ) r16_log_2
  real ( kind = 16 ) r16_sign
  real ( kind = 16 ) s
  integer ( kind = 4 ) temp
  real ( kind = 16 ) x

  s = r16_sign ( x )
  temp = int ( r16_log_2 ( abs ( x ) ) )
  fac = 2.0_16 ** ( temp - place + 1 )
  r16_chop = s * real ( int ( abs ( x ) / fac ), kind = 16 ) * fac

  return
end
function r16_csqrt ( x )

!*****************************************************************************80
!
!! R16_CSQRT returns the complex square root of an R16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number whose square root is desired.
!
!    Output, complex ( kind = 16 ) R16_CSQRT, the square root of X:
!
  implicit none

  real ( kind = 16 ) argument
  real ( kind = 16 ) magnitude
  real ( kind = 16 ), parameter :: r16_pi &
    = 3.1415926535897932384626433832795028_16
  complex ( kind = 16 ) r16_csqrt
  real ( kind = 16 ) x

  if ( 0.0_16 < x ) then
    magnitude = x
    argument = 0.0_16
  else if ( 0.0_16 == x ) then
    magnitude = 0.0_16
    argument = 0.0_16
  else if ( x < 0.0_16 ) then
    magnitude = - x
    argument = r16_pi
  end if

  magnitude = sqrt ( magnitude )
  argument = argument / 2.0_16

  r16_csqrt = magnitude * cmplx ( cos ( argument ), &
                                  sin ( argument ), kind = 16 )

  return
end
function r16_cube_root ( x )

!*****************************************************************************80
!
!! R16_CUBE_ROOT returns the cube root of an R16.
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
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number whose cube root is desired.
!
!    Output, real ( kind = 16 ) R16_CUBE_ROOT, the cube root of X.
!
  implicit none

  real ( kind = 16 ) r16_cube_root
  real ( kind = 16 ) value
  real ( kind = 16 ) x

  if ( 0.0_16 < x ) then
    value = x ** ( 1.0_16 / 3.0_16 )
  else if ( x == 0.0_16 ) then
    value = 0.0_16
  else
    value = -( abs ( x ) ) ** ( 1.0_16 / 3.0_16 )
  end if

  r16_cube_root = value

  return
end
function r16_diff ( x, y, n )

!*****************************************************************************80
!
!! R16_DIFF computes the difference of two R16's to a specified accuracy.
!
!  Discussion:
!
!    The user controls how many binary digits of accuracy
!    are to be used.
!
!    N determines the accuracy of the value of the result.  If N = 10,
!    for example, only 11 binary places will be used in the arithmetic.
!    In general, only N+1 binary places will be used.
!
!    N may be zero.  However, a negative value of N should
!    not be used, since this will cause both X and Y to look like 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, Y, the two values whose difference is desired.
!
!    Input, integer ( kind = 4 ) N, the number of binary digits to use.
!
!    Output, real ( kind = 16 ) R16_DIFF, the value of X-Y.
!
  implicit none

  real ( kind = 16 ) cx
  real ( kind = 16 ) cy
  integer ( kind = 4 ) n
  real ( kind = 16 ) pow2
  real ( kind = 16 ) r16_diff
  real ( kind = 16 ) size
  real ( kind = 16 ) x
  real ( kind = 16 ) y

  if ( x == y ) then
    r16_diff = 0.0_16
    return
  end if

  pow2 = 2.0_16 ** n
!
!  Compute the magnitude of X and Y, and take the larger of the
!  two.  At least one of the two values is not zero!
!
  size = max ( abs ( x ), abs ( y ) )
!
!  Make normalized copies of X and Y.  One of the two values will
!  actually be equal to 1.
!
  cx = x / size
  cy = y / size
!
!  Here's where rounding comes in.  We know that the larger of the
!  the two values equals 1.  We multiply both values by 2**N,
!  where N+1 is the number of binary digits of accuracy we want
!  to use, truncate the values, and divide back by 2**N.
!
  cx = real ( int ( cx * pow2 + sign ( 0.5_16, cx ) ), kind = 16 ) / pow2
  cy = real ( int ( cy * pow2 + sign ( 0.5_16, cy ) ), kind = 16 ) / pow2
!
!  Take the difference now.
!
  r16_diff = cx - cy
!
!  Undo the scaling.
!
  r16_diff = r16_diff * size

  return
end
subroutine r16_digit ( x, idigit, digit )

!*****************************************************************************80
!
!! R16_DIGIT returns a particular decimal digit of an R16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number whose NDIG-th decimal digit
!    is desired.  If X is zero, all digits will be returned as 0.
!
!    Input, integer ( kind = 4 ) IDIGIT, the position of the desired decimal
!    digit.  A value of 1 means the leading digit, a value of 2 the second digit
!    and so on.
!
!    Output, integer ( kind = 4 ) DIGIT, the value of the IDIGIT-th decimal
!    digit of X.
!
  implicit none

  integer ( kind = 4 ) digit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idigit
  integer ( kind = 4 ) ival
  real ( kind = 16 ) x
  real ( kind = 16 ) xcopy

  if ( x == 0.0_16 ) then
    digit = 0
    return
  end if

  if ( idigit <= 0 ) then
    digit = 0
    return
  end if
!
!  Set XCOPY = X, and then force XCOPY to lie between 1 and 10.
!
  xcopy = abs ( x )

  do while ( xcopy < 1.0_16 )
    xcopy = xcopy * 10.0_16
  end do

  do while ( 10.0_16 <= xcopy )
    xcopy = xcopy / 10.0_16
  end do

  do i = 1, idigit
    ival = int ( xcopy )
    xcopy = ( xcopy - ival ) * 10.0_16
  end do

  digit = ival

  return
end
function r16_epsilon ( )

!*****************************************************************************80
!
!! R16_EPSILON returns the R16 roundoff unit.
!
!  Discussion:
!
!    An R16 is a real ( kind = 16 ) value.
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 16 ) R16_EPSILON, the round-off unit.
!
  implicit none

  real ( kind = 16 ) one
  real ( kind = 16 ) r16_add
  real ( kind = 16 ) r16_epsilon
  real ( kind = 16 ) r16_sub
  real ( kind = 16 ) temp
  real ( kind = 16 ) test
  real ( kind = 16 ) value

  one = real ( 1, kind = 16 )

  value = one
  temp = value / 2.0_16
  test = r16_add ( one, temp )

  do while ( one < test )
    value = temp
    temp = value / 2.0_16
    test = r16_add ( one, temp )
  end do

  r16_epsilon = value

  return
end
function r16_exp ( x )

!*****************************************************************************80
!
!! R16_EXP computes the exponential of an R16, avoiding overflow and underflow.
!
!  Discussion:
!
!    My experience with the G95 compiler has included many unpleasant
!    floating point exceptions when very small arguments are given to
!    the exponential function.
!
!    This routine is designed to avoid such problems.
!
!    Ideally, the rule would be:
!
!                    X <= log ( TINY ) => R16_EXP ( X ) = 0
!    log ( HUGE ) <= X                 => R16_EXP ( X ) = HUGE
!
!    However, the G95 math library seems to produce infinity for
!    EXP ( LOG ( HUGE ( X ) ), rather than HUGE ( X ), so we've
!    included a fudge factor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the argument of the exponential function.
!
!    Output, real ( kind = 16 ) R16_EXP, the value of exp ( X ).
!
  implicit none

  real ( kind = 16 ), parameter :: log_max = 709.711_16
  real ( kind = 16 ), parameter :: log_min = -708.467_16
  real ( kind = 16 ) r16_exp
  real ( kind = 16 ) x

  if ( x <= log_min ) then
    r16_exp = 0.0_16
  else if ( x < log_max ) then
    r16_exp = exp ( x )
  else
    r16_exp = huge ( x )
  end if

  return
end
function r16_factorial ( n )

!*****************************************************************************80
!
!! R16_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = 16 ) R16_FACTORIAL, the factorial of N.
!
  implicit none

  real ( kind = 16 ) r16_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  r16_factorial = 1.0_16

  do i = 1, n
    r16_factorial = r16_factorial * real ( i, kind = 16 )
  end do

  return
end
function r16_factorial2 ( n )

!*****************************************************************************80
!
!! R16_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    N!!
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the double factorial
!    function.  If N is less than 1,the value is returned as 1.0.
!
!    Output, real ( kind = 16 ) R16_FACTORIAL2, the value.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 16 ) r16_factorial2
  real ( kind = 16 ) r16_n

  if ( n < 1 ) then
    r16_factorial2 = 1.0_16
    return
  end if

  r16_n = real ( n, kind = 16 )
  r16_factorial2 = 1.0_16

  do while ( 1.0_16 < r16_n )
    r16_factorial2 = r16_factorial2 * r16_n
    r16_n = r16_n - 2.0_16
  end do

  return
end
function r16_floor ( r )

!*****************************************************************************80
!
!! R16_FLOOR rounds an R16 "down" (towards -oo) to the next integer.
!
!  Example:
!
!    R     Value
!
!   -1.1  -2
!   -1.0  -1
!   -0.9  -1
!    0.0   0
!    5.0   5
!    5.1   5
!    5.9   5
!    6.0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) R, the value to be rounded down.
!
!    Output, integer ( kind = 4 ) R16_FLOOR, the rounded value.
!
  implicit none

  real ( kind = 16 ) r
  integer ( kind = 4 ) r16_floor
  integer ( kind = 4 ) value

  value = int ( r )
  if ( r < real ( value, kind = 16 ) ) then
    value = value - 1
  end if

  r16_floor = value

  return
end
function r16_fraction ( x )

!*****************************************************************************80
!
!! R16_FRACTION returns the fraction part of an R16.
!
!  Discussion:
!
!    If we regard a real number as
!
!      R = SIGN * ( WHOLE + FRACTION )
!
!    where
!
!      SIGN is +1 or -1,
!      WHOLE is a nonnegative integer
!      FRACTION is a nonnegative real number strictly less than 1,
!
!    then this routine returns the value of FRACTION.
!
!  Example:
!
!     R      FRACTION
!
!    0.00      0.00
!    1.01      0.01
!    2.02      0.02
!   19.73      0.73
!   -4.34      0.34
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the argument.
!
!    Output, real ( kind = 16 ) R16_FRACTION, the fraction part of X.
!
  implicit none

  real ( kind = 16 ) r16_fraction
  real ( kind = 16 ) x

  r16_fraction = abs ( x ) - real ( int ( abs ( x ) ), kind = 16 )

  return
end
function r16_huge ( )

!*****************************************************************************80
!
!! R16_HUGE returns a very large R16.
!
!  Discussion:
!
!    An R16 is a real ( kind = 16 ) value.
!
!    The value returned by this function is NOT required to be the
!    maximum representable R16.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 16 ) R16_HUGE, a "huge" value.
!
  implicit none

  real ( kind = 16 ) r16_huge

  r16_huge = 1.0E+30_16

  return
end
function r16_in_01 ( a )

!*****************************************************************************80
!
!! R16_IN_01 is TRUE if an R16 is in the range [0,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A, the value.
!
!    Output, logical R16_IN_01, is TRUE if 0 <= A <= 1.
!
  implicit none

  real ( kind = 16 ) a
  logical r16_in_01
  logical value

  if ( a < 0.0_16 .or. 1.0_16 < a ) then
    value = .false.
  else
    value = .true.
  end if

  r16_in_01 = value

  return
end
function r16_is_int ( r )

!*****************************************************************************80
!
!! R16_IS_INT determines if an R16 represents an integer value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) R, the number to be checked.
!
!    Output, logical R16_IS_INT, is TRUE if R is an integer value.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  real ( kind = 16 ) r
  logical r16_is_int
  logical value

  if ( real ( i4_huge, kind = 16 ) < r ) then
    value = .false.
  else if ( r < - real ( i4_huge, kind = 16 ) ) then
    value = .false.
  else if ( r == real ( int ( r ), kind = 16 ) ) then
    value = .true.
  else
    value = .false.
  end if

  r16_is_int = value

  return
end
function r16_log_2 ( x )

!*****************************************************************************80
!
!! R16_LOG_2 returns the logarithm base 2 of an R16.
!
!  Discussion:
!
!    value = Log ( |X| ) / Log ( 2.0 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number whose base 2 logarithm is desired.
!    X should not be 0.
!
!    Output, real ( kind = 16 ) R16_LOG_2, the logarithm base 2 of the absolute
!    value of X.  It should be true that |X| = 2**R16_LOG_2.
!
  implicit none

  real ( kind = 16 ) r16_log_2
  real ( kind = 16 ) x

  if ( x == 0.0_16 ) then
    r16_log_2 = - huge ( x )
  else
    r16_log_2 = log ( abs ( x ) ) / log ( 2.0_16 )
  end if

  return
end
function r16_log_10 ( x )

!*****************************************************************************80
!
!! R16_LOG_10 returns the logarithm base 10 of an R16.
!
!  Discussion:
!
!    value = Log10 ( |X| )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number whose base 2 logarithm is desired.
!    X should not be 0.
!
!    Output, real ( kind = 16 ) R16_LOG_10, the logarithm base 10 of the
!    absolute value of X.  It should be true that |X| = 10**R16_LOG_10.
!
  implicit none

  real ( kind = 16 ) r16_log_10
  real ( kind = 16 ) x

  if ( x == 0.0_16 ) then
    r16_log_10 = - huge ( x )
  else
    r16_log_10 = log10 ( abs ( x ) )
  end if

  return
end
function r16_log_b ( x, b )

!*****************************************************************************80
!
!! R16_LOG_B returns the logarithm base B of an R16.
!
!  Discussion:
!
!    value = log ( |X| ) / log ( |B| )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number whose base B logarithm is desired.
!    X should not be 0.
!
!    Input, real ( kind = 16 ) B, the base, which should not be 0, 1 or -1.
!
!    Output, real ( kind = 16 ) R16_LOG_B, the logarithm base B of the absolute
!    value of X.  It should be true that |X| = |B|**R16_LOG_B.
!
  implicit none

  real ( kind = 16 ) b
  real ( kind = 16 ) r16_log_b
  real ( kind = 16 ) x

  if ( b == 0.0_16 .or. b == 1.0_16 .or. b == - 1.0_16 ) then
    r16_log_b = - huge ( x )
  else if ( abs ( x ) == 0.0_16 ) then
    r16_log_b = - huge ( x )
  else
    r16_log_b = log ( abs ( x ) ) / log ( abs ( b ) )
  end if

  return
end
subroutine r16_mant ( x, s, r, l )

!*****************************************************************************80
!
!! R16_MANT computes the "mantissa" or "fraction part" of an R16.
!
!  Discussion:
!
!    X = S * R * 2**L
!
!    S is +1 or -1,
!    R is an real value between 1.0 and 2.0,
!    L is an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number to be decomposed.
!
!    Output, integer ( kind = 4 ) S, the "sign" of the number.
!    S will be -1 if X is less than 0, and +1 if X is greater
!    than or equal to zero.
!
!    Output, real ( kind = 16 ) R, the mantissa of X.  R will be greater
!    than or equal to 1, and strictly less than 2.  The one
!    exception occurs if X is zero, in which case R will also
!    be zero.
!
!    Output, integer ( kind = 4 ) L, the integer part of the logarithm
!    (base 2) of X.
!
  implicit none

  integer ( kind = 4 ) l
  real ( kind = 16 ) r
  integer ( kind = 4 ) s
  real ( kind = 16 ) x
!
!  Determine the sign.
!
  if ( x < 0.0_16 ) then
    s = - 1
  else
    s = + 1
  end if
!
!  Set R to the absolute value of X, and L to zero.
!  Then force R to lie between 1 and 2.
!
  if ( x < 0.0_16 ) then
    r = - x
  else
    r = + x
  end if

  l = 0
!
!  Time to bail out if X is zero.
!
  if ( x == 0.0_16 ) then
    return
  end if

  do while ( 2.0_16 <= r )
    r = r / 2.0_16
    l = l + 1
  end do

  do while ( r < 1.0_16 )
    r = r * 2.0_16
    l = l - 1
  end do

  return
end
function r16_mod ( x, y )

!*****************************************************************************80
!
!! R16_MOD returns the remainder of R16 division.
!
!  Discussion:
!
!    If
!      REM = R16_MOD ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM has the same sign as X, and abs ( REM ) < Y.
!
!  Example:
!
!        X         Y     R16_MOD  R16_MOD Factorization
!
!      107        50       7      107 =  2 *  50 + 7
!      107       -50       7      107 = -2 * -50 + 7
!     -107        50      -7     -107 = -2 *  50 - 7
!     -107       -50      -7     -107 =  2 * -50 - 7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number to be divided.
!
!    Input, real ( kind = 16 ) Y, the number that divides X.
!
!    Output, real ( kind = 16 ) R16_MOD, the remainder when X is divided by Y.
!
  implicit none

  real ( kind = 16 ) r16_mod
  real ( kind = 16 ) x
  real ( kind = 16 ) y

  if ( y == 0.0_16 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R16_MOD - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R16_MOD ( X, Y ) called with Y = ', y
    stop
  end if

  r16_mod = x - real ( int ( x / y ), kind = 16 ) * y

  if ( x < 0.0_16 .and. 0.0_16 < r16_mod ) then
    r16_mod = r16_mod - abs ( y )
  else if ( 0.0_16 < x .and. r16_mod < 0.0_16 ) then
    r16_mod = r16_mod + abs ( y )
  end if

  return
end
function r16_modp ( x, y )

!*****************************************************************************80
!
!! R16_MODP returns the nonnegative remainder of R16 division.
!
!  Discussion:
!
!    An R16 is a real ( kind = 16 ) value.
!
!    If
!      REM = R16_MODP ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, R16_MODP(A,360.0) is between 0 and 360, always.
!
!  Example:
!
!        X         Y     MOD R16_MODP R16_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number to be divided.
!
!    Input, real ( kind = 16 ) Y, the number that divides X.
!
!    Output, real ( kind = 16 ) R16_MODP, the nonnegative remainder
!    when X is divided by Y.
!
  implicit none

  real ( kind = 16 ) r16_modp
  real ( kind = 16 ) x
  real ( kind = 16 ) y

  if ( y == 0.0_16 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R16_MODP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R16_MODP ( X, Y ) called with Y = ', y
    stop
  end if

  r16_modp = mod ( x, y )

  if ( r16_modp < 0.0_16 ) then
    r16_modp = r16_modp + abs ( y )
  end if

  return
end
function r16_mop ( i )

!*****************************************************************************80
!
!! R16_MOP returns the I-th power of -1 as an R16 value.
!
!  Discussion:
!
!    An R16 is a real ( kind = 16 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, real ( kind = 16 ) R16_MOP, the I-th power of -1.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 16 ) r16_mop

  if ( mod ( i, 2 ) == 0 ) then
    r16_mop = + 1.0_16
  else
    r16_mop = - 1.0_16
  end if

  return
end
function r16_nint ( x )

!*****************************************************************************80
!
!! R16_NINT returns the nearest integer to an R16.
!
!  Example:
!
!        X        R16_NINT
!
!      1.3         1
!      1.4         1
!      1.5         1 or 2
!      1.6         2
!      0.0         0
!     -0.7        -1
!     -1.1        -1
!     -1.6        -2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the value.
!
!    Output, integer ( kind = 4 ) R16_NINT, the nearest integer to X.
!
  implicit none

  integer ( kind = 4 ) r16_nint
  integer ( kind = 4 ) s
  real ( kind = 16 ) x

  if ( x < 0.0_16 ) then
    s = - 1
  else
    s = + 1
  end if

  r16_nint = s * int ( abs ( x ) + 0.5_16 )

  return
end
function r16_normal ( a, b, seed )

!*****************************************************************************80
!
!! R16_NORMAL returns a scaled pseudonormal R16.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A, the mean of the PDF.
!
!    Input, real ( kind = 16 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 16 ) R16_NORMAL, a sample of the normal PDF.
!
  implicit none

  real ( kind = 16 ) a
  real ( kind = 16 ) b

  real ( kind = 16 ) r1
  real ( kind = 16 ) r2
  real ( kind = 16 ) r16_normal
  real ( kind = 16 ), parameter :: r16_pi &
    = 3.1415926535897932384626433832795028_16
  real ( kind = 16 ) r16_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed2 = 0
  integer ( kind = 4 ), save :: used = 0
  real ( kind = 16 ) x
  real ( kind = 16 ), save :: y = 0.0_16
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal and its corresponding seed.
!
  if ( mod ( used, 2 ) == 0 ) then

    r1 = r16_uniform_01 ( seed )

    if ( r1 == 0.0_16 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R16_NORMAL - Fatal error!'
      write ( *, '(a)' ) '  R16_UNIFORM_01 returned a value of 0.'
      stop
    end if

    seed2 = seed
    r2 = r16_uniform_01 ( seed2 )

    x = sqrt ( - 2.0_16 * log ( r1 ) ) * cos ( 2.0_16 * r16_pi * r2 )
    y = sqrt ( - 2.0_16 * log ( r1 ) ) * sin ( 2.0_16 * r16_pi * r2 )
!
!  On odd calls, return the second normal and its corresponding seed.
!
  else

    seed = seed2
    x = y

  end if

  used = used + 1

  r16_normal = a + b * x

  return
end
function r16_normal_01 ( seed )

!*****************************************************************************80
!
!! R16_NORMAL_01 returns a unit pseudonormal R16.
!
!  Discussion:
!
!    An R16 is a real ( kind = 16 ) value.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
!    This means that on every other call, essentially, the input value of
!    SEED is ignored, since the code saves the second normal random value.
!
!    If you didn't know this, you might be confused since, usually, the
!    output of a random number generator can be completely controlled by
!    the input value of the SEED.  If I were more careful, I could rewrite
!    this routine so that it would distinguish between cases where the input
!    value of SEED is the output value from the previous call (all is well)
!    and those cases where it is not (the user has decided to do something
!    new.  Restart the uniform random number sequence.)  But I'll leave
!    that for later.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 16 ) R16_NORMAL_01, a sample of the standard
!    normal PDF.
!
  implicit none

  real ( kind = 16 ), parameter :: pi = 3.141592653589793_16
  real ( kind = 16 ) r1
  real ( kind = 16 ) r2
  real ( kind = 16 ) r16_normal_01
  real ( kind = 16 ) r16_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed2 = 0
  integer ( kind = 4 ), save :: used = 0
  real ( kind = 16 ) x
  real ( kind = 16 ), save :: y = 0.0_16
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal and its corresponding seed.
!
  if ( mod ( used, 2 ) == 0 ) then

    r1 = r16_uniform_01 ( seed )

    if ( r1 == 0.0_16 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R16_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R16_UNIFORM_01 returned a value of 0.'
      stop
    end if

    seed2 = seed
    r2 = r16_uniform_01 ( seed2 )

    x = sqrt ( - 2.0_16 * log ( r1 ) ) * cos ( 2.0_16 * pi * r2 )
    y = sqrt ( - 2.0_16 * log ( r1 ) ) * sin ( 2.0_16 * pi * r2 )
!
!  On odd calls, return the second normal and its corresponding seed.
!
  else

    seed = seed2
    x = y

  end if

  used = used + 1

  r16_normal_01 = x

  return
end
function r16_pi ( )

!*****************************************************************************80
!
!! R16_PI returns the value of pi as an R16.
!
!  Discussion:
!
!    An R16 is a real ( kind = 16 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 16 ) R16_PI, the value of pi.
!
  implicit none

  real ( kind = 16 ) r16_pi

  r16_pi = 3.1415926535897932384626433832795028_16

  return
end
function r16_power ( r, p )

!*****************************************************************************80
!
!! R16_POWER computes the P-th power of an R16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) R, the base.
!
!    Input, integer ( kind = 4 ) P, the power, which may be negative.
!
!    Output, real ( kind = 16 ) R16_POWER, the value of the P-th power of R.
!
  implicit none

  integer ( kind = 4 ) p
  real ( kind = 16 ) r
  real ( kind = 16 ) r16_power
  real ( kind = 16 ) value
!
!  Special case.  R^0 = 1.
!
  if ( p == 0 ) then

    value = 1.0_16
!
!  Special case.  Positive powers of 0 are 0.
!  For negative powers of 0, we go ahead and compute R^P,
!  relying on the software to complain.
!
  else if ( r == 0.0_16 ) then

    if ( 0 < p ) then
      value = 0.0_16
    else
      value = r**p
    end if

  else if ( 1 <= p ) then
    value = r**p
  else
    value = 1.0_16 / r**(-p)
  end if

  r16_power = value

  return
end
subroutine r16_power_fast ( r, p, rp, mults )

!*****************************************************************************80
!
!! R16_POWER_FAST computes an integer power of an R16.
!
!  Discussion:
!
!    Obviously, R**P can be computed using P-1 multiplications.
!
!    However, R**P can also be computed using at most 2*LOG2(P) multiplications.
!    To do the calculation this way, let N = LOG2(P).
!    Compute A, A**2, A**4, ..., A**N by N-1 successive squarings.
!    Start the value of R**P at A, and each time that there is a 1 in
!    the binary expansion of P, multiply by the current result of the squarings.
!
!    This algorithm is not optimal.  For small exponents, and for special
!    cases, the result can be computed even more quickly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) R, the base.
!
!    Input, integer ( kind = 4 ) P, the power, which may be negative.
!
!    Output, real ( kind = 16 ) RP, the value of R**P.
!
!    Output, integer ( kind = 4 ) MULTS, the number of multiplications
!    and divisions.
!
  implicit none

  integer ( kind = 4 ) mults
  integer ( kind = 4 ) p
  integer ( kind = 4 ) p_mag
  integer ( kind = 4 ) p_sign
  real ( kind = 16 ) r
  real ( kind = 16 ) r2
  real ( kind = 16 ) rp

  mults = 0
!
!  Special bases.
!
  if ( r == 1.0_16 ) then
    rp = 1.0_16
    return
  end if

  if ( r == -1.0_16 ) then

    if ( mod ( p, 2 ) == 1 ) then
      rp = -1.0_16
    else
      rp = 1.0_16
    end if

    return

  end if

  if ( r == 0.0_16 ) then

    if ( p <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R16_POWER_FAST - Fatal error!'
      write ( *, '(a)' ) '  Base R is zero, and exponent is negative.'
      write ( *, '(a,i8)' ) '  Exponent P = ', p
      stop
    end if

    rp = 0.0_16
    return

  end if
!
!  Special powers.
!
  if ( p == -1 ) then
    rp = 1.0_16 / r
    mults = mults + 1
    return
  else if ( p == 0 ) then
    rp = 1.0_16
    return
  else if ( p == 1 ) then
    rp = r
    return
  end if
!
!  Some work to do.
!
  p_mag = abs ( p )
  p_sign = sign ( 1, p )

  rp = 1.0_16
  r2 = r

  do while ( 0 < p_mag )

    if ( mod ( p_mag, 2 ) == 1 ) then
      rp = rp * r2
      mults = mults + 1
    end if

    p_mag = p_mag / 2
    r2 = r2 * r2
    mults = mults + 1

  end do

  if ( p_sign == -1 ) then
    rp = 1.0_16 / rp
    mults = mults + 1
  end if

  return
end
function r16_pythag ( a, b )

!*****************************************************************************80
!
!! R16_PYTHAG computes sqrt ( A * A + B * B ), avoiding overflow and underflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A, B, the values for which sqrt ( A * A + B * B )
!    is desired.
!
!    Output, real ( kind = 16 ) R16_PYTHAG, the value of sqrt ( A * A + B * B ).
!
  implicit none

  real ( kind = 16 ) a
  real ( kind = 16 ) a_abs
  real ( kind = 16 ) b
  real ( kind = 16 ) b_abs
  real ( kind = 16 ) r16_pythag

  a_abs = abs ( a )
  b_abs = abs ( b )

  if ( b_abs < a_abs ) then
    r16_pythag = a_abs * sqrt ( 1.0_16 + ( b_abs / a_abs ) * ( b_abs / a_abs ) )
  else if ( b_abs == 0.0_16 ) then
    r16_pythag = 0.0_16
  else if ( a_abs <= b_abs ) then
    r16_pythag = b_abs * sqrt ( 1.0_16 + ( a_abs / b_abs ) * ( a_abs / b_abs ) )
  end if

  return
end
subroutine r16_round2 ( nplace, x, xround )

!*****************************************************************************80
!
!! R16_ROUND2 rounds an R16 to a specified number of binary digits.
!
!  Discussion:
!
!    Assume that the input quantity X has the form
!
!      X = S * J * 2^L
!
!    where S is plus or minus 1, L is an integer, and J is a binary
!    mantissa which is either exactly zero, or greater than or equal
!    to 0.5 and strictly less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * 2^L
!
!    where S and L are unchanged, and K is a binary mantissa which
!    agrees with J in the first NPLACE binary digits and is zero
!    thereafter.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0 or 0.5.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0, 0.25, 0.50,
!    or 0.75.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPLACE, the number of binary digits to
!    preserve.  NPLACE should be 0 or positive.
!
!    Input, real ( kind = 16 ) X, the number to be decomposed.
!
!    Output, real ( kind = 16 ) XROUND, the rounded value of X.
!
  implicit none

  integer ( kind = 4 ) iplace
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nplace
  integer ( kind = 4 ) s
  real ( kind = 16 ) x
  real ( kind = 16 ) xmant
  real ( kind = 16 ) xround
  real ( kind = 16 ) xtemp

  xround = 0.0_16
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0_16 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign S.
!
  if ( 0.0_16 < x ) then
    s = 1
    xtemp = x
  else
    s = -1
    xtemp = -x
  end if
!
!  3: Force XTEMP to lie between 1 and 2, and compute the
!  logarithm L.
!
  l = 0

  do while ( 2.0_16 <= xtemp )
    xtemp = xtemp / 2.0_16
    l = l + 1
  end do

  do while ( xtemp < 1.0_16 )
    xtemp = xtemp * 2.0_16
    l = l - 1
  end do
!
!  4: Strip out the digits of the mantissa as XMANT, and decrease L.
!
  xmant = 0.0_16
  iplace = 0

  do

    xmant = 2.0_16 * xmant

    if ( 1.0_16 <= xtemp ) then
      xmant = xmant + 1.0_16
      xtemp = xtemp - 1.0_16
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0_16 .or. nplace <= iplace ) then
      xround = s * xmant * 2.0_16**l
      exit
    end if

    l = l - 1
    xtemp = xtemp * 2.0_16

  end do

  return
end
subroutine r16_roundb ( base, nplace, x, xround )

!*****************************************************************************80
!
!! R16_ROUNDB rounds an R16 to a given number of digits in a given base.
!
!  Discussion:
!
!    The code does not seem to do a good job of rounding when
!    the base is negative!
!
!    Assume that the input quantity X has the form
!
!      X = S * J * BASE^L
!
!    where S is plus or minus 1, L is an integer, and J is a
!    mantissa base BASE which is either exactly zero, or greater
!    than or equal to (1/BASE) and less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * BASE^L
!
!    where S and L are unchanged, and K is a mantissa base BASE
!    which agrees with J in the first NPLACE digits and is zero
!    thereafter.
!
!    Note that because of rounding, for most bases, most numbers
!    with a fractional quantities cannot be stored exactly in the
!    computer, and hence will have trailing "bogus" digits.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0,
!    1/BASE, 2/BASE, ..., (BASE-1)/BASE.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0,
!    BASE/BASE^2, (BASE+1)/BASE^2, ...,
!    BASE^2-2/BASE^2, BASE^2-1/BASE^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BASE, the base of the arithmetic.
!    BASE must not be zero.  Theoretically, BASE may be negative.
!
!    Input, integer ( kind = 4 ) NPLACE, the number of digits base BASE to
!    preserve.  NPLACE should be 0 or positive.
!
!    Input, real ( kind = 16 ) X, the number to be decomposed.
!
!    Output, real ( kind = 16 ) XROUND, the rounded value of X.
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) iplace
  integer ( kind = 4 ) is
  integer ( kind = 4 ) js
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nplace
  real ( kind = 16 ) x
  real ( kind = 16 ) xmant
  real ( kind = 16 ) xround
  real ( kind = 16 ) xtemp

  xround = 0.0_16
!
!  0: Error checks.
!
  if ( base == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R16_ROUNDB - Fatal error!'
    write ( *, '(a)' ) '  The base BASE cannot be zero.'
    stop
  end if
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0_16 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign IS.
!
  if ( 0.0_16 < x ) then
    is = 1
    xtemp = x
  else
    is = -1
    xtemp = -x
  end if
!
!  3: Force XTEMP to lie between 1 and ABS(BASE), and compute the
!  logarithm L.
!
  l = 0

  do while ( abs ( base ) <= abs ( xtemp ) )

    xtemp = xtemp / real ( base, kind = 16 )

    if ( xtemp < 0.0_16 ) then
      is = -is
      xtemp = -xtemp
    end if

    l = l + 1

  end do

  do while ( abs ( xtemp ) < 1.0_16 )

    xtemp = xtemp * base

    if ( xtemp < 0.0_16 ) then
      is = -is
      xtemp = -xtemp
    end if

    l = l - 1

  end do
!
!  4: Now strip out the digits of the mantissa as XMANT, and
!  decrease L.
!
  xmant = 0.0_16
  iplace = 0
  js = is

  do

    xmant = base * xmant

    if ( xmant < 0.0_16 ) then
      js = - js
      xmant = - xmant
    end if

    if ( 1.0_16 <= xtemp ) then
      xmant = xmant + int ( xtemp )
      xtemp = xtemp - int ( xtemp )
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0_16 .or. nplace <= iplace ) then
      xround = js * xmant * real ( base, kind = 16 )**l
      exit
    end if

    l = l - 1
    xtemp = xtemp * base

    if ( xtemp < 0.0_16 ) then
      is = - is
      xtemp = - xtemp
    end if

  end do

  return
end
subroutine r16_roundx ( nplace, x, xround )

!*****************************************************************************80
!
!! R16_ROUNDX rounds an R16.
!
!  Discussion:
!
!    Assume that the input quantity X has the form
!
!      X = S * J * 10^L
!
!    where S is plus or minus 1, L is an integer, and J is a decimal
!    mantissa which is either exactly zero, or greater than or equal
!    to 0.1 and less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * 10^L
!
!    where S and L are unchanged, and K is a decimal mantissa which
!    agrees with J in the first NPLACE decimal digits and is zero
!    thereafter.
!
!    Note that because of rounding, most decimal fraction quantities
!    cannot be stored exactly in the computer, and hence will have
!    trailing "bogus" digits.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0, 0.1,
!    0.2, ..., or 0.9.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0, 0.01, 0.02,
!    0.03, ..., 0.98, 0.99.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPLACE, the number of decimal digits to
!    preserve.  NPLACE should be 0 or positive.
!
!    Input, real ( kind = 16 ) X, the number to be decomposed.
!
!    Output, real ( kind = 16 ) XROUND, the rounded value of X.
!
  implicit none

  integer ( kind = 4 ) iplace
  integer ( kind = 4 ) is
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nplace
  real ( kind = 16 ) x
  real ( kind = 16 ) xmant
  real ( kind = 16 ) xround
  real ( kind = 16 ) xtemp

  xround = 0.0_16
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0_16 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign IS.
!
  if ( 0.0_16 < x ) then
    is = 1
    xtemp = x
  else
    is = -1
    xtemp = - x
  end if
!
!  3: Force XTEMP to lie between 1 and 10, and compute the
!  logarithm L.
!
  l = 0

  do while ( 10.0_16 <= x )
    xtemp = xtemp / 10.0_16
    l = l + 1
  end do

  do while ( xtemp < 1.0_16 )
    xtemp = xtemp * 10.0_16
    l = l - 1
  end do
!
!  4: Now strip out the digits of the mantissa as XMANT, and
!  decrease L.
!
  xmant = 0.0_16
  iplace = 0

  do

    xmant = 10.0_16 * xmant

    if ( 1.0_16 <= xtemp ) then
      xmant = xmant + int ( xtemp )
      xtemp = xtemp - int ( xtemp )
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0_16 .or. nplace <= iplace ) then
      xround = is * xmant * ( 10.0_16**l )
      exit
    end if

    l = l - 1
    xtemp = xtemp * 10.0_16

  end do

  return
end
function r16_sign ( x )

!*****************************************************************************80
!
!! R16_SIGN returns the sign of an R16.
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
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number whose sign is desired.
!
!    Output, real ( kind = 16 ) R16_SIGN, the sign of X:
!
  implicit none

  real ( kind = 16 ) r16_sign
  real ( kind = 16 ) x

  if ( x < 0.0_16 ) then
    r16_sign = -1.0_16
  else
    r16_sign = +1.0_16
  end if

  return
end
function r16_sign_opposite ( r1, r2 )

!*****************************************************************************80
!
!! R16_SIGN_OPPOSITE is TRUE if two R16's are not of the same sign.
!
!  Discussion:
!
!    This test could be coded numerically as
!
!      if ( r1 * r2 <= 0.0 ) then ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) R1, R2, the values to check.
!
!    Output, logical R16_SIGN_OPPOSITE, is TRUE if ( R1 <= 0 and 0 <= R2 )
!    or ( R2 <= 0 and 0 <= R1 ).
!
  implicit none

  real ( kind = 16 ) r1
  real ( kind = 16 ) r2
  logical r16_sign_opposite

  r16_sign_opposite = ( r1 <= 0.0_16 .and. 0.0_16 <= r2 ) .or. &
                      ( r2 <= 0.0_16 .and. 0.0_16 <= r1 )

  return
end
function r16_sign_opposite_strict ( r1, r2 )

!*****************************************************************************80
!
!! R16_SIGN_OPPOSITE_STRICT is TRUE if two R16's are strictly of opposite sign.
!
!  Discussion:
!
!    This test could be coded numerically as
!
!      if ( r1 * r2 < 0.0 ) then ...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) R1, R2, the values to check.
!
!    Output, logical R16_SIGN_OPPOSITE_STRICT, is TRUE if ( R1 < 0 and 0 < R2 )
!    or ( R2 < 0 and 0 < R1 ).
!
  implicit none

  real ( kind = 16 ) r1
  real ( kind = 16 ) r2
  logical r16_sign_opposite_strict

  r16_sign_opposite_strict = ( r1 < 0.0_16 .and. 0.0_16 < r2 ) .or. &
                             ( r2 < 0.0_16 .and. 0.0_16 < r1 )

  return
end
function r16_sub ( x, y )

!*****************************************************************************80
!
!! R16_SUB returns the difference of two R16's.
!
!  Discussion:
!
!    An R16 is a real ( kind = 16 ) value.
!
!    FORTRAN90 supplies the - operator, which should generally be used instead
!    of this function!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, Y, the numbers to be subtracted.
!
!    Output, real ( kind = 16 ) R16_SUB, the difference.
!
  implicit none

  real ( kind = 16 ) r16_sub
  real ( kind = 16 ) x
  real ( kind = 16 ) y

  r16_sub = x - y

  return
end
subroutine r16_swap ( x, y )

!*****************************************************************************80
!
!! R16_SWAP swaps two R16's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 16 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 16 ) x
  real ( kind = 16 ) y
  real ( kind = 16 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r16_swap3 ( x, y, z )

!*****************************************************************************80
!
!! R16_SWAP3 swaps three R16's.
!
!  Example:
!
!    Input:
!
!      X = 1, Y = 2, Z = 3
!
!    Output:
!
!      X = 2, Y = 3, Z = 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 16 ) X, Y, Z, three values to be swapped.
!
  implicit none

  real ( kind = 16 ) w
  real ( kind = 16 ) x
  real ( kind = 16 ) y
  real ( kind = 16 ) z

  w = x
  x = y
  y = z
  z = w

  return
end
function r16_tiny ( )

!*****************************************************************************80
!
!! R16_TINY returns a very small but positive R16.
!
!  Discussion:
!
!    FORTRAN90 provides a built-in routine TINY ( X ) that
!    is more suitable for this purpose, returning the smallest positive
!    but normalized real number.
!
!    This routine does NOT try to provide an accurate value for TINY.
!    Instead, it simply returns a "reasonable" value, that is, a rather
!    small, but representable, real number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 16 ) R16_TINY, a "tiny" value.
!
  implicit none

  real ( kind = 16 ) r16_tiny

  r16_tiny = 1.0E-30_16

  return
end
subroutine r16_to_r16_discrete ( r, rmin, rmax, nr, rd )

!*****************************************************************************80
!
!! R16_TO_R16_DISCRETE maps R to RD in [RMIN, RMAX] with NR possible values.
!
!  Formula:
!
!    if ( R < RMIN ) then
!      RD = RMIN
!    else if ( RMAX < R ) then
!      RD = RMAX
!    else
!      T = nint ( ( NR - 1 ) * ( R - RMIN ) / ( RMAX - RMIN ) )
!      RD = RMIN + T * ( RMAX - RMIN ) / real ( NR - 1 )
!
!    In the special case where NR = 1, when
!
!      XD = 0.5 * ( RMAX + RMIN )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) R, the number to be converted.
!
!    Input, real ( kind = 16 ) RMAX, RMIN, the maximum and minimum
!    values for RD.
!
!    Input, integer ( kind = 4 ) NR, the number of allowed values for XD.
!    NR should be at least 1.
!
!    Output, real ( kind = 16 ) RD, the corresponding discrete value.
!
  implicit none

  integer ( kind = 4 ) f
  integer ( kind = 4 ) nr
  real ( kind = 16 ) r
  real ( kind = 16 ) rd
  real ( kind = 16 ) rmax
  real ( kind = 16 ) rmin
!
!  Check for errors.
!
  if ( nr < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R16_TO_R16_DISCRETE - Fatal error!'
    write ( *, '(a,i8)' ) '  NR = ', nr
    write ( *, '(a)' ) '  but NR must be at least 1.'
    stop
  end if

  if ( nr == 1 ) then
    rd = 0.5_16 * ( rmin + rmax )
    return
  end if

  if ( rmax == rmin ) then
    rd = rmax
    return
  end if

  f = nint ( real ( nr, kind = 16 ) * ( rmax - r ) / ( rmax - rmin ) )
  f = max ( f, 0 )
  f = min ( f, nr )

  rd = ( real (      f, kind = 16 ) * rmin   &
       + real ( nr - f, kind = 16 ) * rmax ) &
       / real ( nr,     kind = 16 )

  return
end
subroutine r16_to_dhms ( r, d, h, m, s )

!*****************************************************************************80
!
!! R16_TO_DHMS converts decimal days into days, hours, minutes, seconds.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) R, a decimal number representing a time
!    period measured in days.
!
!    Output, integer ( kind = 4 ) D, H, M, S, the equivalent number of days,
!    hours, minutes and seconds.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  real ( kind = 16 ) r
  real ( kind = 16 ) r_copy
  integer ( kind = 4 ) s

  r_copy = abs ( r )

  d = int ( r_copy )

  r_copy = r_copy - d
  r_copy = 24.0_16 * r_copy
  h = int ( r_copy )

  r_copy = r_copy - h
  r_copy = 60.0_16 * r_copy
  m = int ( r_copy )

  r_copy = r_copy - m
  r_copy = 60.0_16 * r_copy
  s = int ( r_copy )

  if ( r < 0.0_16 ) then
    d = -d
    h = -h
    m = -m
    s = -s
  end if

  return
end
subroutine r16_to_i4 ( x, xmin, xmax, ixmin, ixmax, ix )

!*****************************************************************************80
!
!! R16_TO_I4 maps X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
!
!  Formula:
!
!    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
!    IX := min ( IX, max ( IXMIN, IXMAX ) )
!    IX := max ( IX, min ( IXMIN, IXMAX ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the number to be converted.
!
!    Input, real ( kind = 16 ) XMIN, XMAX, the range.  XMAX and
!    XMIN must not be equal.  It is not necessary that XMIN be less than XMAX.
!
!    Input, integer ( kind = 4 ) IXMIN, IXMAX, the allowed range of the output
!    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
!    It is not necessary that IXMIN be less than IXMAX.
!
!    Output, integer ( kind = 4 ) IX, the value in the range [IXMIN,IXMAX] that
!    corresponds to X.
!
  implicit none

  integer ( kind = 4 ) ix
  integer ( kind = 4 ) ixmax
  integer ( kind = 4 ) ixmin
  real ( kind = 16 ) temp
  real ( kind = 16 ) x
  real ( kind = 16 ) xmax
  real ( kind = 16 ) xmin

  if ( xmax == xmin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R16_TO_I4 - Fatal error!'
    write ( *, '(a)' ) '  XMAX = XMIN, making a zero divisor.'
    write ( *, '(a,g14.6)' ) '  XMAX = ', xmax
    write ( *, '(a,g14.6)' ) '  XMIN = ', xmin
    stop
  end if

  temp = &
      ( ( xmax - x        ) * real ( ixmin, kind = 16 )  &
      + (        x - xmin ) * real ( ixmax, kind = 16 ) ) &
      / ( xmax     - xmin )

  if ( 0.0_16 <= temp ) then
    temp = temp + 0.5_16
  else
    temp = temp - 0.5_16
  end if

  ix = int ( temp )

  return
end
function r16_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R16_UNIFORM returns a scaled pseudorandom R16.
!
!  Discussion:
!
!    An R16 is a real ( kind = 16 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 16 ) R16_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 16 ) a
  real ( kind = 16 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 16 ) r16_uniform
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R16_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r16_uniform = a + ( b - a ) * real ( seed, kind = 16 ) * 4.656612875E-10_16

  return
end
function r16_uniform_01 ( seed )

!*****************************************************************************80
!
!! R16_UNIFORM_01 returns a unit pseudorandom R16.
!
!  Discussion:
!
!    An R16 is a real ( kind = 16 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r16_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R16_UNIFORM_01
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
!    23 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 16 ) R16_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 16 ) r16_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R16_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r16_uniform_01 = real ( seed, kind = 16 ) * 4.656612875E-10_16

  return
end
subroutine r16_unswap3 ( x, y, z )

!*****************************************************************************80
!
!! R16_UNSWAP3 unswaps three R16's.
!
!  Example:
!
!    Input:
!
!      X = 2, Y = 3, Z = 1
!
!    Output:
!
!      X = 1, Y = 2, Z = 3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 16 ) X, Y, Z, three values to be swapped.
!
  implicit none

  real ( kind = 16 ) w
  real ( kind = 16 ) x
  real ( kind = 16 ) y
  real ( kind = 16 ) z

  w = z
  z = y
  y = x
  x = w

  return
end
function r16_walsh_1d ( x, digit )

!*****************************************************************************80
!
!! R16_WALSH_1D evaluates the Walsh function.
!
!  Discussion:
!
!    Consider the binary representation of X, and number the digits
!    in descending order, from leading to lowest, with the units digit
!    being numbered 0.
!
!    The Walsh function W(J)(X) is equal to the J-th binary digit of X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) X, the argument of the Walsh function.
!
!    Input, integer ( kind = 4 ) DIGIT, the index of the Walsh function.
!
!    Output, real ( kind = 16 ) R16_WALSH_1D, the value of the Walsh function.
!
  implicit none

  integer ( kind = 4 ) digit
  integer ( kind = 4 ) n
  real ( kind = 16 ) r16_walsh_1d
  real ( kind = 16 ) x
  real ( kind = 16 ) x_copy
!
!  Hide the effect of the sign of X.
!
  x_copy = abs ( x )
!
!  If DIGIT is positive, divide by 2 DIGIT times.
!  If DIGIT is negative, multiply by 2 (-DIGIT) times.
!
  x_copy = x_copy / 2.0_16 ** digit
!
!  Make it an integer.
!  Because it's positive, and we're using INT, we don't change the
!  units digit.
!
  n = int ( x_copy )
!
!  Is the units digit odd or even?
!
  if ( mod ( n, 2 ) == 0 ) then
    r16_walsh_1d = 0.0_16
  else
    r16_walsh_1d = 1.0_16
  end if

  return
end
subroutine r16mat_border_add ( m, n, table, table2 )

!*****************************************************************************80
!
!! R16MAT_BORDER_ADD adds a "border" to an R16MAT.
!
!  Discussion:
!
!    We suppose the input data gives values of a quantity on nodes
!    in the interior of a 2D grid, and we wish to create a new table
!    with additional positions for the nodes that would be on the
!    border of the 2D grid.
!
!                  0 0 0 0 0 0
!      * * * *     0 * * * * 0
!      * * * * --> 0 * * * * 0
!      * * * *     0 * * * * 0
!                  0 0 0 0 0 0
!
!    The illustration suggests the situation in which a 3 by 4 array
!    is input, and a 5 by 6 array is to be output.
!
!    The old data is shifted to its correct positions in the new array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 16 ) TABLE(M,N), the table data.
!
!    Output, real ( kind = 16 ) TABLE2(M+2,N+2), the augmented table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) table(m,n)
  real ( kind = 16 ) table2(m+2,n+2)

  table2(1,1:n+2) = 0.0_16
  table2(m+2,1:n+2) = 0.0_16
  table2(2:m+1,1) = 0.0_16
  table2(2:m+1,n+2) = 0.0_16

  table2(2:m+1,2:n+1) = table(1:m,1:n)

  return
end
subroutine r16mat_border_cut ( m, n, table, table2 )

!*****************************************************************************80
!
!! R16MAT_BORDER_CUT cuts the "border" of an R16MAT.
!
!  Discussion:
!
!    We suppose the input data gives values of a quantity on nodes
!    on a 2D grid, and we wish to create a new table corresponding only
!    to those nodes in the interior of the 2D grid.
!
!      0 0 0 0 0 0
!      0 * * * * 0    * * * *
!      0 * * * * 0 -> * * * *
!      0 * * * * 0    * * * *
!      0 0 0 0 0 0
!
!    The illustration suggests the situation in which a 5 by 6 array
!    is input, and a 3 by 4 array is to be output.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 16 ) TABLE(M,N), the table data.
!
!    Output, real ( kind = 16 ) TABLE2(M-2,N-2), the new table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) table(m,n)
  real ( kind = 16 ) table2(m-2,n-2)

  if ( m <= 2 .or. n <= 2 ) then
    return
  end if

  table2(1:m-2,1:n-2) = table(2:m-1,2:n-1)

  return
end
subroutine r16mat_cholesky_factor ( n, a, c )

!*****************************************************************************80
!
!! R16MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric matrix.
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is a lower triangular matrix L such that:
!
!      A = L * L'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 16 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 16 ) C(N,N), the N by N lower triangular
!    Cholesky factor.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) c(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 16 ) sum2

  c(1:n,1:n) = a(1:n,1:n)

  do j = 1, n

    c(1:j-1,j) = 0.0_16

    do i = j, n

      sum2 = c(j,i) - dot_product ( c(j,1:j-1), c(i,1:j-1) )

      if ( i == j ) then
        if ( sum2 <= 0.0_16 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R16MAT_CHOLESKY_FACTOR - Fatal error!'
          write ( *, '(a)' ) '  Matrix is not positive definite.'
          stop
        else
          c(i,j) = sqrt ( sum2 )
        end if
      else
        if ( c(j,j) /= 0.0_16 ) then
          c(i,j) = sum2 / c(j,j)
        else
          c(i,j) = 0.0_16
        end if
      end if

    end do

  end do

  return
end
subroutine r16mat_cholesky_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R16MAT_CHOLESKY_SOLVE solves a Cholesky factored linear system A * x = b.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 16 ) A(N,N), the N by N Cholesky factor of the
!    system matrix.
!
!    Input, real ( kind = 16 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 16 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) b(n)
  real ( kind = 16 ) x(n)
!
!  Solve L * y = b.
!
  call r16mat_l_solve ( n, a, b, x )
!
!  Solve L' * x = y.
!
  call r16mat_lt_solve ( n, a, x, x )

  return
end
subroutine r16mat_choresky_factor ( n, a, c )

!*****************************************************************************80
!
!! R16MAT_CHORESKY_FACTOR computes the "Choresky" factor of a symmetric matrix.
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is an upper triangular matrix R such that:
!
!      A = R * R'
!
!    Note that the usual Cholesky factor is a LOWER triangular matrix L
!    such that
!
!      A = L * L'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 16 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 16 ) C(N,N), the N by N upper triangular
!    "Choresky" factor.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) c(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 16 ) sum2

  c(n:1:-1,n:1:-1) = a(1:n,1:n)

  do j = 1, n

    c(1:j-1,j) = 0.0_16

    do i = j, n

      sum2 = c(j,i) - dot_product ( c(j,1:j-1), c(i,1:j-1) )

      if ( i == j ) then
        if ( sum2 <= 0.0_16 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R16MAT_CHORESKY_FACTOR - Fatal error!'
          write ( *, '(a)' ) '  Matrix is not positive definite.'
          stop
        else
          c(i,j) = sqrt ( sum2 )
        end if
      else
        if ( c(j,j) /= 0.0_16 ) then
          c(i,j) = sum2 / c(j,j)
        else
          c(i,j) = 0.0_16
        end if
      end if

    end do

  end do

  c(n:1:-1,n:1:-1) = c(1:n,1:n)

  return
end
subroutine r16mat_copy ( m, n, a, b )

!*****************************************************************************80
!
!! R16MAT_COPY copies an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix to be copied.
!
!    Output, real ( kind = 16 ) B(M,N), a copy of the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  real ( kind = 16 ) b(m,n)

  b(1:m,1:n) = a(1:m,1:n)

  return
end
subroutine r16mat_det ( n, a, det )

!*****************************************************************************80
!
!! R16MAT_DET computes the determinant of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
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
!    Input, real ( kind = 16 ) A(N,N), the matrix whose determinant is desired.
!
!    Output, real ( kind = 16 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) b(n,n)
  real ( kind = 16 ) det
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) piv(1)
  real ( kind = 16 ) t

  b(1:n,1:n) = a(1:n,1:n)

  det = 1.0_16

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

    if ( b(k,k) /= 0.0_16 ) then

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
function r16mat_det_2d ( a )

!*****************************************************************************80
!
!! R16MAT_DET_2D computes the determinant of a 2 by 2 R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The formula for the determinant of a 2 by 2 matrix is
!
!      a11 * a22 - a12 * a21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A(2,2), the matrix whose determinant is desired.
!
!    Output, real ( kind = 16 ) R16MAT_DET_2D, the determinant of the matrix.
!
  implicit none

  real ( kind = 16 ) a(2,2)
  real ( kind = 16 ) r16mat_det_2d

  r16mat_det_2d = a(1,1) * a(2,2) - a(1,2) * a(2,1)

  return
end
function r16mat_det_3d ( a )

!*****************************************************************************80
!
!! R16MAT_DET_3D computes the determinant of a 3 by 3 R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The formula for the determinant of a 3 by 3 matrix is
!
!        a11 * a22 * a33 - a11 * a23 * a32
!      + a12 * a23 * a31 - a12 * a21 * a33
!      + a13 * a21 * a32 - a13 * a22 * a31
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A(3,3), the matrix whose determinant is desired.
!
!    Output, real ( kind = 16 ) R16MAT_DET_3D, the determinant of the matrix.
!
  implicit none

  real ( kind = 16 ) a(3,3)
  real ( kind = 16 ) r16mat_det_3d

  r16mat_det_3d = &
         a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
       + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
       + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

  return
end
function r16mat_det_4d ( a )

!*****************************************************************************80
!
!! R16MAT_DET_4D computes the determinant of a 4 by 4 R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A(4,4), the matrix whose determinant is desired.
!
!    Output, real ( kind = 16 ) R16MAT_DET_4D, the determinant of the matrix.
!
  implicit none

  real ( kind = 16 ) a(4,4)
  real ( kind = 16 ) r16mat_det_4d

  r16mat_det_4d = &
         a(1,1) * ( &
             a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
       - a(1,2) * ( &
             a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
       + a(1,3) * ( &
             a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
       - a(1,4) * ( &
             a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
           + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

  return
end
function r16mat_det_5d ( a )

!*****************************************************************************80
!
!! R16MAT_DET_5D computes the determinant of a 5 by 5 R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A(5,5), the matrix whose determinant is desired.
!
!    Output, real ( kind = 16 ) R16MAT_DET_5D, the determinant of the matrix.
!
  implicit none

  real ( kind = 16 ) a(5,5)
  real ( kind = 16 ) b(4,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 16 ) r16mat_det_4d
  real ( kind = 16 ) r16mat_det_5d
!
!  Expand the determinant into the sum of the determinants of the
!  five 4 by 4 matrices created by dropping row 1, and column k.
!
  r16mat_det_5d = 0.0_16

  do k = 1, 5

    do i = 1, 4
      do j = 1, 4

        if ( j < k ) then
          inc = 0
        else
          inc = 1
        end if

        b(i,j) = a(i+1,j+inc)

      end do
    end do

    r16mat_det_5d = r16mat_det_5d &
      + (-1)**( k + 1 ) * a(1,k) * r16mat_det_4d ( b )

  end do

  return
end
subroutine r16mat_diag_add_scalar ( n, a, s )

!*****************************************************************************80
!
!! R16MAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns.
!
!    Input/output, real ( kind = 16 ) A(N,N), the N by N matrix to be modified.
!
!    Input, real ( kind = 16 ) S, the value to be added to the diagonal
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 16 ) s

  do i = 1, n
    a(i,i) = a(i,i) + s
  end do

  return
end
subroutine r16mat_diag_add_vector ( n, a, v )

!*****************************************************************************80
!
!! R16MAT_DIAG_ADD_VECTOR adds a vector to the diagonal of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix.
!
!    Input/output, real ( kind = 16 ) A(N,N), the N by N matrix.
!
!    Input, real ( kind = 16 ) V(N), the vector to be added to the diagonal
!    of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 16 ) v(n)

  do i = 1, n
    a(i,i) = a(i,i) + v(i)
  end do

  return
end
subroutine r16mat_diag_get_vector ( n, a, v )

!*****************************************************************************80
!
!! R16MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix.
!
!    Input, real ( kind = 16 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 16 ) V(N), the diagonal entries
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 16 ) v(n)

  do i = 1, n
    v(i) = a(i,i)
  end do

  return
end
subroutine r16mat_diag_set_scalar ( n, a, s )

!*****************************************************************************80
!
!! R16MAT_DIAG_SET_SCALAR sets the diagonal of an R16MAT to a scalar value.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns.
!
!    Input/output, real ( kind = 16 ) A(N,N), the N by N matrix to be modified.
!
!    Input, real ( kind = 16 ) S, the value to be assigned to the diagonal
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 16 ) s

  do i = 1, n
    a(i,i) = s
  end do

  return
end
subroutine r16mat_diag_set_vector ( n, a, v )

!*****************************************************************************80
!
!! R16MAT_DIAG_SET_VECTOR sets the diagonal of an R16MAT to a vector.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns.
!
!    Input/output, real ( kind = 16 ) A(N,N), the N by N matrix.
!
!    Input, real ( kind = 16 ) V(N), the vector to be assigned to the
!    diagonal of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 16 ) v(n)

  do i = 1, n
    a(i,i) = v(i)
  end do

  return
end
subroutine r16mat_expand_linear ( m, n, x, mfat, nfat, xfat )

!*****************************************************************************80
!
!! R16MAT_EXPAND_LINEAR linearly interpolates new data into an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    In this routine, the expansion is specified by giving the number
!    of intermediate values to generate between each pair of original
!    data rows and columns.
!
!    The interpolation is not actually linear.  It uses the functions
!
!      1, x, y, and xy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    input data.
!
!    Input, real ( kind = 16 ) X(M,N), the original data.
!
!    Input, integer ( kind = 4 ) MFAT, NFAT, the number of data values
!    to interpolate between each row, and each column, of original data values.
!
!    Output, real ( kind = 16 ) XFAT(M2,N2), the fattened data, where
!    M2 = (M-1)*(MFAT+1)+1,
!    N2 = (N-1)*(NFAT+1)+1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mfat
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nfat

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jjj
  integer ( kind = 4 ) jp1
  real ( kind = 16 ) s
  real ( kind = 16 ) t
  real ( kind = 16 ) x(m,n)
  real ( kind = 16 ) x00
  real ( kind = 16 ) x01
  real ( kind = 16 ) x10
  real ( kind = 16 ) x11
  real ( kind = 16 ) xfat((m-1)*(mfat+1)+1,(n-1)*(nfat+1)+1)

  do i = 1, m

    if ( i < m ) then
      ihi = mfat
    else
      ihi = 0
    end if

    do j = 1, n

      if ( j < n ) then
        jhi = nfat
      else
        jhi = 0
      end if

      if ( i < m ) then
        ip1 = i + 1
      else
        ip1 = i
      end if

      if ( j < n ) then
        jp1 = j + 1
      else
        jp1 = j
      end if

      x00 = x(i,j)
      x10 = x(ip1,j)
      x01 = x(i,jp1)
      x11 = x(ip1,jp1)

      do ii = 0, ihi

        s = real ( ii, kind = 16 ) &
          / real ( ihi + 1, kind = 16 )

        do jj = 0, jhi

          t = real ( jj, kind = 16 ) &
            / real ( jhi + 1, kind = 16 )

          iii = 1 + ( i - 1 ) * ( mfat + 1 ) + ii
          jjj = 1 + ( j - 1 ) * ( nfat + 1 ) + jj

          xfat(iii,jjj) = &
                                            x00   &
              + s     * (       x10       - x00 ) &
              + t     * (             x01 - x00 ) &
              + s * t * ( x11 - x10 - x01 + x00 )

        end do

      end do

    end do

  end do

  return
end
subroutine r16mat_expand_linear2 ( m, n, a, m2, n2, a2 )

!*****************************************************************************80
!
!! R16MAT_EXPAND_LINEAR2 expands an R16MAT by linear interpolation.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    In this version of the routine, the expansion is indicated
!    by specifying the dimensions of the expanded array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), a "small" M by N array.
!
!    Input, integer ( kind = 4 ) M2, N2, the number of rows and columns in A2.
!
!    Output, real ( kind = 16 ) A2(M2,N2), the expanded array, which
!    contains an interpolated version of the data in A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2

  real ( kind = 16 ) a(m,n)
  real ( kind = 16 ) a2(m2,n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 16 ) r
  real ( kind = 16 ) r1
  real ( kind = 16 ) r2
  real ( kind = 16 ) s
  real ( kind = 16 ) s1
  real ( kind = 16 ) s2

  do i = 1, m2

    if ( m2 == 1 ) then
      r = 0.5_16
    else
      r = real ( i - 1, kind = 16 ) &
        / real ( m2 - 1, kind = 16 )
    end if

    i1 = 1 + int ( r * real ( m - 1, kind = 16 ) )
    i2 = i1 + 1

    if ( m < i2 ) then
      i1 = m - 1
      i2 = m
    end if

    r1 = real ( i1 - 1, kind = 16 ) &
       / real ( m - 1, kind = 16 )

    r2 = real ( i2 - 1, kind = 16 ) &
       / real ( m - 1, kind = 16 )

    do j = 1, n2

      if ( n2 == 1 ) then
        s = 0.5_16
      else
        s = real ( j - 1, kind = 16 ) &
          / real ( n2 - 1, kind = 16 )
      end if

      j1 = 1 + int ( s * real ( n - 1, kind = 16 ) )
      j2 = j1 + 1

      if ( n < j2 ) then
        j1 = n - 1
        j2 = n
      end if

      s1 = real ( j1 - 1, kind = 16 ) &
         / real ( n - 1, kind = 16 )

      s2 = real ( j2 - 1, kind = 16 ) &
         / real ( n - 1, kind = 16 )

      a2(i,j) = &
        ( ( r2 - r ) * ( s2 - s ) * a(i1,j1) &
        + ( r - r1 ) * ( s2 - s ) * a(i2,j1) &
        + ( r2 - r ) * ( s - s1 ) * a(i1,j2) &
        + ( r - r1 ) * ( s - s1 ) * a(i2,j2) ) &
        / ( ( r2 - r1 ) * ( s2 - s1 ) )

    end do

  end do

  return
end
subroutine r16mat_givens_post ( n, a, row, col, g )

!*****************************************************************************80
!
!! R16MAT_GIVENS_POST computes the Givens postmultiplier rotation matrix.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The Givens post-multiplier matrix G(ROW,COL) has the property that
!    the (ROW,COL)-th entry of A*G is zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices A and G.
!
!    Input, real ( kind = 16 ) A(N,N), the matrix to be operated upon.
!
!    Input, integer ( kind = 4 ) ROW, COL, the row and column of the
!    entry of A*G which is to be zeroed out.
!
!    Output, real ( kind = 16 ) G(N,N), the Givens rotation matrix.
!    G is an orthogonal matrix, that is, the inverse of
!    G is the transpose of G.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) col
  real ( kind = 16 ) g(n,n)
  integer ( kind = 4 ) row
  real ( kind = 16 ) theta

  call r16mat_identity ( n, g )

  theta = atan2 ( a(row,col), a(row,row) )

  g(row,row) =  cos ( theta )
  g(row,col) = -sin ( theta )
  g(col,row) =  sin ( theta )
  g(col,col) =  cos ( theta )

  return
end
subroutine r16mat_givens_pre ( n, a, row, col, g )

!*****************************************************************************80
!
!! R16MAT_GIVENS_PRE computes the Givens premultiplier rotation matrix.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The Givens premultiplier rotation matrix G(ROW,COL) has the
!    property that the (ROW,COL)-th entry of G*A is zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices A and G.
!
!    Input, real ( kind = 16 ) A(N,N), the matrix to be operated upon.
!
!    Input, integer ( kind = 4 ) ROW, COL, the row and column of the
!    entry of the G*A which is to be zeroed out.
!
!    Output, real ( kind = 16 ) G(N,N), the Givens rotation matrix.
!    G is an orthogonal matrix, that is, the inverse of
!    G is the transpose of G.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) col
  real ( kind = 16 ) g(n,n)
  integer ( kind = 4 ) row
  real ( kind = 16 ) theta

  call r16mat_identity ( n, g )

  theta = atan2 ( a(row,col), a(col,col) )

  g(row,row) =  cos ( theta )
  g(row,col) = -sin ( theta )
  g(col,row) =  sin ( theta )
  g(col,col) =  cos ( theta )

  return
end
subroutine r16mat_hess ( fx, n, x, h )

!*****************************************************************************80
!
!! R16MAT_HESS approximates a Hessian matrix via finite differences.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    H(I,J) = d2 F / d X(I) d X(J)
!
!    The values returned by this routine will be only approximate.
!    In some cases, they will be so poor that they are useless.
!    However, one of the best applications of this routine is for
!    checking your own Hessian calculations, since as Heraclitus
!    said, you'll never get the same result twice when you differentiate
!    a complicated expression by hand.
!
!    The user function routine, here called "FX", should have the form:
!
!      subroutine fx ( n, x, f )
!      integer n
!      real ( kind = 16 ) f
!      real ( kind = 16 ) x(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FX, the name of the user function routine.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 16 ) X(N), the values of the variables.
!
!    Output, real ( kind = 16 ) H(N,N), the approximated N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) eps
  real ( kind = 16 ) f00
  real ( kind = 16 ) fmm
  real ( kind = 16 ) fmp
  real ( kind = 16 ) fpm
  real ( kind = 16 ) fpp
  external fx
  real ( kind = 16 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 16 ) s(n)
  real ( kind = 16 ) x(n)
  real ( kind = 16 ) xi
  real ( kind = 16 ) xj
!
!  Choose the stepsizes.
!
  eps = ( epsilon ( eps ) )**0.33_16

  do i = 1, n
    s(i) = eps * max ( abs ( x(i) ), 1.0_16 )
  end do
!
!  Calculate the diagonal elements.
!
  do i = 1, n

    xi = x(i)

    call fx ( n, x, f00 )

    x(i) = xi + s(i)
    call fx ( n, x, fpp )

    x(i) = xi - s(i)
    call fx ( n, x, fmm )

    h(i,i) = ( ( fpp - f00 ) + ( fmm - f00 ) ) / s(i)**2

    x(i) = xi

  end do
!
!  Calculate the off diagonal elements.
!
  do i = 1, n

    xi = x(i)

    do j = i + 1, n

      xj = x(j)

      x(i) = xi + s(i)
      x(j) = xj + s(j)
      call fx ( n, x, fpp )

      x(i) = xi + s(i)
      x(j) = xj - s(j)
      call fx ( n, x, fpm )

      x(i) = xi - s(i)
      x(j) = xj + s(j)
      call fx ( n, x, fmp )

      x(i) = xi - s(i)
      x(j) = xj - s(j)
      call fx ( n, x, fmm )

      h(j,i) = ( ( fpp - fpm ) + ( fmm - fmp ) ) / ( 4.0_16 * s(i) * s(j) )

      h(i,j) = h(j,i)

      x(j) = xj

    end do

    x(i) = xi

  end do

  return
end
subroutine r16mat_house_axh ( n, a, v, ah )

!*****************************************************************************80
!
!! R16MAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The Householder matrix H(V) is defined by
!
!      H(V) = I - 2 * v * v' / ( v' * v )
!
!    This routine is not particularly efficient.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 16 ) A(N,N), the matrix to be postmultiplied.
!
!    Input, real ( kind = 16 ) V(N), a vector defining a Householder matrix.
!
!    Output, real ( kind = 16 ) AH(N,N), the product A*H.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) ah(n,n)
  real ( kind = 16 ) ah_temp(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 16 ) v(n)
  real ( kind = 16 ) v_normsq

  v_normsq = sum ( v(1:n)**2 )
!
!  Compute A*H' = A*H
!
  do i = 1, n
    do j = 1, n
      ah_temp(i,j) = a(i,j)
      do k = 1, n
        ah_temp(i,j) = ah_temp(i,j) - 2.0_16 * a(i,k) * v(k) * v(j) / v_normsq
      end do
    end do
  end do
!
!  Copy the temporary result into AH.
!  Doing it this way means the user can identify the input arguments A and AH.
!
  ah(1:n,1:n) = ah_temp(1:n,1:n)

  return
end
subroutine r16mat_house_form ( n, v, h )

!*****************************************************************************80
!
!! R16MAT_HOUSE_FORM constructs a Householder matrix from its compact form.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    H(v) = I - 2 * v * v' / ( v' * v )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 16 ) V(N), the vector defining the Householder matrix.
!
!    Output, real ( kind = 16 ) H(N,N), the Householder matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) beta
  real ( kind = 16 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 16 ) v(n)
!
!  Compute the L2 norm of V.
!
  beta = sum ( v(1:n)**2 )
!
!  Form the matrix H.
!
  call r16mat_identity ( n, h )

  do i = 1, n
    do j = 1, n
      h(i,j) = h(i,j) - 2.0_16 * v(i) * v(j) / beta
    end do
  end do

  return
end
subroutine r16mat_house_hxa ( n, a, v, ha )

!*****************************************************************************80
!
!! R16MAT_HOUSE_HXA computes H*A where H is a compact Householder matrix.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The Householder matrix H(V) is defined by
!
!      H(V) = I - 2 * v * v' / ( v' * v )
!
!    This routine is not particularly efficient.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 16 ) A(N,N), the matrix to be premultiplied.
!
!    Input, real ( kind = 16 ) V(N), a vector defining a Householder matrix.
!
!    Output, real ( kind = 16 ) HA(N,N), the product H*A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) ha(n,n)
  real ( kind = 16 ) ha_temp(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 16 ) v(n)
  real ( kind = 16 ) v_normsq

  v_normsq = sum ( v(1:n)**2 )
!
!  Compute A*H' = A*H
!
  do i = 1, n
    do j = 1, n
      ha_temp(i,j) = a(i,j)
      do k = 1, n
        ha_temp(i,j) = ha_temp(i,j) - 2.0_16 * v(i) * v(k) * a(k,j) / v_normsq
      end do
    end do
  end do
!
!  Copy the temporary result into HA.
!  Doing it this way means the user can identify the input arguments A and HA.
!
  ha(1:n,1:n) = ha_temp(1:n,1:n)

  return
end
subroutine r16mat_house_post ( n, a, row, col, h )

!*****************************************************************************80
!
!! R16MAT_HOUSE_POST computes a Householder post-multiplier matrix.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    H(ROW,COL) has the property that the ROW-th column of
!    A*H(ROW,COL) is zero from entry COL+1 to the end.
!
!    In the most common case, where a QR factorization is being computed,
!    ROW = COL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!
!    Input, real ( kind = 16 ) A(N,N), the matrix whose Householder matrix
!    is to be computed.
!
!    Input, integer ( kind = 4 ) ROW, COL, specify the location of the
!    entry of the matrix A which is to be preserved.  The entries in
!    the same row, but higher column, will be zeroed out if
!    A is postmultiplied by H.
!
!    Output, real ( kind = 16 ) H(N,N), the Householder matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) col
  real ( kind = 16 ) h(n,n)
  integer ( kind = 4 ) row
  real ( kind = 16 ) v(n)
  real ( kind = 16 ) w(n)
!
!  Set up the vector V.
!
  w(1:col-1) = 0.0_16
  w(col:n) = a(row,col:n)

  call r16vec_house_column ( n, w, col, v )
!
!  Form the matrix H(V).
!
  call r16mat_house_form ( n, v, h )

  return
end
subroutine r16mat_house_pre ( n, a, row, col, h )

!*****************************************************************************80
!
!! R16MAT_HOUSE_PRE computes a Householder pre-multiplier matrix.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    H(ROW,COL) has the property that the COL-th column of
!    H(ROW,COL)*A is zero from entry ROW+1 to the end.
!
!    In the most common case, where a QR factorization is being computed,
!    ROW = COL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!
!    Input, real ( kind = 16 ) A(N,N), the matrix whose Householder matrix
!    is to be computed.
!
!    Input, integer ( kind = 4 ) ROW, COL, specify the location of the
!    entry of the matrix A which is to be preserved.  The entries in
!    the same column, but higher rows, will be zeroed out if A is
!    premultiplied by H.
!
!    Output, real ( kind = 16 ) H(N,N), the Householder matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) col
  real ( kind = 16 ) h(n,n)
  integer ( kind = 4 ) row
  real ( kind = 16 ) v(n)
  real ( kind = 16 ) w(n)
!
!  Set up the vector V.
!
  w(1:row-1) = 0.0_16
  w(row:n) = a(row:n,col)

  call r16vec_house_column ( n, w, row, v )
!
!  Form the matrix H(V).
!
  call r16mat_house_form ( n, v, h )

  return
end
subroutine r16mat_identity ( n, a )

!*****************************************************************************80
!
!! R16MAT_IDENTITY stores the identity matrix in an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 16 ) A(N,N), the N by N identity matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i

  a(1:n,1:n) = 0.0_16

  do i = 1, n
    a(i,i) = 1.0_16
  end do

  return
end
function r16mat_in_01 ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_IN_01 is TRUE if the entries of an R16MAT are in the range [0,1].
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix.
!
!    Output, logical R16MAT_IN_01, is TRUE if every entry of A is
!    between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  logical r16mat_in_01

  if ( any ( a(1:m,1:n) < 0.0_16 .or. 1.0_16 < a(1:m,1:n) ) ) then
    r16mat_in_01 = .false.
  else
    r16mat_in_01 = .true.
  end if

  return
end
subroutine r16mat_indicator ( m, n, table )

!*****************************************************************************80
!
!! R16MAT_INDICATOR sets up an "indicator" R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The value of each entry suggests its location, as in:
!
!      11  12  13  14
!      21  22  23  24
!      31  32  33  34
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, real ( kind = 16 ) TABLE(M,N), the table.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  real ( kind = 16 ) table(m,n)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, m
    do j = 1, n
      table(i,j) = real ( fac * i + j, kind = 16 )
    end do
  end do

  return
end
subroutine r16mat_inverse_2d ( a, b, det )

!*****************************************************************************80
!
!! R16MAT_INVERSE_2D inverts a 2 by 2 R16MAT using Cramer's rule.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    If the determinant is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If the determinant is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A(2,2), the matrix to be inverted.
!
!    Output, real ( kind = 16 ) B(2,2), the inverse of the matrix A.
!
!    Output, real ( kind = 16 ) DET, the determinant of the matrix A.
!
  implicit none

  real ( kind = 16 ) a(2,2)
  real ( kind = 16 ) b(2,2)
  real ( kind = 16 ) det
  real ( kind = 16 ) r16mat_det_2d
!
!  Compute the determinant of A.
!
  det = r16mat_det_2d ( a )

  if ( det == 0.0_16 ) then

    b(1:2,1:2) = 0.0_16

  else

    b(1,1) =  a(2,2) / det
    b(1,2) = -a(1,2) / det
    b(2,1) = -a(2,1) / det
    b(2,2) =  a(1,1) / det

  end if

  return
end
subroutine r16mat_inverse_3d ( a, b, det )

!*****************************************************************************80
!
!! R16MAT_INVERSE_3D inverts a 3 by 3 R16MAT using Cramer's rule.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    If the determinant is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If the determinant is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A(3,3), the matrix to be inverted.
!
!    Output, real ( kind = 16 ) B(3,3), the inverse of the matrix A.
!
!    Output, real ( kind = 16 ) DET, the determinant of the matrix A.
!
  implicit none

  real ( kind = 16 ) a(3,3)
  real ( kind = 16 ) b(3,3)
  real ( kind = 16 ) det
  real ( kind = 16 ) r16mat_det_3d
!
!  Compute the determinant of A.
!
  det = r16mat_det_3d ( a )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0_16 ) then
    b(1:3,1:3) = 0.0_16
    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit
!  formula.
!
  b(1,1) =  ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
  b(1,2) = -( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
  b(1,3) =  ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det

  b(2,1) = -( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
  b(2,2) =  ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
  b(2,3) = -( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det

  b(3,1) =  ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
  b(3,2) = -( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
  b(3,3) =  ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det

  return
end
subroutine r16mat_inverse_4d ( a, b, det )

!*****************************************************************************80
!
!! R16MAT_INVERSE_4D inverts a 4 by 4 R16MAT using Cramer's rule.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    If the determinant is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If the determinant is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A(4,4), the matrix to be inverted.
!
!    Output, real ( kind = 16 ) B(4,4), the inverse of the matrix A.
!
!    Output, real ( kind = 16 ) DET, the determinant of the matrix A.
!
  implicit none

  real ( kind = 16 ) a(4,4)
  real ( kind = 16 ) b(4,4)
  real ( kind = 16 ) det
  real ( kind = 16 ) r16mat_det_4d
!
!  Compute the determinant of A.
!
  det = r16mat_det_4d ( a )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0_16 ) then

    b(1:4,1:4) = 0.0_16

    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit formula.
!
  b(1,1) = +( &
        + a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(2,3) * ( a(3,4) * a(4,2) - a(3,2) * a(4,4) ) &
        + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        ) / det

  b(2,1) = -( &
        + a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(2,3) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
        ) / det

  b(3,1) = +( &
        + a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        + a(2,2) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(4,1) = -( &
        + a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        + a(2,2) * ( a(3,3) * a(4,1) - a(3,1) * a(4,3) ) &
        + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(1,2) = -( &
        + a(1,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(1,3) * ( a(3,4) * a(4,2) - a(3,2) * a(4,4) ) &
        + a(1,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        ) / det

  b(2,2) = +( &
        + a(1,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(1,3) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(1,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
        ) / det

  b(3,2) = -( &
        + a(1,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        + a(1,2) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(1,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(4,2) = +( &
        + a(1,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        + a(1,2) * ( a(3,3) * a(4,1) - a(3,1) * a(4,3) ) &
        + a(1,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(1,3) = +( &
        + a(1,2) * ( a(2,3) * a(4,4) - a(2,4) * a(4,3) ) &
        + a(1,3) * ( a(2,4) * a(4,2) - a(2,2) * a(4,4) ) &
        + a(1,4) * ( a(2,2) * a(4,3) - a(2,3) * a(4,2) ) &
        ) / det

  b(2,3) = -( &
        + a(1,1) * ( a(2,3) * a(4,4) - a(2,4) * a(4,3) ) &
        + a(1,3) * ( a(2,4) * a(4,1) - a(2,1) * a(4,4) ) &
        + a(1,4) * ( a(2,1) * a(4,3) - a(2,3) * a(4,1) ) &
        ) / det

  b(3,3) = +( &
        + a(1,1) * ( a(2,2) * a(4,4) - a(2,4) * a(4,2) ) &
        + a(1,2) * ( a(2,4) * a(4,1) - a(2,1) * a(4,4) ) &
        + a(1,4) * ( a(2,1) * a(4,2) - a(2,2) * a(4,1) ) &
        ) / det

  b(4,3) = -( &
        + a(1,1) * ( a(2,2) * a(4,3) - a(2,3) * a(4,2) ) &
        + a(1,2) * ( a(2,3) * a(4,1) - a(2,1) * a(4,3) ) &
        + a(1,3) * ( a(2,1) * a(4,2) - a(2,2) * a(4,1) ) &
        ) / det

  b(1,4) = -( &
        + a(1,2) * ( a(2,3) * a(3,4) - a(2,4) * a(3,3) ) &
        + a(1,3) * ( a(2,4) * a(3,2) - a(2,2) * a(3,4) ) &
        + a(1,4) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
        ) / det

  b(2,4) = +( &
        + a(1,1) * ( a(2,3) * a(3,4) - a(2,4) * a(3,3) ) &
        + a(1,3) * ( a(2,4) * a(3,1) - a(2,1) * a(3,4) ) &
        + a(1,4) * ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) &
        ) / det

  b(3,4) = -( &
        + a(1,1) * ( a(2,2) * a(3,4) - a(2,4) * a(3,2) ) &
        + a(1,2) * ( a(2,4) * a(3,1) - a(2,1) * a(3,4) ) &
        + a(1,4) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) &
        ) / det

  b(4,4) = +( &
        + a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
        + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
        + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) &
        ) / det

  return
end
subroutine r16mat_jac ( m, n, eps, fx, x, fprime )

!*****************************************************************************80
!
!! R16MAT_JAC estimates a dense jacobian matrix of the function FX.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    FPRIME(I,J) = d F(I) / d X(J).
!
!    The jacobian is assumed to be dense, and the LINPACK/LAPACK
!    double precision general matrix storage mode ("DGE") is used.
!
!    Forward differences are used, requiring N+1 function evaluations.
!
!    Values of EPS have typically been chosen between
!    sqrt ( EPSMCH ) and sqrt ( sqrt ( EPSMCH ) ) where EPSMCH is the
!    machine tolerance.
!
!    If EPS is too small, then F(X+EPS) will be the same as
!    F(X), and the jacobian will be full of zero entries.
!
!    If EPS is too large, the finite difference estimate will
!    be inaccurate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 16 ) EPS, a tolerance to be used for shifting the
!    X values during the finite differencing.  No single value
!    of EPS will be reliable for all vectors X and functions FX.
!
!    Input, external FX, the name of the user written
!    routine which evaluates the function at a given point X.
!
!    FX should have the form:
!
!      subroutine fx ( m, n, x, f )
!      integer m
!      integer n
!      real ( kind = 16 ) f(m)
!      real ( kind = 16 ) x(n)
!      f(1:m) = ...
!      return
!      end
!
!    Input, real ( kind = 16 ) X(N), the point where the jacobian
!    is to be estimated.
!
!    Output, real ( kind = 16 ) FPRIME(M,N), the M by N estimated jacobian
!    matrix.
!

  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) del
  real ( kind = 16 ) eps
  real ( kind = 16 ) fprime(m,n)
  external fx
  integer ( kind = 4 ) j
  real ( kind = 16 ) x(n)
  real ( kind = 16 ) xsave
  real ( kind = 16 ) work1(m)
  real ( kind = 16 ) work2(m)
!
!  Evaluate the function at the base point, X.
!
  call fx ( m, n, x, work2 )
!
!  Now, one by one, vary each component J of the base point X, and
!  estimate DF(I)/DX(J) = ( F(X+) - F(X) )/ DEL.
!
  do j = 1, n

    xsave = x(j)
    del = eps * ( 1.0_16 + abs ( x(j) ) )
    x(j) = x(j) + del
    call fx ( m, n, x, work1 )
    x(j) = xsave
    fprime(1:m,j) = ( work1(1:m) - work2(1:m) ) / del

  end do

  return
end
subroutine r16mat_l_inverse ( n, a, b )

!*****************************************************************************80
!
!! R16MAT_L_INVERSE inverts a lower triangular R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    A lower triangular matrix is a matrix whose only nonzero entries
!    occur on or below the diagonal.
!
!    The inverse of a lower triangular matrix is a lower triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 16 ) A(N,N), the lower triangular matrix.
!
!    Output, real ( kind = 16 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n

    do i = 1, n

      if ( i < j ) then
        b(i,j) = 0.0_16
      else if ( j == i ) then
        b(i,j) = 1.0_16 / a(i,j)
      else
        b(i,j) = - dot_product ( a(i,1:i-1), b(1:i-1,j) ) / a(i,i)
      end if

    end do
  end do

  return
end
subroutine r16mat_l_print ( m, n, a, title )

!*****************************************************************************80
!
!! R16MAT_L_PRINT prints a lower triangular R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Example:
!
!    M = 5, N = 5
!    A = (/ 11, 21, 31, 41, 51, 22, 32, 42, 52, 33, 43, 53, 44, 54, 55 /)
!
!    11
!    21 22
!    31 32 33
!    41 42 43 44
!    51 52 53 54 55
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(*), the M by N matrix.  Only the lower
!    triangular elements are stored, in column major order.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 16 ) a(*)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) indx(10)
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  integer   ( kind = 4 ) jmax
  integer   ( kind = 4 ) nn
  integer   ( kind = 4 ) size
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  jmax = min ( n, m )

  if ( m <= n ) then
    size = ( m * ( m + 1 ) ) / 2
  else if ( n < m ) then
    size = ( n * ( n + 1 ) ) / 2 + ( m - n ) * n
  end if

  if ( all ( a(1:size) == aint ( a(1:size) ) ) ) then

    nn = 10

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(a8,10i8)' ) '  Col   ', ( j, j = jlo, jhi )
      write ( *, '(a6)' ) '  Row '
      do i = jlo, m
        jhi = min ( jlo + nn - 1, i, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2
        end do
        write ( *, '(i8,10i8)' ) i, int ( a(indx(1:jhi+1-jlo)) )
      end do
    end do

  else if ( maxval ( abs ( a(1:size) ) ) < 1000000.0_16 ) then

    nn = 5

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(8x,5(i8,6x))' ) ( j, j = jlo, jhi )
      write ( *, '(a)' ) ' '
      do i = jlo, m
        jhi = min ( jlo + nn - 1, i, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2
        end do
        write ( *, '(i8,5f14.6)' ) i, a(indx(1:jhi+1-jlo))
      end do
    end do

  else

    nn = 5

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(8x,5(i8,6x))' ) ( j, j = jlo, jhi )
      write ( *, '(a)' ) ' '
      do i = jlo, m
        jhi = min ( jlo + nn - 1, i, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2
        end do
        write ( *, '(i8,5g14.6)' ) i, a(indx(1:jhi+1-jlo))
      end do
    end do

  end if

  return
end
subroutine r16mat_l_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R16MAT_L_SOLVE solves a lower triangular linear system.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 16 ) A(N,N), the N by N lower triangular matrix.
!
!    Input, real ( kind = 16 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 16 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 16 ) x(n)
!
!  Solve L * x = b.
!
  do i = 1, n
    x(i) = ( b(i) - dot_product ( a(i,1:i-1), x(1:i-1) ) ) / a(i,i)
  end do

  return
end
subroutine r16mat_l1_inverse ( n, a, b )

!*****************************************************************************80
!
!! R16MAT_L1_INVERSE inverts a unit lower triangular R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    A unit lower triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's above the main diagonal.
!
!    The inverse of a unit lower triangular matrix is also
!    a unit lower triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call r16mat_l1_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 16 ) A(N,N), the unit lower triangular matrix.
!
!    Output, real ( kind = 16 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n

    do j = 1, n

      if ( i < j ) then
        b(i,j) = 0.0_16
      else if ( j == i ) then
        b(i,j) = 1.0_16
      else
        b(i,j) = - dot_product ( a(i,1:i-1), b(1:i-1,j) )
      end if

    end do
  end do

  return
end
subroutine r16mat_lt_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R16MAT_LT_SOLVE solves a transposed lower triangular linear system.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    Given the lower triangular matrix A, the linear system to be solved is:
!
!      A' * x = b
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns
!    of the matrix.
!
!    Input, real ( kind = 16 ) A(N,N), the N by N lower triangular matrix.
!
!    Input, real ( kind = 16 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 16 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 16 ) x(n)
!
!  Solve L'*x = b.
!
  do i = n, 1, -1
    x(i) = ( b(i) - dot_product ( x(i+1:n), a(i+1:n,i) ) ) / a(i,i)
  end do

  return
end
subroutine r16mat_lu ( m, n, a, l, p, u )

!*****************************************************************************80
!
!! R16MAT_LU computes the LU factorization of a rectangular R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The routine is given an M by N matrix A, and produces
!
!      L, an M by M unit lower triangular matrix,
!      U, an M by N upper triangular matrix, and
!      P, an M by M permutation matrix P,
!
!    so that
!
!      A = P' * L * U.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the M by N matrix to be factored.
!
!    Output, real ( kind = 16 ) L(M,M), the M by M unit lower triangular factor.
!
!    Output, real ( kind = 16 ) P(M,M), the M by M permutation matrix.
!
!    Output, real ( kind = 16 ) U(M,N), the M by N upper triangular factor.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipiv
  integer ( kind = 4 ) j
  real ( kind = 16 ) l(m,m)
  real ( kind = 16 ) p(m,m)
  real ( kind = 16 ) pivot
  real ( kind = 16 ) u(m,n)

!  Initialize:
!
!    U:=A
!    L:=Identity
!    P:=Identity
!
  u(1:m,1:n) = a(1:m,1:n)

  call r16mat_identity ( m, l )

  p(1:m,1:m) = l(1:m,1:m)
!
!  On step J, find the pivot row, IPIV, and the pivot value PIVOT.
!
  do j = 1, min ( m - 1, n )

    pivot = 0.0_16
    ipiv = 0

    do i = j, m

      if ( pivot < abs ( u(i,j) ) ) then
        pivot = abs ( u(i,j) )
        ipiv = i
      end if

    end do
!
!  Unless IPIV is zero, swap rows J and IPIV.
!
    if ( ipiv /= 0 ) then

      call r16row_swap ( m, n, u, j, ipiv )

      call r16row_swap ( m, m, l, j, ipiv )

      call r16row_swap ( m, m, p, j, ipiv )
!
!  Zero out the entries in column J, from row J+1 to M.
!
      do i = j + 1, m

        if ( u(i,j) /= 0.0_16 ) then

          l(i,j) = u(i,j) / u(j,j)

          u(i,j) = 0.0_16

          u(i,j+1:n) = u(i,j+1:n) - l(i,j) * u(j,j+1:n)

        end if

      end do

    end if

  end do

  return
end
function r16mat_max ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_MAX returns the maximum entry of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the M by N matrix.
!
!    Output, real ( kind = 16 ) R16MAT_MAX, the maximum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  real ( kind = 16 ) r16mat_max

  r16mat_max = maxval ( a(1:m,1:n) )

  return
end
subroutine r16mat_max_index ( m, n, a, i, j )

!*****************************************************************************80
!
!! R16MAT_MAX_INDEX returns the location of the maximum entry of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the maximum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj

  i = - 1
  j = - 1

  do jj = 1, n
    do ii = 1, m
      if ( ii == 1 .and. jj == 1 ) then
        i = ii
        j = jj
      else if ( a(i,j) < a(ii,jj) ) then
        i = ii
        j = jj
      end if
    end do
  end do

  return
end
function r16mat_maxcol_minrow ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_MAXCOL_MINROW gets the maximum column minimum row of an M by N R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    R16MAT_MAXCOL_MINROW = max ( 1 <= I <= N ) ( min ( 1 <= J <= M ) A(I,J) )
!
!    For a given matrix, R16MAT_MAXCOL_MINROW <= R16MAT_MINROW_MAXCOL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix.
!
!    Output, real ( kind = 16 ) R16MAT_MAXCOL_MINROW, the maximum column
!    minimum row entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) i
  real ( kind = 16 ) r16mat_maxcol_minrow
  real ( kind = 16 ) r16mat_minrow

  r16mat_maxcol_minrow = 0.0_16

  do i = 1, m

    r16mat_minrow = minval ( a(i,1:n) )

    if ( i == 1 ) then
      r16mat_maxcol_minrow = r16mat_minrow
    else
      r16mat_maxcol_minrow = max ( r16mat_maxcol_minrow, r16mat_minrow )
    end if

  end do

  return
end
function r16mat_maxrow_mincol ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_MAXROW_MINCOL gets the maximum row minimum column of an M by N R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    R16MAT_MAXROW_MINCOL = max ( 1 <= J <= N ) ( min ( 1 <= I <= M ) A(I,J) )
!
!    For a given matrix, R16MAT_MAXROW_MINCOL <= R16MAT_MINCOL_MAXROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix.
!
!    Output, real ( kind = 16 ) R16MAT_MAXROW_MINCOL, the maximum row
!    minimum column entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) j
  real ( kind = 16 ) r16mat_maxrow_mincol
  real ( kind = 16 ) r16mat_mincol

  r16mat_maxrow_mincol = 0.0_16

  do j = 1, n

    r16mat_mincol = minval ( a(1:m,j) )

    if ( j == 1 ) then
      r16mat_maxrow_mincol = r16mat_mincol
    else
      r16mat_maxrow_mincol = max ( r16mat_maxrow_mincol, r16mat_mincol )
    end if

  end do

  return
end
function r16mat_min ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_MIN returns the minimum entry of an M by N R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix.
!
!    Output, real ( kind = 16 ) R16MAT_MIN, the minimum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  real ( kind = 16 ) r16mat_min

  r16mat_min = minval ( a(1:m,1:n) )

  return
end
subroutine r16mat_min_index ( m, n, a, i, j )

!*****************************************************************************80
!
!! R16MAT_MIN_INDEX returns the location of the minimum entry of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the minimum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj

  i = -1
  j = -1

  do jj = 1, n
    do ii = 1, m
      if ( ii == 1 .and. jj == 1 ) then
        i = ii
        j = jj
      else if ( a(ii,jj) < a(i,j) ) then
        i = ii
        j = jj
      end if
    end do
  end do

  return
end
function r16mat_mincol_maxrow ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_MINCOL_MAXROW gets the minimum column maximum row of an M by N R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    R16MAT_MINCOL_MAXROW = min ( 1 <= I <= N ) ( max ( 1 <= J <= M ) A(I,J) )
!
!    For a given matrix, R16MAT_MAXROW_MINCOL <= R16MAT_MINCOL_MAXROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix.
!
!    Output, real ( kind = 16 ) R16MAT_MINCOL_MAXROW, the minimum column
!    maximum row entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) i
  real ( kind = 16 ) r16mat_mincol_maxrow
  real ( kind = 16 ) r16mat_maxrow

  r16mat_mincol_maxrow = 0.0_16

  do i = 1, m

    r16mat_maxrow = maxval ( a(i,1:n) )

    if ( i == 1 ) then
      r16mat_mincol_maxrow = r16mat_maxrow
    else
      r16mat_mincol_maxrow = min ( r16mat_mincol_maxrow, r16mat_maxrow )
    end if

  end do

  return
end
function r16mat_minrow_maxcol ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_MINROW_MAXCOL gets the minimum row maximum column of an M by N R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    R16MAT_MINROW_MAXCOL = min ( 1 <= J <= N ) ( max ( 1 <= I <= M ) A(I,J) )
!
!    For a given matrix, R16MAT_MAXCOL_MINROW <= R16MAT_MINROW_MAXCOL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix.
!
!    Output, real ( kind = 16 ) R16MAT_MINROW_MAXCOL, the minimum row
!    maximum column entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) j
  real ( kind = 16 ) r16mat_minrow_maxcol
  real ( kind = 16 ) r16mat_maxcol

  r16mat_minrow_maxcol = 0.0_16

  do j = 1, n

    r16mat_maxcol = maxval ( a(1:m,j) )

    if ( j == 1 ) then
      r16mat_minrow_maxcol = r16mat_maxcol
    else
      r16mat_minrow_maxcol = min ( r16mat_minrow_maxcol, r16mat_maxcol )
    end if

  end do

  return
end
subroutine r16mat_mm ( n1, n2, n3, a, b, c )

!*****************************************************************************80
!
!! R16MAT_MM multiplies two R16MAT's.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    In FORTRAN90, this operation is more efficiently done by the
!    command:
!
!      C(1:N1,1:N3) = MATMUL ( A(1:N1,1;N2), B(1:N2,1:N3) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the order of the matrices.
!
!    Input, real ( kind = 16 ) A(N1,N2), B(N2,N3), the matrices to multiply.
!
!    Output, real ( kind = 16 ) C(N1,N3), the product matrix C = A * B.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  real ( kind = 16 ) a(n1,n2)
  real ( kind = 16 ) b(n2,n3)
  real ( kind = 16 ) c(n1,n3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  do i = 1, n1
    do j = 1, n3
      c(i,j) = 0.0_16
      do k = 1, n2
        c(i,j) = c(i,j) + a(i,k) * b(k,j)
      end do
    end do
  end do

  return
end
subroutine r16mat_mtv ( m, n, a, x, y )

!*****************************************************************************80
!
!! R16MAT_MTV multiplies a transposed matrix times a vector
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix.
!
!    Input, real ( kind = 16 ) A(M,N), the M by N matrix.
!
!    Input, real ( kind = 16 ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = 16 ) Y(N), the product A'*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  real ( kind = 16 ) x(m)
  real ( kind = 16 ) y(n)

  y(1:n) = matmul ( transpose ( a(1:m,1:n) ), x(1:m) )

  return
end
subroutine r16mat_mv ( m, n, a, x, y )

!*****************************************************************************80
!
!! R16MAT_MV multiplies a matrix times a vector.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    In FORTRAN90, this operation can be more efficiently carried
!    out by the command
!
!      Y(1:M) = MATMUL ( A(1:M,1:N), X(1:N) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix.
!
!    Input, real ( kind = 16 ) A(M,N), the M by N matrix.
!
!    Input, real ( kind = 16 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 16 ) Y(M), the product A*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 16 ) x(n)
  real ( kind = 16 ) y(m)

  do i = 1, m
    y(i) = 0.0_16
    do j = 1, n
      y(i) = y(i) + a(i,j) * x(j)
    end do
  end do

  return
end
subroutine r16mat_nint ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_NINT rounds the entries of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!
!    Input/output, real ( kind = 16 ) A(M,N), the matrix to be NINT'ed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)

  a(1:m,1:n) = real ( nint ( a(1:m,1:n) ), kind = 16 )

  return
end
function r16mat_norm_eis ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_NORM_EIS returns the EISPACK norm of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The EISPACK norm is defined as:
!
!      R16MAT_NORM_EIS =
!        sum ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix whose EISPACK norm is desired.
!
!    Output, real ( kind = 16 ) R16MAT_NORM_EIS, the EISPACK norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  real ( kind = 16 ) r16mat_norm_eis

  r16mat_norm_eis = sum ( abs ( a(1:m,1:n) ) )

  return
end
function r16mat_norm_fro ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_NORM_FRO returns the Frobenius norm of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The Frobenius norm is defined as
!
!      R16MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J) * A(I,J) )
!
!    The matrix Frobenius norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r16vec_norm_l2 ( A * x ) <= r16mat_norm_fro ( A ) * r16vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix whose Frobenius
!    norm is desired.
!
!    Output, real ( kind = 16 ) R16MAT_NORM_FRO, the Frobenius norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  real ( kind = 16 ) r16mat_norm_fro

  r16mat_norm_fro = sqrt ( sum ( a(1:m,1:n)**2 ) )

  return
end
function r16mat_norm_l1 ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_NORM_L1 returns the matrix L1 norm of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The matrix L1 norm is defined as:
!
!      R16MAT_NORM_L1 = max ( 1 <= J <= N )
!        sum ( 1 <= I <= M ) abs ( A(I,J) ).
!
!    The matrix L1 norm is derived from the vector L1 norm, and
!    satisifies:
!
!      r16vec_norm_l1 ( A * x ) <= r16mat_norm_l1 ( A ) * r16vec_norm_l1 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix whose L1 norm is desired.
!
!    Output, real ( kind = 16 ) R16MAT_NORM_L1, the L1 norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) j
  real ( kind = 16 ) r16mat_norm_l1

  r16mat_norm_l1 = 0.0_16

  do j = 1, n
    r16mat_norm_l1 = max ( r16mat_norm_l1, sum ( abs ( a(1:m,j) ) ) )
  end do

  return
end
function r16mat_norm_l2 ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_NORM_L2 returns the matrix L2 norm of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The matrix L2 norm is defined as:
!
!      R16MAT_NORM_L2 = sqrt ( max ( 1 <= I <= M ) LAMBDA(I) )
!
!    where LAMBDA contains the eigenvalues of A * A'.
!
!    The matrix L2 norm is derived from the vector L2 norm, and
!    satisifies:
!
!      r16vec_norm_l2 ( A * x ) <= r16mat_norm_l2 ( A ) * r16vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix whose L2 norm is desired.
!
!    Output, real ( kind = 16 ) R16MAT_NORM_L2, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  real ( kind = 16 ) b(m,m)
  real ( kind = 16 ) diag(m)
  real ( kind = 16 ) r16mat_norm_l2
!
!  Compute B = A * A'.
!
  b(1:m,1:m) = matmul ( a(1:m,1:n), transpose ( a(1:m,1:n) ) )
!
!  Diagonalize B.
!
  call r16mat_symm_jacobi ( m, b )
!
!  Find the maximum eigenvalue, and take its square root.
!
  call r16mat_diag_get_vector ( m, b, diag )

  r16mat_norm_l2 = sqrt ( maxval ( diag(1:m) ) )

  return
end
function r16mat_norm_li ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_NORM_LI returns the matrix L-oo norm of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The matrix L-oo norm is defined as:
!
!      R16MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
!
!    The matrix L-oo norm is derived from the vector L-oo norm,
!    and satisifies:
!
!      r16vec_norm_li ( A * x ) <= r16mat_norm_li ( A ) * r16vec_norm_li ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix whose L-oo
!    norm is desired.
!
!    Output, real ( kind = 16 ) R16MAT_NORM_LI, the L-oo norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) i
  real ( kind = 16 ) r16mat_norm_li

  r16mat_norm_li = 0.0_16

  do i = 1, m
    r16mat_norm_li = max ( r16mat_norm_li, sum ( abs ( a(i,1:n) ) ) )
  end do

  return
end
subroutine r16mat_nullspace ( m, n, a, nullspace_size, nullspace )

!*****************************************************************************80
!
!! R16MAT_NULLSPACE computes the nullspace of a matrix.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
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
!    08 August 2010
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
!    Input, real ( kind = 16 ) A(M,N), the matrix to be analyzed.
!
!    Input, integer ( kind = 4 ) NULLSPACE_SIZE, the size of the nullspace.
!
!    Output, real ( kind = 16 ) NULLSPACE(N,NULLSPACE_SIZE), vectors that
!    span the nullspace.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nullspace_size

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) col(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  real ( kind = 16 ) nullspace(n,nullspace_size)
  integer ( kind = 4 ) row(m)
  real ( kind = 16 ) rref(m,n)
!
!  Make a copy of A.
!
  rref(1:m,1:n) = a(1:m,1:n)
!
!  Get the reduced row echelon form of A.
!
  call r16mat_rref ( m, n, rref )
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
      if ( rref(i,j) == 1.0_16 ) then
        row(i) = j
        col(j) = j
        exit
      end if
    end do
  end do

  nullspace(1:n,1:nullspace_size) = 0.0_16

  j2 = 0
!
!  If column J does not contain a leading 1, then it contains
!  information about a null vector.
!
  do j = 1, n

    if ( col(j) < 0 ) then

      j2 = j2 + 1

      do i = 1, m
        if ( rref(i,j) /= 0.0_16 ) then
          i2 = row(i)
          nullspace(i2,j2) = - rref(i,j)
        end if
      end do

      nullspace(j,j2) = 1.0_16

    end if

  end do

  return
end
subroutine r16mat_nullspace_size ( m, n, a, nullspace_size )

!*****************************************************************************80
!
!! R16MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
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
!    08 August 2010
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
!    Input, real ( kind = 16 ) A(M,N), the matrix to be analyzed.
!
!    Output, integer ( kind = 4 ) NULLSPACE_SIZE, the estimated size
!    of the nullspace.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) leading
  integer ( kind = 4 ) nullspace_size
  real ( kind = 16 ) rref(m,n)
!
!  Get the reduced row echelon form of A.
!
  rref(1:m,1:n) = a(1:m,1:n)

  call r16mat_rref ( m, n, rref )
!
!  Count the leading 1's in A.
!
  leading = 0
  do i = 1, m
    do j = 1, n
      if ( rref(i,j) == 1.0_16 ) then
        leading = leading + 1
        exit
      end if
    end do
  end do

  nullspace_size = n - leading

  return
end
subroutine r16mat_orth_uniform ( n, seed, a )

!*****************************************************************************80
!
!! R16MAT_ORTH_UNIFORM returns a random orthogonal R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
!    National Academy of Sciences of Belarus, for convincingly
!    pointing out the severe deficiencies of an earlier version of
!    this routine.
!
!    Essentially, the computation involves saving the Q factor of the
!    QR factorization of a matrix whose entries are normally distributed.
!    However, it is only necessary to generate this matrix a column at
!    a time, since it can be shown that when it comes time to annihilate
!    the subdiagonal elements of column K, these (transformed) elements of
!    column K are still normally distributed random values.  Hence, there
!    is no need to generate them at the beginning of the process and
!    transform them K-1 times.
!
!    For computational efficiency, the individual Householder transformations
!    could be saved, as recommended in the reference, instead of being
!    accumulated into an explicit matrix format.
!
!  Properties:
!
!    The inverse of A is equal to A'.
!
!    A * A'  = A' * A = I.
!
!    Columns and rows of A have unit Euclidean norm.
!
!    Distinct pairs of columns of A are orthogonal.
!
!    Distinct pairs of rows of A are orthogonal.
!
!    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.
!
!    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.
!
!    The determinant of A is +1 or -1.
!
!    All the eigenvalues of A have modulus 1.
!
!    All singular values of A are 1.
!
!    All entries of A are between -1 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pete Stewart,
!    Efficient Generation of Random Orthogonal Matrices With an Application
!    to Condition Estimators,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 3, June 1980, pages 403-409.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 16 ) A(N,N), the orthogonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 16 ) r16_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 16 ) v(n)
  real ( kind = 16 ) x(n)
!
!  Start with A = the identity matrix.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 1.0_16
      else
        a(i,j) = 0.0_16
      end if
    end do
  end do
!
!  Now behave as though we were computing the QR factorization of
!  some other random matrix.  Generate the N elements of the first column,
!  compute the Householder matrix H1 that annihilates the subdiagonal elements,
!  and set A := A * H1' = A * H.
!
!  On the second step, generate the lower N-1 elements of the second column,
!  compute the Householder matrix H2 that annihilates them,
!  and set A := A * H2' = A * H2 = H1 * H2.
!
!  On the N-1 step, generate the lower 2 elements of column N-1,
!  compute the Householder matrix HN-1 that annihilates them, and
!  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
!  This is our random orthogonal matrix.
!
  do j = 1, n-1
!
!  Set the vector that represents the J-th column to be annihilated.
!
    x(1:j-1) = 0.0_16

    do i = j, n
      x(i) = r16_normal_01 ( seed )
    end do
!
!  Compute the vector V that defines a Householder transformation matrix
!  H(V) that annihilates the subdiagonal elements of X.
!
    call r16vec_house_column ( n, x, j, v )
!
!  Postmultiply the matrix A by H'(V) = H(V).
!
    call r16mat_house_axh ( n, a, v, a )

  end do

  return
end
subroutine r16mat_plot ( m, n, a, title )

!*****************************************************************************80
!
!! R16MAT_PLOT "plots" an R16MAT, with an optional title.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 16 ) a(m,n)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character              r16mat_plot_symbol
  character ( len = 70 ) string
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do jlo = 1, n, 70
    jhi = min ( jlo + 70-1, n )
    write ( *, '(a)' ) ' '
    write ( *, '(8x,2x,70i1)' ) ( mod ( j, 10 ), j = jlo, jhi )
    write ( *, '(a)' ) ' '

    do i = 1, m
      do j = jlo, jhi
        string(j+1-jlo:j+1-jlo) = r16mat_plot_symbol ( a(i,j) )
      end do
      write ( *, '(i8,2x,a)' ) i, string(1:jhi+1-jlo)
    end do
  end do

  return
end
function r16mat_plot_symbol ( r )

!*****************************************************************************80
!
!! R16MAT_PLOT_SYMBOL returns a symbol for an element of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) R, a value whose symbol is desired.
!
!    Output, character R16MAT_PLOT_SYMBOL, is
!    '-' if R is negative,
!    '0' if R is zero,
!    '+' if R is positive.
!
  implicit none

  character              r16mat_plot_symbol
  real      ( kind = 16 ) r

  if ( r < 0.0_16 ) then
    r16mat_plot_symbol = '-'
  else if ( r == 0.0_16 ) then
    r16mat_plot_symbol = '0'
  else if ( 0.0_16 < r ) then
    r16mat_plot_symbol = '+'
  end if

  return
end
subroutine r16mat_poly_char ( n, a, p )

!*****************************************************************************80
!
!! R16MAT_POLY_CHAR computes the characteristic polynomial of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, real ( kind = 16 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 16 ) P(0:N), the coefficients of the characteristic
!    polynomial of A.  P(N) contains the coefficient of X**N
!    (which will be 1), P(I) contains the coefficient of X**I,
!    and P(0) contains the constant term.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real ( kind = 16 ) p(0:n)
  real ( kind = 16 ) r16mat_trace
  real ( kind = 16 ) trace
  real ( kind = 16 ) work1(n,n)
  real ( kind = 16 ) work2(n,n)
!
!  Initialize WORK1 to the identity matrix.
!
  call r16mat_identity ( n, work1 )

  p(n) = 1.0_16

  do order = n - 1, 0, -1
!
!  Work2 = A * WORK1.
!
    work2(1:n,1:n) = matmul ( a(1:n,1:n), work1(1:n,1:n) )
!
!  Take the trace.
!
    trace = r16mat_trace ( n, work2 )
!
!  P(ORDER) = -Trace ( WORK2 ) / ( N - ORDER )
!
    p(order) = - trace / real ( n - order, kind = 16 )
!
!  WORK1 := WORK2 + P(ORDER) * Identity.
!
    work1(1:n,1:n) = work2(1:n,1:n)

    do i = 1, n
      work1(i,i) = work1(i,i) + p(order)
    end do

  end do

  return
end
subroutine r16mat_power ( n, a, npow, b )

!*****************************************************************************80
!
!! R16MAT_POWER computes a nonnegative power of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The algorithm is:
!
!      B = I
!      do NPOW times:
!        B = A * B
!      end
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 16 ) A(N,N), the matrix to be raised to a power.
!
!    Input, integer ( kind = 4 ) NPOW, the power to which A is to be raised.
!    NPOW must be nonnegative.
!
!    Output, real ( kind = 16 ) B(N,N), the value of A**NPOW.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) b(n,n)
  integer ( kind = 4 ) ipow
  integer ( kind = 4 ) npow

  if ( npow < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R16MAT_POWER - Fatal error!'
    write ( *, '(a)' ) '  Input value of NPOW < 0.'
    write ( *, '(a,i8)' ) '  NPOW = ', npow
    stop
  end if

  call r16mat_identity ( n, b )

  do ipow = 1, npow
    b(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )
  end do

  return
end
subroutine r16mat_power_method ( n, a, r, v )

!*****************************************************************************80
!
!! R16MAT_POWER_METHOD applies the power method to an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    If the power method has not converged, then calling the routine
!    again immediately with the output from the previous call will
!    continue the iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 16 ) A(N,N), the matrix.
!
!    Output, real ( kind = 16 ) R, the estimated eigenvalue.
!
!    Input/output, real ( kind = 16 ) V(N), on input, an estimate
!    for the eigenvector.  On output, an improved estimate for the
!    eigenvector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) av(n)
  real ( kind = 16 ) eps
  integer ( kind = 4 ) it
  real ( kind = 16 ), parameter :: it_eps = 0.0001_16
  integer ( kind = 4 ), parameter :: it_max = 100
  integer ( kind = 4 ), parameter :: it_min = 10
  integer ( kind = 4 ) j
  real ( kind = 16 ) r
  real ( kind = 16 ) r2
  real ( kind = 16 ) r_old
  real ( kind = 16 ) v(n)

  eps = sqrt ( epsilon ( 1.0_16 ) )

  r = sqrt ( sum ( v(1:n)**2 ) )

  if ( r == 0.0_16 ) then
    v(1:n) = 1.0_16
    r = sqrt ( real ( n, kind = 16 ) )
  end if

  v(1:n) = v(1:n) / r

  do it = 1, it_max

    av(1:n) = matmul ( a(1:n,1:n), v(1:n) )

    r_old = r
    r = sqrt ( sum ( av(1:n)**2 ) )

    if ( it_min < it ) then
      if ( abs ( r - r_old ) <= it_eps * ( 1.0_16 + abs ( r ) ) ) then
        exit
      end if
    end if

    v(1:n) = av(1:n)

    if ( r /= 0.0_16 ) then
      v(1:n) = v(1:n) / r
    end if
!
!  Perturb V a bit, to avoid cases where the initial guess is exactly
!  the eigenvector of a smaller eigenvalue.
!
    if ( it < it_max / 2 ) then
      j = 1 + mod ( it - 1, n )
      v(j) = v(j) + eps * ( 1.0_16 + abs ( v(j) ) )
      r2 = sqrt ( sum ( v(1:n)**2 ) )
      v(1:n) = v(1:n) / r2
    end if

  end do

  return
end
subroutine r16mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R16MAT_PRINT prints an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 16 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  character ( len = * )  title

  call r16mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r16mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R16MAT_PRINT_SOME prints some of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 16 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 16 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r16mat_print2 ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_PRINT2 prints an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input, real ( kind = 16 ) A(M,N), the M by N matrix to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  real ( kind = 16 ) amax
  real ( kind = 16 ) amin
  integer ( kind = 4 ) i
  character ( len = 10 ) iform
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  logical integ
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) lmax
  integer ( kind = 4 ) npline
  real ( kind = 16 ) r16_log_10
!
!  Check if all entries are integral.
!
  integ = .true.

  do i = 1, m
    do j = 1, n

      if ( integ ) then
        if ( a(i,j) /= real ( int ( a(i,j) ), kind = 16 ) ) then
          integ = .false.
        end if
      end if

    end do
  end do
!
!  Find the maximum and minimum entries.
!
  amax = maxval ( a(1:m,1:n) )
  amin = minval ( a(1:m,1:n) )
!
!  Use the information about the maximum size of an entry to
!  compute an intelligent format for use with integer entries.
!
!  Later, we might also do this for real matrices.
!
  lmax = int ( r16_log_10 ( amax ) )

  if ( integ ) then
    npline = 79 / ( lmax + 3 )
    write ( iform, '(''('',i2,''I'',i2,'')'')' ) npline, lmax+3
  else
    npline = 5
    iform = ' '
  end if
!
!  Print a scalar quantity.
!
  if ( m == 1 .and. n == 1 ) then

    if ( integ ) then
      write ( *, iform ) int ( a(1,1) )
    else
      write ( *, '(2x,g14.6)' ) a(1,1)
    end if
!
!  Column vector of length M,
!
  else if ( n == 1 ) then

    do ilo = 1, m, npline

      ihi = min ( ilo+npline-1, m )

      if ( integ ) then
        write ( *, iform ) ( int ( a(i,1) ), i = ilo, ihi )
      else
        write ( *, '(2x,5g14.6)' ) a(ilo:ihi,1)
      end if

    end do
!
!  Row vector of length N,
!
  else if ( m == 1 ) then

    do jlo = 1, n, npline

      jhi = min ( jlo+npline-1, n )

      if ( integ ) then
        write ( *, iform ) int ( a(1,jlo:jhi) )
      else
        write ( *, '(2x,5g14.6)' ) a(1,jlo:jhi)
      end if

    end do
!
!  M by N Array
!
  else

    do jlo = 1, n, npline

      jhi = min ( jlo+npline-1, n )

      if ( npline < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8,a,i8)' ) 'Matrix columns ', jlo, ' to ', jhi
        write ( *, '(a)' ) ' '
      end if

      do i = 1, m

        if ( integ ) then
          write ( *, iform ) int ( a(i,jlo:jhi) )
        else
          write ( *, '(2x,5g14.6)' ) a(i,jlo:jhi)
        end if

      end do
    end do

  end if

  return
end
subroutine r16mat_ref ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_REF computes the row echelon form of a matrix.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
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
!     1.0  3.0  0.0  2.0  6.0  3.0  1.0
!     0.0  0.0  0.0  1.0  2.0  4.5  1.5
!     0.0  0.0  0.0  0.0  0.0  1.0  0.3
!     0.0  0.0  0.0  0.0  0.0  0.0  0.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
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
!    Input/output, real ( kind = 16 ) A(M,N).  On input, the matrix to be
!    analyzed.  On output, the REF form of the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lead
  integer ( kind = 4 ) r
  real ( kind = 16 ) temp

  lead = 1

  do r = 1, m

    if ( n < lead ) then
      exit
    end if

    i = r

    do while ( a(i,lead) == 0.0_16 )

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

    do i = r + 1, m
      a(i,1:n) = a(i,1:n) - a(i,lead) * a(r,1:n)
    end do

    lead = lead + 1

  end do

  return
end
subroutine r16mat_rref ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_RREF computes the reduced row echelon form of a matrix.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
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
!    08 August 2010
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
!    Input/output, real ( kind = 16 ) A(M,N).  On input, the matrix to be
!    analyzed.  On output, the RREF form of the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lead
  integer ( kind = 4 ) r
  real ( kind = 16 ) temp

  lead = 1

  do r = 1, m

    if ( n < lead ) then
      exit
    end if

    i = r

    do while ( a(i,lead) == 0.0_16 )

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
subroutine r16mat_solve ( n, rhs_num, a, info )

!*****************************************************************************80
!
!! R16MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) RHS_NUM, the number of right hand sides.
!    RHS_NUM must be at least 0.
!
!    Input/output, real ( kind = 16 ) A(N,N+RHS_NUM), contains in rows and
!    columns 1 to N the coefficient matrix, and in columns N+1 through
!    N+rhs_num, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) rhs_num

  real ( kind = 16 ) a(n,n+rhs_num)
  real ( kind = 16 ) apivot
  real ( kind = 16 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot
  integer ( kind = 4 ) j
  real ( kind = 16 ) t(n+rhs_num)

  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j + 1, n
      if ( abs ( apivot ) < abs ( a(i,j) ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0_16 ) then
      info = j
      return
    end if
!
!  The pivot row moves into the J-th row.
!
    if ( ipivot /= j ) then
      t(       1:n+rhs_num) = a(ipivot,1:n+rhs_num)
      a(ipivot,1:n+rhs_num) = a(j,     1:n+rhs_num)
      a(j,     1:n+rhs_num) = t(       1:n+rhs_num)
    end if
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0_16
    a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then
        factor = a(i,j)
        a(i,j) = 0.0_16
        a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)
      end if

    end do

  end do

  return
end
subroutine r16mat_solve_2d ( a, b, det, x )

!*****************************************************************************80
!
!! R16MAT_SOLVE_2D solves a 2 by 2 linear system using Cramer's rule.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    If the determinant DET is returned as zero, then the matrix A is
!    singular, and does not have an inverse.  In that case, X is
!    returned as the zero vector.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A(2,2), the matrix.
!
!    Input, real ( kind = 16 ) B(2), the right hand side.
!
!    Output, real ( kind = 16 ) DET, the determinant of the matrix A.
!
!    Output, real ( kind = 16 ) X(2), the solution of the system,
!    if DET is nonzero.
!
  implicit none

  real ( kind = 16 ) a(2,2)
  real ( kind = 16 ) b(2)
  real ( kind = 16 ) det
  real ( kind = 16 ) x(2)
!
!  Compute the determinant.
!
  det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0_16 ) then
    x(1:2) = 0.0_16
    return
  end if
!
!  Compute the solution.
!
  x(1) = (  a(2,2) * b(1) - a(1,2) * b(2) ) / det
  x(2) = ( -a(2,1) * b(1) + a(1,1) * b(2) ) / det

  return
end
subroutine r16mat_solve_3d ( a, b, det, x )

!*****************************************************************************80
!
!! R16MAT_SOLVE_3D solves a 3 by 3 linear system using Cramer's rule.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    If the determinant DET is returned as zero, then the matrix A is
!    singular, and does not have an inverse.  In that case, X is
!    returned as the zero vector.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 16 ) A(3,3), the matrix.
!
!    Input, real ( kind = 16 ) B(3), the right hand side.
!
!    Output, real ( kind = 16 ) DET, the determinant of the matrix A.
!
!    Output, real ( kind = 16 ) X(3), the solution of the system,
!    if DET is nonzero.
!
  implicit none

  real ( kind = 16 ) a(3,3)
  real ( kind = 16 ) b(3)
  real ( kind = 16 ) det
  real ( kind = 16 ) x(3)
!
!  Compute the determinant.
!
  det =  a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
       + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
       + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0_16 ) then
    x(1:3) = 0.0_16
    return
  end if
!
!  Compute the solution.
!
  x(1) = (   ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) * b(1) &
           - ( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) * b(2) &
           + ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) * b(3) ) / det

  x(2) = ( - ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) * b(1) &
           + ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) * b(2) &
           - ( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) * b(3) ) / det

  x(3) = (   ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) * b(1) &
           - ( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) * b(2) &
           + ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) * b(3) ) / det

  return
end
subroutine r16mat_solve2 ( n, a, b, x, ierror )

!*****************************************************************************80
!
!! R16MAT_SOLVE2 computes the solution of an N by N linear system.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The linear system may be represented as
!
!      A*X = B
!
!    If the linear system is singular, but consistent, then the routine will
!    still produce a solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input/output, real ( kind = 16 ) A(N,N).
!    On input, A is the coefficient matrix to be inverted.
!    On output, A has been overwritten.
!
!    Input/output, real ( kind = 16 ) B(N).
!    On input, B is the right hand side of the system.
!    On output, B has been overwritten.
!
!    Output, real ( kind = 16 ) X(N), the solution of the linear system.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error detected.
!    1, consistent singularity.
!    2, inconsistent singularity.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) amax
  real ( kind = 16 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) ipiv(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 16 ) x(n)

  ierror = 0

  ipiv(1:n) = 0
  x(1:n) = 0.0_16
!
!  Process the matrix.
!
  do k = 1, n
!
!  In column K:
!    Seek the row IMAX with the properties that:
!      IMAX has not already been used as a pivot;
!      A(IMAX,K) is larger in magnitude than any other candidate.
!
    amax = 0.0_16
    imax = 0
    do i = 1, n
      if ( ipiv(i) == 0 ) then
        if ( amax < abs ( a(i,k) ) ) then
          imax = i
          amax = abs ( a(i,k) )
        end if
      end if
    end do
!
!  If you found a pivot row IMAX, then,
!    eliminate the K-th entry in all rows that have not been used for pivoting.
!
    if ( imax /= 0 ) then

      ipiv(imax) = k
      a(imax,k+1:n) = a(imax,k+1:n) / a(imax,k)
      b(imax) = b(imax) / a(imax,k)
      a(imax,k) = 1.0_16

      do i = 1, n

        if ( ipiv(i) == 0 ) then
          a(i,k+1:n) = a(i,k+1:n) - a(i,k) * a(imax,k+1:n)
          b(i) = b(i) - a(i,k) * b(imax)
          a(i,k) = 0.0_16
        end if

      end do

    end if

  end do
!
!  Now, every row with nonzero IPIV begins with a 1, and
!  all other rows are all zero.  Begin solution.
!
  do j = n, 1, -1

    imax = 0
    do k = 1, n
      if ( ipiv(k) == j ) then
        imax = k
      end if
    end do

    if ( imax == 0 ) then

      x(j) = 0.0_16

      if ( b(j) == 0.0_16 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R16MAT_SOLVE2 - Warning:'
        write ( *, '(a,i8)' ) '  Consistent singularity, equation = ', j
      else
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R16MAT_SOLVE2 - Error:'
        write ( *, '(a,i8)' ) '  Inconsistent singularity, equation = ', j
      end if

    else

      x(j) = b(imax)

      do i = 1, n
        if ( i /= imax ) then
          b(i) = b(i) - a(i,j) * x(j)
        end if
      end do

    end if

  end do

  return
end
subroutine r16mat_symm_eigen ( n, x, q, a )

!*****************************************************************************80
!
!! R16MAT_SYMM_EIGEN returns a symmetric matrix with given eigensystem.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The user must supply the desired eigenvalue vector, and the desired
!    eigenvector matrix.  The eigenvector matrix must be orthogonal.  A
!    suitable random orthogonal matrix can be generated by R16MAT_ORTH_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 16 ) X(N), the desired eigenvalues for the matrix.
!
!    Input, real ( kind = 16 ) Q(N,N), the eigenvector matrix of A.
!
!    Output, real ( kind = 16 ) A(N,N), a symmetric matrix with
!    eigenvalues X and eigenvectors the columns of Q.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 16 ) q(n,n)
  real ( kind = 16 ) x(n)
!
!  Set A = Q * Lambda * Q'.
!
  a(1:n,1:n) = 0.0_16

  do i = 1, n
    do j = 1, n
      do k = 1, n
        a(i,j) = a(i,j) + q(i,k) * x(k) * q(j,k)
      end do
    end do
  end do

  return
end
subroutine r16mat_symm_jacobi ( n, a )

!*****************************************************************************80
!
!! R16MAT_SYMM_JACOBI applies Jacobi eigenvalue iteration to a symmetric matrix.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    This code was modified so that it treats as zero the off-diagonal
!    elements that are sufficiently close to, but not exactly, zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input/output, real ( kind = 16 ) A(N,N), a symmetric N by N matrix.
!    On output, the matrix has been overwritten by an approximately
!    diagonal matrix, with the eigenvalues on the diagonal.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) c
  real ( kind = 16 ) r16mat_norm_fro
  real ( kind = 16 ), parameter :: eps = 0.00001_16
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 16 ) norm_fro
  real ( kind = 16 ) s
  real ( kind = 16 ) sum2
  real ( kind = 16 ) t
  real ( kind = 16 ) t1
  real ( kind = 16 ) t2
  real ( kind = 16 ) u

  norm_fro = r16mat_norm_fro ( n, n, a )

  it = 0

  do

    it = it + 1

    do i = 1, n
      do j = 1, i - 1

        if ( eps * norm_fro < abs ( a(i,j) ) + abs ( a(j,i) ) ) then

          u = ( a(j,j) - a(i,i) ) / ( a(i,j) + a(j,i) )

          t = sign ( 1.0_16, u ) / ( abs ( u ) + sqrt ( u * u + 1.0_16 ) )
          c = 1.0_16 / sqrt ( t * t + 1.0_16 )
          s = t * c
!
!  A -> A * Q.
!
          do k = 1, n
            t1 = a(i,k)
            t2 = a(j,k)
            a(i,k) = t1 * c - t2 * s
            a(j,k) = t1 * s + t2 * c
          end do
!
!  A -> QT * A
!
          do k = 1, n
            t1 = a(k,i)
            t2 = a(k,j)
            a(k,i) = c * t1 - s * t2
            a(k,j) = s * t1 + c * t2
          end do

        end if
      end do
    end do
!
!  Test the size of the off-diagonal elements.
!
    sum2 = 0.0_16
    do i = 1, n
      do j = 1, i - 1
        sum2 = sum2 + abs ( a(i,j) )
      end do
    end do

    if ( sum2 <= eps * ( norm_fro + 1.0_16 ) ) then
      exit
    end if

    if ( it_max <= it ) then
      exit
    end if

  end do

  return
end
subroutine r16mat_to_r16plu ( n, a, pivot, lu, info )

!*****************************************************************************80
!
!! R16MAT_TO_R16PLU factors a general R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    This routine is a simplified version of the LINPACK routine DGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 16 ) A(N,N), the matrix to be factored.
!
!    Output, integer ( kind = 4 ) PIVOT(N), a vector of pivot indices.
!
!    Output, real ( kind = 16 ) LU(N,N), an upper triangular matrix U and
!    the multipliers L which were used to obtain it.  The factorization
!    can be written A = L * U, where L is a product of permutation and
!    unit lower triangular matrices and U is upper triangular.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 16 ) lu(n,n)
  real ( kind = 16 ) temp

  lu(1:n,1:n) = a(1:n,1:n)

  info = 0

  do k = 1, n - 1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k + 1, n
      if ( abs ( lu(l,k) ) < abs ( lu(i,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( lu(l,k) == 0.0_16 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R16MAT_TO_R16PLU - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      temp    = lu(l,k)
      lu(l,k) = lu(k,k)
      lu(k,k) = temp
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    lu(k+1:n,k) = -lu(k+1:n,k) / lu(k,k)
!
!  Row elimination with column indexing.
!
    do j = k + 1, n

      if ( l /= k ) then
        temp    = lu(l,j)
        lu(l,j) = lu(k,j)
        lu(k,j) = temp
      end if

      lu(k+1:n,j) = lu(k+1:n,j) + lu(k+1:n,k) * lu(k,j)

    end do

  end do

  pivot(n) = n

  if ( lu(n,n) == 0.0_16 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R16MAT_TO_R16PLU - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
function r16mat_trace ( n, a )

!*****************************************************************************80
!
!! R16MAT_TRACE computes the trace of an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The trace of a square matrix is the sum of the diagonal elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, real ( kind = 16 ) A(N,N), the matrix whose trace is desired.
!
!    Output, real ( kind = 16 ) R16MAT_TRACE, the trace of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 16 ) r16mat_trace

  r16mat_trace = 0.0_16
  do i = 1, n
    r16mat_trace = r16mat_trace + a(i,i)
  end do

  return
end
subroutine r16mat_transpose_in_place ( n, a )

!*****************************************************************************80
!
!! R16MAT_TRANSPOSE_IN_PLACE transposes a square matrix in place.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns
!    of the matrix A.
!
!    Input/output, real ( kind = 16 ) A(N,N), the matrix to be transposed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 16 ) t

  do j = 1, n
    do i = 1, j - 1
      t      = a(i,j)
      a(i,j) = a(j,i)
      a(j,i) = t
    end do
  end do

  return
end
subroutine r16mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R16MAT_TRANSPOSE_PRINT prints an R16MAT, transposed.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 16 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  character ( len = * )  title

  call r16mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r16mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R16MAT_TRANSPOSE_PRINT_SOME prints some of an R16MAT, transposed.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 16 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r16mat_u_inverse ( n, a, b )

!*****************************************************************************80
!
!! R16MAT_U_INVERSE inverts an upper triangular R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    An upper triangular matrix is a matrix whose only nonzero entries
!    occur on or above the diagonal.
!
!    The inverse of an upper triangular matrix is an upper triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 16 ) A(N,N), the upper triangular matrix.
!
!    Output, real ( kind = 16 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = n, 1, -1

    do i = n, 1, -1

      if ( j < i ) then
        b(i,j) = 0.0_16
      else if ( i == j ) then
        b(i,j) = 1.0_16 / a(i,j)
      else
        b(i,j) = - dot_product ( a(i,i+1:j), b(i+1:j,j) ) / a(i,i)
      end if

    end do
  end do

  return
end
subroutine r16mat_u1_inverse ( n, a, b )

!*****************************************************************************80
!
!! R16MAT_U1_INVERSE inverts a unit upper triangular R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    A unit upper triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's below the main diagonal.
!
!    The inverse of a unit upper triangular matrix is also
!    a unit upper triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call r16mat_u1_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 16 ) A(N,N), the unit upper triangular matrix.
!
!    Output, real ( kind = 16 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  real ( kind = 16 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = n, 1, -1

    do i = n, 1, -1

      if ( j < i ) then
        b(i,j) = 0.0_16
      else if ( i == j ) then
        b(i,j) = 1.0_16
      else
        b(i,j) = - dot_product ( a(i,i+1:j), b(i+1:j,j) )
      end if

    end do
  end do

  return
end
subroutine r16mat_uniform ( m, n, a, b, seed, r )

!*****************************************************************************80
!
!! R16MAT_UNIFORM fills an R16MAT with scaled pseudorandom numbers.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer  ( kind = 4 )M, N, the number of rows and columns in
!    the array.
!
!    Input, real ( kind = 16 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 16 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a
  real ( kind = 16 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 16 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = a + ( b - a ) * real ( seed, kind = 16 ) * 4.656612875E-10_16

    end do
  end do

  return
end
subroutine r16mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R16MAT_UNIFORM_01 fills an R16MAT with unit pseudorandom numbers.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 16 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 16 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 16 ) * 4.656612875E-10_16

    end do
  end do

  return
end
subroutine r16mat_vand2 ( n, x, a )

!*****************************************************************************80
!
!! R16MAT_VAND2 returns the N by N row Vandermonde matrix A.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!    The row Vandermonde matrix returned by this routine reads "across"
!    rather than down.  In particular, each row begins with a 1, followed by
!    some value X, followed by successive powers of X.
!
!  Formula:
!
!    A(I,J) = X(I)**(J-1)
!
!  Properties:
!
!    A is nonsingular if, and only if, the X values are distinct.
!
!    The determinant of A is
!
!      det(A) = product ( 2 <= I <= N ) (
!        product ( 1 <= J <= I-1 ) ( ( X(I) - X(J) ) ) ).
!
!    The matrix A is generally ill-conditioned.
!
!  Example:
!
!    N = 5, X = (2, 3, 4, 5, 6)
!
!    1 2  4   8   16
!    1 3  9  27   81
!    1 4 16  64  256
!    1 5 25 125  625
!    1 6 36 216 1296
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix desired.
!
!    Input, real ( kind = 16 ) X(N), the values that define A.
!
!    Output, real ( kind = 16 ) A(N,N), the N by N row Vandermonde matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 16 ) x(n)

  do i = 1, n
    do j = 1, n

      if ( j == 1 .and. x(i) == 0.0_16 ) then
        a(i,j) = 1.0_16
      else
        a(i,j) = x(i)**(j-1)
      end if

    end do
  end do

  return
end
subroutine r16mat_zero ( m, n, a )

!*****************************************************************************80
!
!! R16MAT_ZERO zeroes an R16MAT.
!
!  Discussion:
!
!    An R16MAT is an array of R16 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Output, real ( kind = 16 ) A(M,N), the matrix of zeroes.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 16 ) a(m,n)

  a(1:m,1:n) = 0.0_16

  return
end
subroutine r16vec_heap_a ( n, a )

!*****************************************************************************80
!
!! R16VEC_HEAP_A reorders an R16VEC into an ascending heap.
!
!  Discussion:
!
!    An R16VEC is a vector of R16's.
!
!    An ascending heap is an array A with the property that, for every index J,
!    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, real ( kind = 16 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  real ( kind = 16 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the smaller of the two values,
!  and update M if necessary.
!
        if ( a(m+1) < a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the small descendant is smaller than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( key <= a(m) ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot.
!
    a(ifree) = key

  end do

  return
end
subroutine r16vec_heap_d ( n, a )

!*****************************************************************************80
!
!! R16VEC_HEAP_D reorders an R16VEC into an descending heap.
!
!  Discussion:
!
!    An R16VEC is a vector of R16's.
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, real ( kind = 16 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  real ( kind = 16 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine r16vec_print ( n, a, title )

!*****************************************************************************80
!
!! R16VEC_PRINT prints an R16VEC.
!
!  Discussion:
!
!    An R16VEC is a vector of R16's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 16 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g24.16)' ) i, ':', a(i)
  end do

  return
end
subroutine r16vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R16VEC_PRINT_PART prints "part" of an R16VEC.
!
!  Discussion:
!
!    An R16VEC is a vector of R16's.
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 16 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g24.16)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g24.16)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ........................'
    i = n
    write ( *, '(2x,i8,a,1x,g24.16)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g24.16)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g24.16,2x,a)' ) i, ':', a(i), '...more entries...'

  end if

  return
end
subroutine r16vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! R16VEC_PRINT_SOME prints "some" of an R16VEC.
!
!  Discussion:
!
!    An R16VEC is a vector of R16's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 16 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices
!    to print.  The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 16 ) a(n)
  integer( kind = 4 ) i
  integer( kind = 4 ) i_hi
  integer( kind = 4 ) i_lo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,a,1x,g24.16)' ) i, ':', a(i)
  end do

  return
end
subroutine r16vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R16VEC_UNIFORM_01 returns a unit pseudorandom R16VEC.
!
!  Discussion:
!
!    An R16VEC is a vector of R16's.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2010
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
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 16 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 16 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R16VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 16 ) * 4.656612875D-10

  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
