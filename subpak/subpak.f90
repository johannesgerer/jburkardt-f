subroutine angle_shift ( alpha, beta, gamma )

!*****************************************************************************80
!
!! ANGLE_SHIFT shifts angle ALPHA to lie between BETA and BETA+2PI.
!
!  Discussion:
!
!    The input angle ALPHA is shifted by multiples of 2 * PI to lie
!    between BETA and BETA+2*PI.
!
!    The resulting angle GAMMA has all the same trigonometric function
!    values as ALPHA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the angle to be shifted.
!
!    Input, real ( kind = 8 ) BETA, defines the lower endpoint of
!    the angle range.
!
!    Output, real ( kind = 8 ) GAMMA, the shifted angle.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  if ( alpha < beta ) then
    gamma = beta - mod ( beta - alpha, 2.0D+00 * pi ) + 2.0D+00 * pi
  else
    gamma = beta + mod ( alpha - beta, 2.0D+00 * pi )
  end if

  return
end
subroutine angle_shift_deg ( alpha, beta, gamma )

!*****************************************************************************80
!
!! ANGLE_SHIFT_DEG shifts angle ALPHA to lie between BETA and BETA+360.
!
!  Discussion:
!
!    The input angle ALPHA is shifted by multiples of 360 to lie
!    between BETA and BETA+360.
!
!    The resulting angle GAMMA has all the same trigonometric function
!    values as ALPHA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the angle to be shifted, in degrees.
!
!    Input, real ( kind = 8 ) BETA, defines the lower endpoint of
!    the angle range.
!
!    Output, real ( kind = 8 ) GAMMA, the shifted angle.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma

  if ( alpha < beta ) then
    gamma = beta - mod ( beta - alpha, 360.0D+00 ) + 360.0D+00
  else
    gamma = beta + mod ( alpha - beta, 360.0D+00 )
  end if

  return
end
subroutine angle_to_rgb ( angle, r, g, b )

!*****************************************************************************80
!
!! ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle in the color hexagon.  
!    The sextants are defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real ( kind = 8 ) R, G, B, RGB specifications for the color 
!    that lies at the given angle, on the perimeter of the color hexagon.  One
!    value will be 1, and one value will be 0.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle2
  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ), parameter :: degrees_to_radians = &
    3.141592653589793D+00 / 180.0D+00
  real ( kind = 8 ) r

  angle = mod ( angle, 360.0D+00 )

  if ( angle < 0.0D+00 ) then
    angle = angle + 360.0D+00
  end if

  if ( angle <= 60.0D00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = 1.0D+00
    g = tan ( angle2 )
    b = 0.0D+00

  else if ( angle <= 120.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = cos ( angle2 ) / sin ( angle2 )
    g = 1.0D+00
    b = 0.0D+00

  else if ( angle <= 180.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = 1.0D+00
    b = tan ( angle2 )

  else if ( angle <= 240.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = cos ( angle2 ) / sin ( angle2 )
    b = 1.0D+00

  else if ( angle <= 300.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = tan ( angle2 )
    g = 0.0D+00
    b = 1.0D+00

  else if ( angle <= 360.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = 1.0D+00
    g = 0.0D+00
    b = cos ( angle2 ) / sin ( angle2 )

  end if

  return
end
subroutine axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, nticks )

!*****************************************************************************80
!
!! AXIS_LIMITS returns "nice" axis limits for a plot.
!
!  Discussion:
!
!    The routine is given information about the range of a variable, and
!    the number of divisions desired.  It returns suggestions for
!    labeling a plotting axis for the variable, including the
!    starting and ending points, the length of a single division,
!    and a suggested tick marking for the axis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XMIN, XMAX, the lower and upper values that
!    must be included on the axis.  XMIN must be less than XMAX.
!
!    Input, integer ( kind = 4 ) NDIVS, the number of divisions desired along
!    the axis.
!
!    Output, real ( kind = 8 ) PXMIN, PXMAX, the recommended lower and upper
!    axis bounds.  It will be the case that PXMIN <= XMIN < XMAX <= PXMAX.
!
!    Output, real ( kind = 8 ) PXDIV, the recommended size of a single division.
!
!    Output, integer ( kind = 4 ) NTICKS, a suggested number of ticks to use,
!    if subdividing each of the NDIVS divisions of the axis.
!
  implicit none

  integer ( kind = 4 ), parameter :: nsteps = 5

  real ( kind = 8 ) best
  real ( kind = 8 ) good
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intlog
  integer ( kind = 4 ), dimension ( nsteps ) :: iticks = (/ 5, 4, 4, 5, 5 /)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ndivs
  integer ( kind = 4 ) nticks
  real ( kind = 8 ) pxmax
  real ( kind = 8 ) pxmax2
  real ( kind = 8 ) pxmin
  real ( kind = 8 ) pxmin2
  real ( kind = 8 ) pxdiv
  real ( kind = 8 ) pxdiv2
  real ( kind = 8 ) r8_log_10
  real ( kind = 8 ) reldif
  real ( kind = 8 ), dimension ( nsteps ) :: steps = (/ &
    1.0D+00,  2.0D+00,  4.0D+00,  5.0D+00, 10.0D+00 /)
  real ( kind = 8 ) temp
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  if ( xmin == xmax ) then
    xmin = xmin - 0.5D+00
    xmax = xmax + 0.5D+00
  else if ( xmax < xmin ) then
    temp = xmin
    xmin = xmax
    xmax = temp
  end if

  if ( ndivs <= 0 ) then
    ndivs = 5
  end if
!
!  Set RELDIF, the size of the X interval divided by the largest X.
!
  if ( xmax /= xmin ) then
    reldif = ( xmax - xmin ) / max ( abs ( xmax ), abs ( xmin ) )
  else
    reldif = 0.0D+00
  end if
!
!  If RELDIF tells us that XMIN and XMAX are extremely close,
!  do some simple things.
!
  if ( reldif < 0.00001D+00 ) then

    if ( xmax == 0.0D+00 ) then

      pxdiv = 1.0D+00

    else

      intlog = int ( r8_log_10 ( xmax ) )

      if ( intlog < 0 ) then
        intlog = intlog - 1
      end if

      pxdiv = 10.0D+00**intlog

      if ( 1.0D+00 < pxdiv ) then
        pxdiv = 1.0D+00
      end if

    end if

    nticks = 5
    pxmin = xmax - real ( ndivs / 2, kind = 8 ) * pxdiv
    pxmax = xmax + real ( ndivs - ( ndivs / 2 ), kind = 8 ) * pxdiv
!
!  But now handle the more general case, when XMIN and XMAX
!  are relatively far apart.
!
  else

    best = -999.0D+00
!
!  On second loop, increase INTLOG by 1.
!
    do j = 1, 2
!
!  Compute INTLOG, roughly the logarithm base 10 of the range
!  divided by the number of divisions.
!
      intlog = int ( r8_log_10 ( ( xmax - xmin ) &
        / real ( ndivs, kind = 8 ) ) ) + ( j - 1 )

      if ( xmax - xmin  < real ( ndivs, kind = 8 ) ) then
        intlog = intlog - 1
      end if
!
!  Now consider taking 1, 2, 4, 5 or 10 steps of size 10**INTLOG:
!
      do i = 1, nsteps
!
!  Compute the size of each step.
!
        pxdiv2 = steps(i) * 10.0D+00**intlog
!
!  Make sure NDIVS steps can reach from XMIN to XMAX, at least.
!
        if ( xmax <= xmin + ndivs * pxdiv2 ) then
!
!  Now decide where to start the axis.
!  Start the axis at PXMIN2, to the left of XMIN, and
!  representing a whole number of steps of size PXDIV2.
!
          if ( 0.0D+00 <= xmin ) then
            ival = int ( xmin / pxdiv2 )
          else
            ival = int ( xmin / pxdiv2 ) - 1
          end if

          pxmin2 = ival * pxdiv2
!
!  PXMAX2 is, of course, NDIVS steps above PXMIN2.
!
          pxmax2 = pxmin2 + ndivs * pxdiv2
!
!  Only consider going on if PXMAX2 is at least XMAX.
!
          if ( xmax <= pxmax2 ) then
!
!  Now judge this grid by the relative amount of wasted axis length.
!
            good = ( xmax - xmin ) / ( pxmax2 - pxmin2 )

            if ( best < good ) then
              best = good
              pxmax = pxmax2
              pxmin = pxmin2
              pxdiv = pxdiv2
              nticks = iticks(i)
            end if

          end if

        end if

      end do

    end do

  end if
!
!  If necessary, adjust the locations of PXMIN and PXMAX so that the
!  interval is more symmetric in containing XMIN through XMAX.
!
  do

    ilo = int ( xmin - pxmin ) / pxdiv
    ihi = int ( pxmax - xmax ) / pxdiv

    if ( ihi < ilo + 2 ) then
      exit
    end if

    pxmin = pxmin - pxdiv
    pxmax = pxmax - pxdiv

  end do

  return
end
subroutine bar_check ( digit, check )

!*****************************************************************************80
!
!! BAR_CHECK computes the check digit for a barcode.
!
!  Formula:
!
!    CHECK = SUM ( I = 1, 11, by 2's ) DIGIT(I)
!       + 3 * SUM ( I = 2, 10, by 2's ) DIGIT(I)
!
!    CHECK = MOD ( 10 - MOD ( CHECK, 10 ), 10 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT(12), entries 1 through 11 of DIGIT 
!    contain the digits of the bar code.  Each entry must be between 0 and 9.
!    The 12th digit should be the check digit.
!
!    Output, integer ( kind = 4 ) CHECK, the correct check digit.  If the bar
!    code is correct, then DIGIT(12) should equal CHECK.
!
  implicit none

  integer ( kind = 4 ) check
  integer ( kind = 4 ) digit(12)

  check = sum ( digit(1:11:2) ) + 3 * sum ( digit(2:10:2) )

  check = mod ( 10 - mod ( check, 10 ), 10 )

  return
end
subroutine bar_code ( digit, bar )

!*****************************************************************************80
!
!! BAR_CODE constructs the 113 character barcode from 11 digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) DIGIT(12).
!    On input, the first 11 entries of DIGIT contain a code to be
!    turned into a barcode.
!    On output, the 12-th entry of DIGIT is a check digit.
!
!    Output, character ( len = 113 ) BAR, the bar code corresponding to the
!    digit information.
!
  implicit none

  character ( len = 113 ) bar
  integer ( kind = 4 ) check
  character ( len = 7 ) codel
  character ( len = 7 ) coder
  integer ( kind = 4 ) digit(12)
  integer ( kind = 4 ) i
!
!  9 character quiet zone.
!
  bar(1:9) = '000000000'
!
!  3 character guard pattern.
!
  bar(10:12) = '101'
!
!  7 character product category.
!
  call bar_digit_code_left ( digit(1), codel )
  bar(13:19) = codel
!
!  35 characters contain the 5 digit manufacturer code.
!
  do i = 1, 5
    call bar_digit_code_left ( digit(i+1), codel )
    bar(20+(i-1)*7:20+(i-1)*7+6) = codel
  end do
!
!  Center guard pattern.
!
  bar(55:59) = '01010'
!
!  35 characters contain the 5 digit product code.
!
  do i = 1, 5
    call bar_digit_code_right ( digit(i+6), coder )
    bar(60+(i-1)*7:60+(i-1)*7+6) = coder
  end do
!
!  Compute the check digit.
!
  call bar_check ( digit, check )
  digit(12) = check

  call bar_digit_code_right ( digit(12), coder )
  bar(95:101) = coder
!
!  Guard pattern.
!
  bar(102:104) = '101'
!
!  Quiet zone.
!
  bar(105:113) = '000000000'

  return
end
subroutine bar_digit_code_left ( digit, codel )

!*****************************************************************************80
!
!! BAR_DIGIT_CODE_LEFT returns the 7 character left bar code for a digit.
!
!  Example:
!
!    DIGIT = 3
!    CODEL = '0111101'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT, the digit, between 0 and 9.
!
!    Output, character ( len = 7 ) CODEL, the left code for the digit.
!
  implicit none

  character ( len = 7 ) codel
  integer ( kind = 4 ) digit

  if ( digit == 0 ) then
    codel = '0001101'
  else if ( digit == 1 ) then
    codel = '0011001'
  else if ( digit == 2 ) then
    codel = '0010011'
  else if ( digit == 3 ) then
    codel = '0111101'
  else if ( digit == 4 ) then
    codel = '0100011'
  else if ( digit == 5 ) then
    codel = '0110001'
  else if ( digit == 6 ) then
    codel = '0101111'
  else if ( digit == 7 ) then
    codel = '0111011'
  else if ( digit == 8 ) then
    codel = '0110111'
  else if ( digit == 9 ) then
    codel = '0001011'
  else
    codel = '???????'
  end if

  return
end
subroutine bar_digit_code_right ( digit, coder )

!*****************************************************************************80
!
!! BAR_DIGIT_CODE_RIGHT returns the 7 character right bar code for a digit.
!
!  Example:
!
!    DIGIT = 3
!    CODER = '1000010'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT, the digit, between 0 and 9.
!
!    Output, character ( len = 7 ) CODER, the right code for the digit.
!
  implicit none

  character ( len = 7 ) coder
  integer ( kind = 4 ) digit

  if ( digit == 0 ) then
    coder = '1110010'
  else if ( digit == 1 ) then
    coder = '1100110'
  else if ( digit == 2 ) then
    coder = '1101100'
  else if ( digit == 3 ) then
    coder = '1000010'
  else if ( digit == 4 ) then
    coder = '1011100'
  else if ( digit == 5 ) then
    coder = '1001110'
  else if ( digit == 6 ) then
    coder = '1010000'
  else if ( digit == 7 ) then
    coder = '1000100'
  else if ( digit == 8 ) then
    coder = '1001000'
  else if ( digit == 9 ) then
    coder = '1110100'
  else
    coder = '???????'
  end if

  return
end
function bmi_english ( w_lb, h_ft, h_in )

!*****************************************************************************80
!
!! BMI_ENGLISH computes the body mass index given English measurements.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W_LB, the body weight in pounds.
!
!    Input, real ( kind = 8 ) H_FT, H_IN, the body height in feet and inches
!
!    Output, real ( kind = 8 ) BMI_ENGLISH, the body mass index.
!
  implicit none

  real ( kind = 8 ) bmi_english
  real ( kind = 8 ) bmi_metric
  real ( kind = 8 ) feet_to_meters
  real ( kind = 8 ) h_ft
  real ( kind = 8 ) h_in
  real ( kind = 8 ) h_m
  real ( kind = 8 ) pounds_to_kilograms
  real ( kind = 8 ) w_kg
  real ( kind = 8 ) w_lb

  w_kg = pounds_to_kilograms ( w_lb )

  h_m = feet_to_meters ( h_ft + ( h_in / 12.0D+00 ) )

  bmi_english = bmi_metric ( w_kg, h_m )

  return
end
function bmi_metric ( w_kg, h_m )

!*****************************************************************************80
!
!! BMI_METRIC computes the body mass index given metric measurements.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W_KG, the body weight in kilograms.
!
!    Input, real ( kind = 8 ) H_M, the body height in meters.
!
!    Output, real ( kind = 8 ) BMI_METRIC, the body mass index.
!
  implicit none

  real ( kind = 8 ) bmi_metric
  real ( kind = 8 ) h_m
  real ( kind = 8 ) w_kg

  bmi_metric = ( w_kg / h_m ) / h_m

  return
end
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT is TRUE if C is a decimal digit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, is TRUE if C is a digit.
!
  implicit none

  character c
  logical ch_is_digit

  if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
  end if

  return
end
function degrees_to_radians ( degrees )

!*****************************************************************************80
!
!! DEGREES_TO_RADIANS converts an angle measure from degrees to radians.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DEGREES, the angle measure in degrees.
!
!    Output, real ( kind = 8 ) DEGREES_TO_RADIANS, the angle measure in radians.
!
  implicit none

  real ( kind = 8 ) degrees
  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  degrees_to_radians = ( degrees / 180.0D+00 ) * pi

  return
end
function e_constant ( )

!*****************************************************************************80
!
!! E_CONSTANT returns the value of E.
!
!  Discussion:
!
!    "E" was named in honor of Euler, but is known as Napier's constant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) E_CONSTANT, the base of the natural 
!    logarithm system.
!
  implicit none

  real ( kind = 8 ) e_constant

  e_constant = 2.718281828459045D+00

  return
end
function euler_constant ( )

!*****************************************************************************80
!
!! EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
!
!  Discussion:
!
!    The Euler-Mascheroni constant is often denoted by a lower-case
!    Gamma.  Gamma is defined as
!
!      Gamma = limit ( M -> oo ) ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EULER_CONSTANT, the value of the 
!    Euler-Mascheroni constant.
!
  implicit none

  real ( kind = 8 ) euler_constant

  euler_constant = 0.5772156649015328D+00

  return
end
subroutine fac_div ( prime_num, npower1, npower2, npower3 )

!*****************************************************************************80
!
!! FAC_DIV divides two quantities represented as prime factors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER1(PRIME_NUM), the powers of primes
!    in the representation of the first quantity.
!
!    Input, integer ( kind = 4 ) NPOWER2(PRIME_NUM), the powers of primes
!    in the representation of the second quantity.
!
!    Output, integer ( kind = 4 ) NPOWER3(PRIME_NUM), the powers of primes
!    in the representation of the quotient.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) npower1(prime_num)
  integer ( kind = 4 ) npower2(prime_num)
  integer ( kind = 4 ) npower3(prime_num)

  npower3(1:prime_num) = npower1(1:prime_num) - npower2(1:prime_num)

  return
end
subroutine fac_gcd ( prime_num, npower1, npower2, npower3 )

!*****************************************************************************80
!
!! FAC_GCD finds the GCD of two products of prime factors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER1(PRIME_NUM), the powers of primes
!    in the representation of the first quantity.  All the powers
!    must be nonnegative.
!
!    Input, integer ( kind = 4 ) NPOWER2(PRIME_NUM), the powers of primes
!    in the representation of the second quantity.  All the powers
!    must be nonnegative.
!
!    Output, integer ( kind = 4 ) NPOWER3(PRIME_NUM), the powers of primes
!    in the representation of the GCD.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) npower1(prime_num)
  integer ( kind = 4 ) npower2(prime_num)
  integer ( kind = 4 ) npower3(prime_num)

  do i = 1, prime_num

    if ( npower1(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAC_GCD - Fatal error!'
      write ( *, '(a)' ) '  One of the powers is negative!'
      stop
    end if

    if ( npower2(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAC_GCD - Fatal error!'
      write ( *, '(a)' ) '  One of the powers is negative!'
      stop
    end if

    npower3(i) = min ( npower1(i), npower2(i) )

  end do

  return
end
subroutine fac_lcm ( prime_num, npower1, npower2, npower3 )

!*****************************************************************************80
!
!! FAC_LCM finds the LCM of two products of prime factors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER1(PRIME_NUM), the powers of primes
!    in the representation of the first quantity.
!
!    Input, integer ( kind = 4 ) NPOWER2(PRIME_NUM), the powers of primes
!    in the representation of the second quantity.
!
!    Output, integer ( kind = 4 ) NPOWER3(PRIME_NUM), the powers of primes
!    in the representation of the LCM.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) npower1(prime_num)
  integer ( kind = 4 ) npower2(prime_num)
  integer ( kind = 4 ) npower3(prime_num)

  do i = 1, prime_num

    if ( npower1(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAC_LCM - Fatal error!'
      write ( *, '(a)' ) '  One of the powers is negative!'
      stop
    end if

    if ( npower2(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAC_LCM - Fatal error!'
      write ( *, '(a)' ) '  One of the powers is negative!'
      stop
    end if

    npower3(i) = max ( npower1(i), npower2(i) )

  end do

  return
end
subroutine fac_mul ( prime_num, npower1, npower2, npower3 )

!*****************************************************************************80
!
!! FAC_MUL multiplies two quantities represented as prime factors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER1(PRIME_NUM), the powers of primes
!    in the representation of the first quantity.
!
!    Input, integer ( kind = 4 ) NPOWER2(PRIME_NUM), the powers of primes
!    in the representation of the second quantity.
!
!    Output, integer ( kind = 4 ) NPOWER3(PRIME_NUM), the powers of primes
!    in the representation of the product.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) npower1(prime_num)
  integer ( kind = 4 ) npower2(prime_num)
  integer ( kind = 4 ) npower3(prime_num)

  npower3(1:prime_num) = npower1(1:prime_num) + npower2(1:prime_num)

  return
end
subroutine fac_print ( prime_num, npower )

!*****************************************************************************80
!
!! FAC_PRINT prints a product of prime factors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER(PRIME_NUM), the powers of primes
!    in the representation of the quantity.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) npower(prime_num)
  integer ( kind = 4 ) prime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Prime     Power'
  write ( *, '(a)' ) ' '
  do i = 1, prime_num
    if ( npower(i) /= 0 ) then
      write ( *, '(i8,2x,i8)' ) prime(i), npower(i)
    end if
  end do

  return
end
subroutine fac_to_i4 ( prime_num, npower, intval )

!*****************************************************************************80
!
!! FAC_TO_I4 converts a product of prime factors into an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER(PRIME_NUM), the powers of primes
!    in the representation of the quantity.  If any of these powers
!    are negative, then INTVAL will be set to 0.
!
!    Output, integer ( kind = 4 ) INTVAL, the integer represented by the 
!    product of the prime factors.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) npower(prime_num)
  integer ( kind = 4 ) prime

  intval = 1
  do i = 1, prime_num

    if ( npower(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAC_TO_I4 - Fatal error!'
      write ( *, '(a)' ) '  One of the powers is negative!'
      stop
    end if

    intval = intval * prime(i)**npower(i)

  end do

  return
end
subroutine fac_to_rat ( prime_num, npower, top, bot )

!*****************************************************************************80
!
!! FAC_TO_RAT converts a prime factorization into a rational value.
!
!  Example:
!
!    Start with the prime factorization representation:
!
!      40/9 = 2^3 * 3^(-2) * 5
!
!    Input:
!
!      NPOWER = ( 3, -2, 1 )
!
!    Output:
!
!      TOP = 40 ( = 2^3 * 5^1 = PRIME(1)^3                 * PRIME(3)^1 )
!      BOT = 9  ( = 3^2       =             PRIME(2)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER(PRIME_NUM).  NPOWER(I) is the power of
!    the I-th prime in the prime factorization.  NPOWER(I) may
!    be positive or negative.
!
!    Output, integer ( kind = 4 ) TOP, BOT, the top and bottom of a rational 
!    value.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) bot
  integer ( kind = 4 ) i
  integer ( kind = 4 ) npower(prime_num)
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) top

  top = 1
  bot = 1
  do i = 1, prime_num
    if ( 0 < npower(i) ) then
      top = top * prime(i)**npower(i)
    else if ( npower(i) < 0 ) then
      bot = bot * prime(i)**(-npower(i))
    end if
  end do

  return
end
function feet_to_meters ( ft )

!*****************************************************************************80
!
!! FEET_TO_METERS converts a measurement in feet to meters.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) FT, the length in feet.
!
!    Output, real ( kind = 8 ) FEET_TO_METERS, the corresponding 
!    length in meters.
!
  implicit none

  real ( kind = 8 ) feet_to_meters
  real ( kind = 8 ) ft

  feet_to_meters = 0.0254D+00 * 12.0D+00 * ft

  return
end
function gauss_sum ( dim_num, n, amplitude, center, width, x )

!*****************************************************************************80
!
!! GAUSS_SUM evaluates a function that is the sum of Gaussians.
!
!  Discussion:
!
!    Gauss_Sum(X) = Sum ( 1 <= J <= Ngauss ) Amplitude(I) * exp ( -Arg )
!
!    where
!
!      Arg = sum ( 1 <= I <= DIM_NUM ) 
!        ( ( ( X(I) - Center(I,J) ) / Width(J) )^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of component Gaussian functions.
!
!    Input, real ( kind = 8 ) AMPLITUDE(N), CENTER(DIM_NUM,N), WIDTH(N),
!    the amplitude, center and width for the component Gaussian functions.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the point at which the function 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) GAUSS_SUM, the value of the function.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real ( kind = 8 ) amplitude(n)
  real ( kind = 8 ) arg
  real ( kind = 8 ) center(dim_num,n)
  real ( kind = 8 ) gauss_sum
  integer ( kind = 4 ) j
  real ( kind = 8 ) width(n)
  real ( kind = 8 ) x(dim_num)

  gauss_sum = 0.0D+00

  do j = 1, n

    arg = sum ( ( ( x(1:dim_num) - center(1:dim_num,j) ) / width(j) )**2 )

    gauss_sum = gauss_sum + amplitude(j) * exp ( -arg )

  end do

  return
end
subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 ) today
  integer ( kind = 4 ) values(8)
  character ( len = 5 ) zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) /  11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) /  30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) /  23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp                                    /   6.0D+00

  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( i4_huge, kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == i4_huge ) then
    seed = seed - 1
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
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2008
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
subroutine grid1 ( dim_num, nstep, x1, x2, x )

!*****************************************************************************80
!
!! GRID1 finds grid points between X1 and X2 in N dimensions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1 and X2.
!
!    Input, integer ( kind = 4 ) NSTEP, the number of points to be generated.
!    NSTEP must be at least 2.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), the first and last
!    points, between which the equally spaced points are
!    to be computed.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,NSTEP), the set of equally spaced
!    points.  Each column of X represents one point, with X(*,1) = X1
!    and X(*,NSTEP) = X2.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(dim_num,nstep)
  real ( kind = 8 ) x1(dim_num)
  real ( kind = 8 ) x2(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID1 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID1 - Fatal error!'
    write ( *, '(a)' ) '  NSTEP < 2.'
    write ( *, '(a,i8)' ) '  NSTEP = ', nstep
    stop
  end if

  do i = 1, nstep
    x(1:dim_num,i) = &
      ( real ( nstep - i,     kind = 8 ) * x1(1:dim_num)   &
      + real (         i - 1, kind = 8 ) * x2(1:dim_num) ) &
      / real ( nstep     - 1, kind = 8 )
  end do

  return
end
subroutine grid1n ( j, dim_num, nstep, x1, x2, x )

!*****************************************************************************80
!
!! GRID1N finds the I-th grid point between X1 and X2 in N dimensions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J, the number of the desired point.
!    Normally J would be between 1 and NSTEP, but that is
!    not necessary.  Note that J = 1 returns X1 and J = NSTEP
!    returns X2.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X, 
!    X1 and X2.
!
!    Input, integer ( kind = 4 ) NSTEP, this is the number of equally
!    spaced points that are between X1 and X2.  NSTEP must
!    be at least 2, because X1 and X2 are always included
!    in the set of points.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), the first and last
!    points, between which the equally spaced points lie.
!
!    Output, real ( kind = 8 ) X(DIM_NUM), the J-th grid point between X1
!    and X2.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep

  integer ( kind = 4 ) j
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) x1(dim_num)
  real ( kind = 8 ) x2(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID1N - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID1N - Fatal error!'
    write ( *, '(a)' ) '  NSTEP <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP = ', nstep
    stop
  end if

  x(1:dim_num) = ( real ( nstep - j,     kind = 8 ) * x1(1:dim_num) &
                 + real (         j - 1, kind = 8 ) * x2(1:dim_num) ) &
                 / real ( nstep     - 1, kind = 8 )

  return
end
subroutine grid2 ( j1, j2, dim_num, nstep, x1, x2, x )

!*****************************************************************************80
!
!! GRID2 computes grid points between X1 and X2 in N dimensions.
!
!  Discussion:
!
!    GRID2 computes grid points between X1 and X2 in N dimensions.
!
!    However, X1 need not be the first point computed, nor X2 the last.
!    The user must specify the steps on which X1 and X2 are passed
!    through.  These steps may even be outside the range of 1 through NSTEP.
!
!    We assume that a set of equally spaced points have
!    been drawn on the line through X1 and X2, and that
!    they have been numbered, with X1 labeled J1 and X2
!    labeled J2.  J1 or J2 may be between 1 and NSTEP,
!    in which case X1 or X2 will actually be returned in the
!    X array, but there is no requirement that J1 or J2
!    satisfy this condition.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J1, J2.  J1 specifies the step on which
!    X1 would be computed, and similarly for J2.  
!    J1 and J2 must be distinct.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1 and X2.
!
!    Input, integer ( kind = 4 ) NSTEP, this is the number of equally
!    spaced points that are to be generated.
!    NSTEP should be at least 1.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), the points that define
!    the line along which the equally spaced points are generated, and
!    which may or may not be included in the set of computed points.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,NSTEP), the set of equally spaced
!    points.  Each column of X represents one point.
!    If 1 <= J1 <= NSTEP, then X(*,J1) = X1, and similarly for J2.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) x(dim_num,nstep)
  real ( kind = 8 ) x1(dim_num)
  real ( kind = 8 ) x2(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID2 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( j1 == j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID2 - Fatal error!'
    write ( *, '(a)' ) '  J1 = J2, leading to zero denominator.'
    write ( *, '(a,i8)' ) '  J1 = ', j1
    write ( *, '(a,i8)' ) '  J2 = ', j2
    stop
  end if

  do j = 1, nstep
    do i = 1, dim_num
      x(i,j) = ( real ( j2 - j,      kind = 8 ) * x1(i)   &
               + real (      j - j1, kind = 8 ) * x2(i) ) &
               / real ( j2     - j1, kind = 8 )
    end do
  end do

  return
end
subroutine grid2n ( j, j1, j2, dim_num, x1, x2, x )

!*****************************************************************************80
!
!! GRID2N computes one grid point between X1 and X2 in N dimensions.
!
!  Discussion:
!
!    However, X1 need not be the first point computed, nor X2 the last.
!    The user must specify the steps on which X1 and X2 are passed through.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J, the J coordinate of the desired point.
!    Note that if J = J1, X will be returned as X1, and if
!    J = J2, X will be returned as X2.
!
!    Input, integer ( kind = 4 ) J1, J2.  J1 specifies the step on which
!    X1 would be computed, and similarly for J2.  That is,
!    we assume that a set of equally spaced points have
!    been drawn on the line through X1 and X2, and that
!    they have been numbered, with X1 labeled J1 and X2
!    labeled J2.  J1 and J2 must be distinct.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1 and X2.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), the points that define
!    the line along which the equally spaced points are
!    generated, and which may or may not be included in the
!    set of computed points.
!
!    Output, real ( kind = 8 ) X(DIM_NUM).  X(I) is the J-th point from the
!    set of equally spaced points.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) x1(dim_num)
  real ( kind = 8 ) x2(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID2N - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( j1 == j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID2N - Fatal error!'
    write ( *, '(a)' ) '  J1 = J2, leading to zero denominator.'
    write ( *, '(a,i8)' ) '  J1 = ', j1
    write ( *, '(a,i8)' ) '  J2 = ', j2
    stop
  end if

  do i = 1, dim_num
    x(i) = ( real ( j2 - j,      kind = 8 ) * x1(i)   &
           + real (      j - j1, kind = 8 ) * x2(i) ) &
           / real ( j2     - j1, kind = 8 )
  end do

  return
end
subroutine grid3 ( dim_num, nstep1, nstep2, x1, x2, x3, x )

!*****************************************************************************80
!
!! GRID3 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
!
!  Discussion:
!
!    The line between X1 and X2 will have NSTEP1 points generated along 
!    it, and the line between X1 and X3 will have NSTEP2 points generated
!    along it.
!
!    Fixing the second and third indices of X represents one point, with
!    the following special values:
!
!      X(*,1,1)      = X1
!      X(*,NSTEP1,1) = X2
!      X(*,1,NSTEP2) = X3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1, 
!    X2 and X3.
!
!    Input, integer ( kind = 4 ) NSTEP1, NSTEP2.  These are the number of
!    equally spaced points to generate in the first and second
!    directions.  NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
!    X3 are always included in the set of points.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,NSTEP1,NSTEP2), the set of equally
!    spaced points.  
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep1
  integer ( kind = 4 ) nstep2

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) psi1
  real ( kind = 8 ) psi2
  real ( kind = 8 ) psi3
  real ( kind = 8 ) x(dim_num,nstep1,nstep2)
  real ( kind = 8 ) x1(dim_num)
  real ( kind = 8 ) x2(dim_num)
  real ( kind = 8 ) x3(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep1 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3 - Fatal error!'
    write ( *, '(a)' ) '  NSTEP1 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3 - Fatal error!'
    write ( *, '(a)' ) '  NSTEP2 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
    stop
  end if

  do j = 1, nstep1

    psi2 = real ( j - 1, kind = 8 ) &
         / real ( nstep1 - 1, kind = 8 )

    do k = 1, nstep2

      psi3 = real ( k - 1, kind = 8 ) &
           / real ( nstep2 - 1, kind = 8 )

      psi1 = 1.0D+00 - psi2 - psi3

      x(1:dim_num,j,k) = psi1 * x1(1:dim_num) &
                       + psi2 * x2(1:dim_num) &
                       + psi3 * x3(1:dim_num)

    end do
  end do

  return
end
subroutine grid3n ( j, k, dim_num, nstep1, nstep2, x1, x2, x3, x )

!*****************************************************************************80
!
!! GRID3N computes a parallelogram grid on 3 points in N dimensions.
!
!  Discussion:
!
!    The line between X1 and X2 will have NSTEP1
!    points generated along it, and the line between X1 and
!    X3 will have NSTEP2 points generated along it.
!
!    The following special values are:
!
!      J       K         X
!
!      1       1         X1
!      NSTEP1  1         X2
!      1       NSTEP2    X3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J, K, the parallelogram coordinates
!    of the point.  J measures steps from X1 to X2, and
!    K measures steps from X1 to X3.  Normally, J would
!    be between 1 and NSTEP1, K between 1 and NSTEP2,
!    but this is not necessary.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1, 
!    X2 and X3.
!
!    Input, integer ( kind = 4 ) NSTEP1, NSTEP2.  These are the number of
!    equally spaced points to generate in the first and second
!    directions.  NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
!    X3 are always included in the set of points.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
!    Output, real ( kind = 8 ) X(DIM_NUM), the point with coordinates (J,K)
!    from the the set of equally  spaced points.  
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep1
  integer ( kind = 4 ) nstep2

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) psi1
  real ( kind = 8 ) psi2
  real ( kind = 8 ) psi3
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) x1(dim_num)
  real ( kind = 8 ) x2(dim_num)
  real ( kind = 8 ) x3(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3N - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep1 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3N - Fatal error!'
    write ( *, '(a)' ) '  NSTEP1 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3N - Fatal error!'
    write ( *, '(a)' ) '  NSTEP2 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
    stop
  end if

  psi2 = real ( j - 1, kind = 8 ) &
       / real ( nstep1 - 1, kind = 8 )

  psi3 = real ( k - 1, kind = 8 ) &
       / real ( nstep2 - 1, kind = 8 )

  psi1 = 1.0D+00 - psi2 - psi3

  x(1:dim_num) = psi1 * x1(1:dim_num) &
               + psi2 * x2(1:dim_num) &
               + psi3 * x3(1:dim_num)

  return
end
subroutine grid4 ( j1, j2, k1, k2, dim_num, nstep1, nstep2, x1, x2, x3, x )

!*****************************************************************************80
!
!! GRID4 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
!
!  Discussion:
!
!    Unlike GRID3, GRID4 does not necessarily place X1 at the
!    "origin" of the parallelogram, with X2 and X3 set at the
!    extreme J and K coordinates.  Instead, the user is free
!    to specify the J and K coordinates of the points, although
!    they are required to lie on a subparallelogram of the
!    larger one.
!
!    The line through X1 and X2 will have NSTEP1
!    points generated along it, and the line through X1 and
!    X3 will have NSTEP2 points generated along it.
!
!    If we imagine that the
!    main parallelogram is drawn first, with coordinate
!    ranges 1 <= J <= NSTEP1 and 1 <= K <= NSTEP2, then
!    these indices determine the (J,K) coordinates of the
!    three points, namely:
!
!      X1 : (J1,K1)
!      X2 : (J2,K1)
!      X3 : (J1,K2)
!
!    Of course, we actually start with the points X1, X2,
!    and X3, and they define a parallelogram and a (J,K)
!    coordinate system over the plane containing them.  We
!    then are free to consider the parallelogram defined
!    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
!    which may or may not contain any of the points X1, X2
!    and X3.
!
!    Assuming that the indices J1, J2, K1 and K2 are "within
!    bounds", the following special values will be computed:
!
!      X(*,J1,K1) = X1
!      X(*,J2,K1) = X2
!      X(*,J1,K2) = X3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J1, J2, K1, K2, the indices.  
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1, 
!    X2 and X3.
!
!    Input, integer ( kind = 4 ) NSTEP1, NSTEP2.  These are the number of
!    equally spaced points to generate in the first and second
!    directions.  NSTEP1 and NSTEP2 should be at least 1.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,NSTEP1,NSTEP2), the set of equally
!    spaced points.  Fixing the second and third indices
!    of X represents one point.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep1
  integer ( kind = 4 ) nstep2

  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) psi1
  real ( kind = 8 ) psi2
  real ( kind = 8 ) psi3
  real ( kind = 8 ) x(dim_num,nstep1,nstep2)
  real ( kind = 8 ) x1(dim_num)
  real ( kind = 8 ) x2(dim_num)
  real ( kind = 8 ) x3(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep1 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4 - Fatal error!'
    write ( *, '(a)' ) '  NSTEP1 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4 - Fatal error!'
    write ( *, '(a)' ) '  NSTEP2 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
    stop
  end if

  if ( j1 == j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4 - Fatal error!'
    write ( *, '(a)' ) '  J1 = J2'
    write ( *, '(a,i8)' ) '  J1 = ', j1
    write ( *, '(a,i8)' ) '  J2 = ', j2
    stop
  end if

  if ( k1 == k2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4 - Fatal error!'
    write ( *, '(a)' ) '  K1 = K2'
    write ( *, '(a,i8)' ) '  K1 = ', k1
    write ( *, '(a,i8)' ) '  K2 = ', k2
    stop
  end if

  do j = 1, nstep1

    psi2 = real (  j - j1, kind = 8 ) &
         / real ( j2 - j1, kind = 8 )

    do k = 1, nstep2

      psi3 = real (  k - k1, kind = 8 ) &
           / real ( k2 - k1, kind = 8 )

      psi1 = 1.0D+00 - psi2 - psi3

      x(1:dim_num,j,k) = psi1 * x1(1:dim_num) &
                       + psi2 * x2(1:dim_num) &
                       + psi3 * x3(1:dim_num)

    end do
  end do

  return
end
subroutine grid4n ( j, j1, j2, k, k1, k2, dim_num, nstep1, nstep2, x1, x2, &
  x3, x )

!*****************************************************************************80
!
!! GRID4N computes a single point on a parallelogram grid in N space.
!
!  Discussion:
!
!    The computation is identical to that of GRID4, except that
!    only one point at a time is computed.
!
!    The line through X1 and X2 will have NSTEP1
!    points generated along it, and the line through X1 and
!    X3 will have NSTEP2 points generated along it.
!
!    The following special values will be computed:
!
!      J  K  X
!
!      J1 K1 X1
!      J2 K2 X2
!      J1 K2 X3
!
!    If we imagine that the main parallelogram is drawn first, with 
!    coordinate ranges 1 <= J <= NSTEP1 and 1 <= K <= NSTEP2, then
!    the indices J and K determine the (J,K) coordinates of the
!    three points X1, X2, and X3, namely:
!
!      X1 : (J1,K1)
!      X2 : (J2,K1)
!      X3 : (J1,K2)
!
!    Of course, we actually start with the points X1, X2,
!    and X3, and they define a parallelogram and an (J,K)
!    coordinate system over the plane containing them.  We
!    then are free to consider the parallelogram defined
!    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
!    which may or may not contain any of the points X1, X2
!    and X3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J, the J coordinate of the point X.
!
!    Input, integer ( kind = 4 ) J1, J2.  See discussion.
!
!    Input, integer ( kind = 4 ) K, the K coordinate of the point X.
!
!    Input, integer ( kind = 4 ) K1, K2.  See discussion.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points 
!    X, X1, X2 and X3.
!
!    Input, integer ( kind = 4 ) NSTEP1, NSTEP2.  These are the number of
!    equally spaced points generated in the first and second
!    directions.
!    NSTEP1 and NSTEP2 should be at least 1.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
!    Output, real ( kind = 8 ) X(DIM_NUM), the point whose parallelogram
!    coordinates are (J,K).
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep1
  integer ( kind = 4 ) nstep2

  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) psi1
  real ( kind = 8 ) psi2
  real ( kind = 8 ) psi3
  real ( kind = 8 ) x(dim_num)
  real ( kind = 8 ) x1(dim_num)
  real ( kind = 8 ) x2(dim_num)
  real ( kind = 8 ) x3(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4N - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep1 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4N - Fatal error!'
    write ( *, '(a)' ) '  NSTEP1 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4N - Fatal error!'
    write ( *, '(a)' ) '  NSTEP2 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
    stop
  end if

  if ( j1 == j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4N - Fatal error!'
    write ( *, '(a)' ) '  J1 = J2'
    write ( *, '(a,i8)' ) '  J1 = ', j1
    write ( *, '(a,i8)' ) '  J2 = ', j2
    stop
  end if

  if ( k1 == k2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4N - Fatal error!'
    write ( *, '(a)' ) '  K1 = K2'
    write ( *, '(a,i8)' ) '  K1 = ', k1
    write ( *, '(a,i8)' ) '  K2 = ', k2
    stop
  end if

  psi2 = real ( j  - j1, kind = 8 ) &
       / real ( j2 - j1, kind = 8 )

  psi3 = real ( k  - k1, kind = 8 ) &
       / real ( k2 - k1, kind = 8 )

  psi1 = 1.0D+00 - psi2 - psi3

  x(1:dim_num) = psi1 * x1(1:dim_num) &
               + psi2 * x2(1:dim_num) &
               + psi3 * x3(1:dim_num)

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_sign ( x )

!*****************************************************************************80
!
!! I4_SIGN evaluates the sign of an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
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
!    Input, integer ( kind = 4 ) X, the number whose sign is desired.
!
!    Output, integer ( kind = 4 ) I4_SIGN, the sign of X:
!
  implicit none

  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    i4_sign = - 1
  else
    i4_sign = + 1
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine i4_to_digits_decimal ( i, n, digit )

!*****************************************************************************80
!
!! I4_TO_DIGITS_DECIMAL determines the last N decimal digits of an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the integer to be analyzed.
!
!    Input, integer ( kind = 4 ) N, the number of digits to determine.
!
!    Output, integer ( kind = 4 ) DIGIT(N), the last N decimal digits of I.
!    DIGIT(I) is the "coefficient" of 10**(I-1).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) digit(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_copy
  integer ( kind = 4 ) j

  i_copy = abs ( i )

  do j = 1, n
    digit(j) = mod ( i_copy, 10 )
    i_copy = ( i_copy - digit(j) ) / 10
  end do

  return
end
subroutine i4_to_fac ( intval, prime_num, npower )

!*****************************************************************************80
!
!! I4_TO_FAC converts an I4 into a product of prime factors.
!
!  Discussion:
!
!    This routine will fail if the input integer is not positive,
!    or if PRIME_NUM is too small to account for the factors of the integer.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The formula is:
!
!      INTVAL = Product ( 1 <= I <= PRIME_NUM ) PRIME(I)**NPOWER(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INTVAL, the integer to be factored.
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the number of prime factors for
!    which storage has been allocated.
!
!    Output, integer ( kind = 4 ) NPOWER(PRIME_NUM), the powers of the primes.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) intcopy
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) npower(prime_num)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) prime

  if ( intval <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_FAC - Fatal error!'
    write ( *, '(a)' ) '  Input integer is not positive.'
    stop
  end if
!
!  Try dividing the remainder by each prime.
!
  intcopy = intval

  do i = 1, prime_num

    npower(i) = 0

    p = prime ( i )

    do while ( mod ( intcopy, p ) == 0 )
      npower(i) = npower(i) + 1
      intcopy = intcopy / p
    end do

  end do

  return
end
function i4_to_isbn ( i )

!*****************************************************************************80
!
!! I4_TO_ISBN converts an I4 to an ISBN digit.
!
!  Discussion:
!
!    Only the integers 0 through 10 can be input.  The representation
!    of 10 is 'X'.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Book Industry Study Group,
!    The Evolution in Product Identification:
!    Sunrise 2005 and the ISBN-13,
!    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, an integer between 0 and 10.
!
!    Output, character I4_TO_ISBN, the ISBN character code of the integer.
!    If I is illegal, then I4_TO_ISBN is set to '?'.
!
  implicit none

  integer ( kind = 4 ) i
  character i4_to_isbn

       if ( i == 0 ) then
    i4_to_isbn = '0'
  else if ( i == 1 ) then
    i4_to_isbn = '1'
  else if ( i == 2 ) then
    i4_to_isbn = '2'
  else if ( i == 3 ) then
    i4_to_isbn = '3'
  else if ( i == 4 ) then
    i4_to_isbn = '4'
  else if ( i == 5 ) then
    i4_to_isbn = '5'
  else if ( i == 6 ) then
    i4_to_isbn = '6'
  else if ( i == 7 ) then
    i4_to_isbn = '7'
  else if ( i == 8 ) then
    i4_to_isbn = '8'
  else if ( i == 9 ) then
    i4_to_isbn = '9'
  else if ( i == 10 ) then
    i4_to_isbn = 'X'
  else
    i4_to_isbn = '?'
  end if

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2006
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
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine i4int_to_r8int ( imin, imax, i, rmin, rmax, r )

!*****************************************************************************80
!
!! I4INT_TO_R8INT maps an I4INT to an R8INT.
!
!  Discussion:
!
!    The formula used is:
!
!      R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IMIN, IMAX, the range.
!
!    Input, integer ( kind = 4 ) I, the integer to be converted.
!
!    Input, real ( kind = 8 ) RMIN, RMAX, the range.
!
!    Output, real ( kind = 8 ) R, the corresponding value in [RMIN,RMAX].
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) r
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmin

  if ( imax == imin ) then

    r = 0.5D+00 * ( rmin + rmax )

  else

    r = ( real ( imax - i,        kind = 8 ) * rmin   &
        + real (        i - imin, kind = 8 ) * rmax ) &
        / real ( imax     - imin, kind = 8 )

  end if

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_min ( n, a, amin )

!*****************************************************************************80
!
!! I4VEC_MIN computes the minimum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
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
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMIN, the value of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) amin

  amin = minval ( a(1:n) )

  return
end
subroutine i4vec_permute ( n, p, a )

!*****************************************************************************80
!
!! I4VEC_PERMUTE permutes an I4VEC in place.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    This routine permutes an array of integer "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = (   1,   2,   3,   4,   5 )
!
!    Output:
!
!      A    = (   2,   4,   5,   1,   3 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be permuted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_temp
  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine i4vec_uniform ( n, a, b, seed, x )

!*****************************************************************************80
!
!! I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) X(N), a vector of numbers between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
    r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
      +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
    value = nint ( r, kind = 4 )

    value = max ( value, min ( a, b ) )
    value = min ( value, max ( a, b ) )

    x(i) = value

  end do

  return
end
subroutine ij_next ( i, j, n )

!*****************************************************************************80
!
!! IJ_NEXT returns the next matrix index.
!
!  Discussion:
!
!    For N = 3, the sequence of indices returned is:
!
!      (1,1), (1,2), (1,3), (2,1), (2,2), (2,3), (3,1), (3,2), (3,3), (0,0).
!
!    Note that once the value (N,N) is returned, the next value returned
!    will be (0,0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On input, the current pair of 
!    indices.  On output, the next pair of indices.  If either index is illegal
!    on input, the output value of (I,J) will be (1,1).
!
!    Input, integer ( kind = 4 ) N, the maximum value for I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n

  if ( n < 1 ) then
    i = 0
    j = 0
    return
  end if

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then
    i = 1
    j = 1
    return
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n ) then
    i = i + 1
    j = 1
  else
    i = 0
    j = 0
  end if

  return
end
subroutine ij_next_gt ( i, j, n )

!*****************************************************************************80
!
!! IJ_NEXT_GT returns the next matrix index, with the constraint that I < J.
!
!  Discussion:
!
!    For N = 5, the sequence of indices returned is:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5), (4,5).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On input, the current pair of 
!    indices.  On output, the next pair of indices.  If either index is illegal 
!    on input, the output value of (I,J) will be (1,2).
!
!    Input, integer ( kind = 4 ) N, the maximum value for I and J.
!    A value of N less than 2 is nonsense.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n

  if ( n < 2 ) then
    i = 0
    j = 0
    return
  end if

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j .or. j <= i ) then
    i = 1
    j = 2
    return
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n - 1 ) then
    i = i + 1
    j = i + 1
  else
    i = 0
    j = 0
  end if

  return
end
subroutine index_box2_next_2d ( n1, n2, ic, jc, i, j, more )

!*****************************************************************************80
!
!! INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
!
!  Discussion:
!
!    The box has center at (IC,JC), and has half-widths N1 and N2.
!    The indices are exactly those which are between (IC-N1,JC-N2) and
!    (IC+N1,JC+N2) with the property that at least one of I and J
!    is an "extreme" value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the half-widths of the box, that is, 
!    the maximum distance allowed between (IC,JC) and (I,J).
!
!    Input, integer  ( kind = 4 )IC, JC, the central cell of the box.
!
!    Input/output, integer ( kind = 4 ) I, J.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I and J are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  logical more
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 ) then
    more = .false.
    return
  end if
!
!  Increment J.
!
  j = j + 1
!
!  Check J.
!
  if ( jc + n2 < j ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. ( i == ic - n1 .or. i == ic + n1 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, k, more )

!*****************************************************************************80
!
!! INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
!
!  Discussion:
!
!    The box has a central cell of (IC,JC,KC), with a half widths of
!    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3) and
!    (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J, and K
!    is an "extreme" value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the "half widths" of the box, 
!    that is, the maximum distances from the central cell allowed for I, J 
!    and K.
!
!    Input, integer ( kind = 4 ) IC, JC, KC, the central cell of the box.
!
!    Input/output, integer ( kind = 4 ) I, J, K.  On input, the previous index 
!    set.  On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I, J, and K are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  logical more
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    k = kc - n3
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 .and. k == kc + n3 ) then
    more = .false.
    return
  end if
!
!  Increment K.
!
  k = k + 1
!
!  Check K.
!
  if ( kc + n3 < k ) then
    k = kc - n3
    j = j + 1
  else if ( k < kc + n3 .and. &
    ( i == ic - n1 .or. &
      i == ic + n1 .or. &
      j == jc - n2 .or. &
      j == jc + n2 ) ) then
    return
  else
    k = kc + n3
    return
  end if
!
!  Check J.
!
  if ( jc + n2 < j ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. &
    ( i == ic - n1 .or. &
      i == ic + n1 .or. &
      k == kc - n3 .or. &
      k == kc + n3 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
function index1_col ( i_min, i, i_max, index_min )

!*****************************************************************************80
!
!! INDEX1_COL indexes a 1D vector by columns.
!
!  Discussion:
!
!    This 1D routine is provided primarily for analogy.
!    Moreover, in 1D there is no difference between row and column indexing.
!
!    Index       Element
!    ---------   --------
!    INDEX_MIN   I_MIN
!    INDEX1_COL  I
!   (INDEX_MAX)  I_MAX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for the first index,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) INDEX_MIN, the index of element I_MIN.
!    Typically, this is 0 or 1.
!
!    Output, integer ( kind = 4 ) INDEX1_COL, the index of element I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) index_min
  integer ( kind = 4 ) index1_col

  index1_col = index_min + ( i - i_min )

  return
end
function index1_row ( i_min, i, i_max, index_min )

!*****************************************************************************80
!
!! INDEX1_ROW indexes a 1D vector by rows.
!
!  Discussion:
!
!    This 1D routine is provided primarily for analogy.
!    Moreover, in 1D there is no difference between row and column indexing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for the first index,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) INDEX_MIN, the index of element I_MIN.
!    Typically, this is 0 or 1.
!
!    Output, integer ( kind = 4 ) INDEX1_ROW, the index of element I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) index_min
  integer ( kind = 4 ) index1_row

  index1_row = index_min + ( i - i_min )

  return
end
function index2_col ( i_min, i, i_max, j_min, j, j_max, index_min )

!*****************************************************************************80
!
!! INDEX2_COL indexes a 2D array by columns.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
!    and increasing the row index first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) INDEX_MIN, the index of element (I_MIN,J_MIN).
!    Typically, this is 0 or 1.
!
!    Output, integer ( kind = 4 ) INDEX2_COL, the index of element (I,J).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) index_min
  integer ( kind = 4 ) index2_col
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min

  index2_col = index_min &
             + (         i - i_min ) &
             + ( i_max + 1 - i_min ) * ( j - j_min )

  return
end
function index2_row ( i_min, i, i_max, j_min, j, j_max, index_min )

!*****************************************************************************80
!
!! INDEX2_ROW indexes a 2D array by row.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
!    and increasing the column index first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) INDEX_MIN, the index of element (I_MIN,J_MIN).
!    Typically, this is 0 or 1.
!
!    Output, integer ( kind = 4 ) INDEX2_ROW, the index of element (I,J).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) index_min
  integer ( kind = 4 ) index2_row
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min

  index2_row = index_min &
             +                         ( j - j_min ) &
             + ( i - i_min ) * ( j_max + 1 - j_min )

  return
end
function index3_col ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, &
  index_min )

!*****************************************************************************80
!
!! INDEX3_COL indexes a 3D array by columns.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry (I_MIN,J_MIN,K_MIN), 
!    and increasing the row index first, then the column index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) K_MIN, K, K_MAX, for plane indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) INDEX_MIN, the index of (I_MIN,J_MIN,K_MIN).
!    Typically, this is 0 or 1.
!
!    Output, integer ( kind = 4 ) INDEX3_COL, the index of element (I,J,K).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) index_min
  integer ( kind = 4 ) index3_col
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) k_min

  index3_col = index_min &
             + (         i - i_min ) &
             + ( i_max + 1 - i_min ) * (         j - j_min ) *  &
             + ( i_max + 1 - i_min ) * ( j_max + 1 - j_min ) * ( k - k_min )

  return
end
function index3_row ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, &
  index_min )

!*****************************************************************************80
!
!! INDEX3_ROW indexes a 3D array by rows.
!
!  Discussion:
!
!    When we say "by rows", we really just mean that entries of the array are 
!    indexed starting at entry (I_MIN,J_MIN,K_MIN), and the increasing the LAST
!    index first, then the next-to-the-last, and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) K_MIN, K, K_MAX, for plane indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) INDEX_MIN, the index of (I_MIN,J_MIN,K_MIN).
!    Typically, this is 0 or 1.
!
!    Output, integer ( kind = 4 ) INDEX3_ROW, the index of element (I,J,K).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) index_min
  integer ( kind = 4 ) index3_row
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) k_min

  index3_row = index_min &
             +                                                 ( k - k_min ) &
             +                         ( j - j_min ) * ( k_max + 1 - k_min ) &
             + ( i - i_min ) * ( j_max + 1 - j_min ) * ( k_max + 1 - k_min )

  return
end
function index4_col ( i1_min, i1, i1_max, i2_min, i2, i2_max, i3_min, i3, &
  i3_max, i4_min, i4, i4_max, index_min )

!*****************************************************************************80
!
!! INDEX4_COL indexes a 4D array by columns.
!
!  Discussion:
!
!    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
!    and increasing the initial index first, then the second, third and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1_MIN, I1, I1_MAX, for index 1,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I2_MIN, I2, I2_MAX, for index 2,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I3_MIN, I3, I3_MAX, for index 3,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I4_MIN, I4, I4_MAX, for index 4,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) INDEX_MIN, the index of 
!    (I1_MIN,I2_MIN,I3_MIN,I4_MIN).  Typically, this is 0 or 1.
!
!    Output, integer ( kind = 4 ) INDEX4_COL, the index of (I1,I2,I3,I4).
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i1_max
  integer ( kind = 4 ) i1_min
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2_max
  integer ( kind = 4 ) i2_min
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i3_max
  integer ( kind = 4 ) i3_min
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_max
  integer ( kind = 4 ) i4_min
  integer ( kind = 4 ) index_min
  integer ( kind = 4 ) index4_col

  index4_col = index_min &
    + (         i1 - i1_min ) &
    + ( i1_max + 1 - i1_min ) * (         i2 - i2_min ) &
    + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) * (         i3 - i3_min ) &
    + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) * ( i3_max + 1 - i3_min ) &
      (         i4 - i4_min )

  return
end
function index4_row ( i1_min, i1, i1_max, i2_min, i2, i2_max, i3_min, i3, &
  i3_max, i4_min, i4, i4_max, index_min )

!*****************************************************************************80
!
!! INDEX4_ROW indexes a 4D array by rows.
!
!  Discussion:
!
!    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
!    and increasing the last index, then the next to last, and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1_MIN, I1, I1_MAX, for index 1,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I2_MIN, I2, I2_MAX, for index 2,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I3_MIN, I3, I3_MAX, for index 3,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I4_MIN, I4, I4_MAX, for index 4,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) INDEX_MIN, the index of element 
!    (I1_MIN,I2_MIN,I3_MIN,I4_MIN).  Typically, this is 0 or 1.
!
!    Output, integer ( kind = 4 ) INDEX4_ROW, the index of (I1,I2,I3,I4).
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i1_max
  integer ( kind = 4 ) i1_min
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2_max
  integer ( kind = 4 ) i2_min
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i3_max
  integer ( kind = 4 ) i3_min
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_max
  integer ( kind = 4 ) i4_min
  integer ( kind = 4 ) index_min
  integer ( kind = 4 ) index4_row

  index4_row = index_min &
    +         ( i4 - i4_min ) &
    +                                                     ( i3 - i3_min ) &
    * ( i4_max + 1 - i4_min ) &
    +                           ( i2 - i2_min ) * ( i3_max + 1 - i3_min ) &
    * ( i4_max + 1 - i4_min ) &
    + ( i1 - i1_min ) * ( i2_max + 1 - i2_min ) * ( i3_max + 1 - i3_min ) &
    * ( i4_max + 1 - i4_min )

  return
end
function indexn_col ( n, i_min, i, i_max, index_min )

!*****************************************************************************80
!
!! INDEXN_COL indexes an ND array by columns.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry 
!      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
!    and increasing the first index up to I_MAX(1), 
!    then the second and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of indices.
!
!    Input, integer ( kind = 4 ) I_MIN(N), the minimum indices.
!
!    Input, integer ( kind = 4 ) I(N), the indices.
!
!    Input, integer ( kind = 4 ) I_MAX(N), for maximum indices.
!
!    Input, integer ( kind = 4 ) INDEX_MIN, the index of 
!    ( I_MIN(1), I_MIN(2),...,I_MIN(N) ).  Typically, this is 0 or 1.
!
!    Output, integer ( kind = 4 ) INDEXN_COL, the index of element I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i(n)
  integer ( kind = 4 ) i_max(n)
  integer ( kind = 4 ) i_min(n)
  integer ( kind = 4 ) index_min
  integer ( kind = 4 ) indexn_col
  integer ( kind = 4 ) j

  indexn_col = ( i(n) - i_min(n) )

  do j = n - 1, 1, - 1
    indexn_col = indexn_col * ( i_max(j) + 1 - i_min(j) ) + ( i(j) - i_min(j) )
  end do
  indexn_col = indexn_col + index_min

  return
end
function indexn_row ( n, i_min, i, i_max, index_min )

!*****************************************************************************80
!
!! INDEXN_ROW indexes an ND array by rows.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry 
!      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
!    and increasing the last index up to I_MAX(N), 
!    then the next-to-last and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of indices.
!
!    Input, integer ( kind = 4 ) I_MIN(N), the minimum indices.
!
!    Input, integer ( kind = 4 ) I(N), the indices.
!
!    Input, integer ( kind = 4 ) I_MAX(N), for maximum indices.
!
!    Input, integer ( kind = 4 ) INDEX_MIN, the index of 
!    ( I_MIN(1), I_MIN(2),...,I_MIN(N) ).  Typically, this is 0 or 1.
!
!    Output, integer ( kind = 4 ) INDEXN_ROW, the index of element I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i(n)
  integer ( kind = 4 ) i_max(n)
  integer ( kind = 4 ) i_min(n)
  integer ( kind = 4 ) index_min
  integer ( kind = 4 ) indexn_row
  integer ( kind = 4 ) j

  indexn_row = ( i(1) - i_min(1) )

  do j = 2, n
    indexn_row = indexn_row * ( i_max(j) + 1 - i_min(j) ) + ( i(j) - i_min(j) )
  end do
  indexn_row = indexn_row + index_min

  return
end
subroutine isbn_check ( isbn, check )

!*****************************************************************************80
!
!! ISBN_CHECK checks an ISBN code.
!
!  Discussion:
!
!    ISBN stands for International Standard Book Number.  A unique ISBN
!    is assigned to each new book.  The ISBN includes 10 digits.  There is
!    an initial digit, then a dash, then a set of digits which are a
!    code for the publisher, another digit, and then the check digit:
!
!      initial-publisher-book-check
!
!    The number of digits used for the publisher and book codes can vary,
!    but the check digit is always one digit, and the total number of
!    digits is always 10.
!
!    The check digit is interesting, because it is a way of trying to
!    make sure that an ISBN has not been incorrectly copied.  Specifically,
!    if the ISBN is correct, then its ten digits will satisfy
!
!       10 * A + 9 * B + 8 * C + 7 * D + 6 * E
!      + 5 * F * 4 * G * 3 * H + 2 * I +     J  = 0 mod 11.
!
!    Here, we've taken 'A' to represent the first digit and 'J' the
!    last (which is the check digit).  In order for the code to work,
!    the value of J must be allowed to be anything from 0 to 10.  In
!    cases where J works out to be 10, the special digit 'X' is used.
!    An 'X' digit can only occur in this last check-digit position
!    on an ISBN.
!
!  Example:
!
!    0-8493-9640-9
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Book Industry Study Group,
!    The Evolution in Product Identification:
!    Sunrise 2005 and the ISBN-13,
!    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
!
!  Parameters:
!
!    Input, character ( len = * ) ISBN, an ISBN code.
!
!    Output, integer ( kind = 4 ) CHECK, the value of the ISBN check sum.
!    If CHECK is zero, the ISBN code is legitimate.
!    If CHECK is -1, then the ISBN code is not legitimate because it does
!    not contain exactly 10 digits.  If CHECK is between 1 and 10, then
!    the ISBN code has the right number of digits, but at least one of
!    the digits is incorrect.
!
  implicit none

  character c
  logical ch_is_digit
  integer ( kind = 4 ) check
  integer ( kind = 4 ) digit(10)
  integer ( kind = 4 ) i
  character ( len = * ) isbn
  integer ( kind = 4 ) isbn_to_i4
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) num_digit
!
!  Determine how many digits have been supplied.
!
  lenc = len_trim ( isbn )

  i = 0
  num_digit = 0

  do

    i = i + 1

    if ( lenc < i ) then
      exit
    end if

    c = isbn(i:i)

    if ( ch_is_digit ( c ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i4 ( c )

    else if ( ( num_digit == 9 .and. isbn(i:i) == 'X' ) .or. &
              ( num_digit == 9 .and. isbn(i:i) == 'x' ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i4 ( c )

    end if

    if ( 10 <= num_digit ) then
      exit
    end if

  end do
!
!  If we didn't get exactly 10 digits, return with an error.
!
  if ( num_digit /= 10 ) then
    check = -1
    return
  end if
!
!  Compute the checksum.
!
  check = 0
  do i = 1, 10
    check = check + ( 11 - i ) * digit(i)
  end do

  check = mod ( check, 11 )

  return
end
subroutine isbn_fill ( isbn )

!*****************************************************************************80
!
!! ISBN_FILL fills in a missing digit in an ISBN code.
!
!  Example:
!
!    Input:
!
!      0-8493-9?40-9
!
!    Output:
!
!      0-8493-9640-9
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Book Industry Study Group,
!    The Evolution in Product Identification:
!    Sunrise 2005 and the ISBN-13,
!    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
!
!  Parameters:
!
!    Input/output, character ( len = * ) ISBN, a partial ISBN code.  On input,
!    a single digit has been replaced by the character '?', signifying
!    that that digit is missing.  The routine replaces the question
!    mark by the correct digit.
!
  implicit none

  character c
  logical ch_is_digit
  integer ( kind = 4 ) check
  integer ( kind = 4 ) digit(10)
  integer ( kind = 4 ) digit_pos
  integer ( kind = 4 ) i
  character i4_to_isbn
  character ( len = * ) isbn
  integer ( kind = 4 ) isbn_pos
  integer ( kind = 4 ) isbn_to_i4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) num_digit

  lenc = len_trim ( isbn )

  i = 0
  isbn_pos = -1
  digit_pos = -1
  num_digit = 0

  do

    i = i + 1

    if ( lenc < i ) then
      exit
    end if

    c = isbn(i:i)

    if ( ch_is_digit ( c ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i4 ( c )

    else if ( ( num_digit == 9 .and. isbn(i:i) == 'X' ) .or. &
              ( num_digit == 9 .and. isbn(i:i) == 'x' ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i4 ( c )

    else if ( c == '?' ) then

      if ( isbn_pos == -1 ) then

        num_digit = num_digit + 1
        digit(num_digit) = 0
        digit_pos = num_digit
        isbn_pos = i

      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ISBN_FILL - Fatal error!'
        write ( *, '(a)' ) '  Only one question mark is allowed!'
        return
      end if

    end if

    if ( 10 <= num_digit ) then
      exit
    end if

  end do

  if ( num_digit /= 10 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ISBN_FILL - Fatal error!'
    write ( *, '(a)' ) '  The input ISBN code did not have 10 digits.'
    return
  end if

  if ( isbn_pos == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ISBN_FILL - Fatal error!'
    write ( *, '(a)' ) '  A question mark is required!'
    return
  end if

  check = 0
  do i = 1, 10
    check = check + ( 11 - i ) * digit(i)
  end do

  check = mod ( check, 11 )

  if ( check == 0 ) then

    k = 0
!
!  Need to solve the modular equation:
!
!    A * X = B mod C
!
!  Below is a stupid way.  One day I will come back and fix this up.
!
  else

    do i = 1, 10
      j = ( 11 - digit_pos ) * i + check
      if ( mod ( j, 11 ) == 0 ) then
        k = i
      end if
    end do

  end if

  isbn(isbn_pos:isbn_pos) = i4_to_isbn ( k )

  return
end
function isbn_to_i4 ( c )

!*****************************************************************************80
!
!! ISBN_TO_I4 converts an ISBN character into an integer.
!
!  Discussion:
!
!    The characters '0' through '9' stand for themselves, but
!    the character 'X' or 'x' stands for 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Book Industry Study Group,
!    The Evolution in Product Identification:
!    Sunrise 2005 and the ISBN-13,
!    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
!
!  Parameters:
!
!    Input, character C, the ISBN character code to be converted.
!
!    Output, integer ( kind = 4 ) ISBN_TO_I4, the numeric value of the character
!    code, between 0 and 10.  This value is returned as -1 if C is
!    not a valid character code.
!
  implicit none

  character c
  integer ( kind = 4 ) isbn_to_i4

       if ( c == '0' ) then
    isbn_to_i4 = 0
  else if ( c == '1' ) then
    isbn_to_i4 = 1
  else if ( c == '2' ) then
    isbn_to_i4 = 2
  else if ( c == '3' ) then
    isbn_to_i4 = 3
  else if ( c == '4' ) then
    isbn_to_i4 = 4
  else if ( c == '5' ) then
    isbn_to_i4 = 5
  else if ( c == '6' ) then
    isbn_to_i4 = 6
  else if ( c == '7' ) then
    isbn_to_i4 = 7
  else if ( c == '8' ) then
    isbn_to_i4 = 8
  else if ( c == '9' ) then
    isbn_to_i4 = 9
  else if ( c == 'X' .or. c == 'x' ) then
    isbn_to_i4 = 10
  else
    isbn_to_i4 = -1
  end if

  return
end
function iset2_compare ( x1, y1, x2, y2 )

!*****************************************************************************80
!
!! ISET2_COMPARE compares two I2 sets.
!
!  Discussion:
!
!    The I2 set (X1,Y1) < (X2,Y2) if
!
!      min ( X1, Y1 ) < min ( X2, Y2 ) or
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) < max ( X2, Y2 )
!
!    The I2 set (X1,Y1) = (X2,Y2) if
!
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) = max ( X2, Y2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X1, Y1, the first I2 set.
!
!    Input, integer ( kind = 4 ) X2, Y2, the second I2 set.
!
!    Output, integer ( kind = 4 ) ISET2_COMPARE: 
!    -1, (X1,Y1) < (X2,Y2);
!     0, (X1,Y1) = (X2,Y2);
!    +1, (X1,Y1) > (X2,Y2).
!
  implicit none

  integer ( kind = 4 ) a1
  integer ( kind = 4 ) a2
  integer ( kind = 4 ) b1
  integer ( kind = 4 ) b2
  integer ( kind = 4 ) iset2_compare
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x1
  integer ( kind = 4 ) x2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  a1 = min ( x1, y1 )
  b1 = max ( x1, y1 )

  a2 = min ( x2, y2 )
  b2 = max ( x2, y2 )

  if ( a1 < a2 ) then
    value = -1
  else if ( a2 < a1 ) then
    value = +1
  else if ( b1 < b2 ) then
    value = -1
  else if ( b2 < b1 ) then
    value = +1
  else
    value = 0
  end if

  iset2_compare = value

  return
end
subroutine iset2_index_insert_unique ( n_max, n, x, y, indx, &
  xval, yval, ival, ierror )

!*****************************************************************************80
!
!! ISET2_INDEX_INSERT_UNIQUE inserts unique values into an indexed sorted list.
!
!  Discussion:
!
!    If the input value does not occur in the list, then N, X, Y and INDX
!    are updated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_MAX, the maximum size of the list.
!
!    Input/output, integer ( kind = 4 ) N, the size of the list.
!
!    Input/output, integer ( kind = 4 ) X(N), Y(N), the list of I2 sets.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, YVAL, the value to be inserted if it is
!    not already in the list.
!
!    Output, integer ( kind = 4 ) IVAL, the index in INDX corresponding to the
!    value XVAL, YVAL.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 if an error 
!    occurred.
!
  implicit none

  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ) xval
  integer ( kind = 4 ) y(n_max)
  integer ( kind = 4 ) yval

  ierror = 0

  if ( n <= 0 ) then

    if ( n_max <= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ISET2_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = min ( xval, yval )
    y(1) = max ( xval, yval )
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL ) already occur in the list?
!
  call iset2_index_search ( n_max, n, x, y, indx, xval, yval, &
    less, equal, more )

  if ( equal == 0 ) then

    if ( n_max <= n ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ISET2_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = min ( xval, yval )
    y(n+1) = max ( xval, yval )
    ival = more
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = equal

  end if

  return
end
subroutine iset2_index_search ( n_max, n, x, y, indx, xval, yval, &
  less, equal, more )

!*****************************************************************************80
!
!! ISET2_INDEX_SEARCH searches for an I2 set value in an indexed sorted list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_MAX, the maximum size of the list.
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, integer ( kind = 4 ) X(N), Y(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, YVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) LESS, EQUAL, MORE, the indexes in INDX of the
!    list entries that are just less than, equal to, and just greater
!    than the test value.  If the test value does not occur in the list,
!    then EQUAL is zero.  If the test value is the minimum entry of the
!    list, then LESS is 0.  If the test value is the greatest entry of
!    the list, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) compare
  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) iset2_compare
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ) xhi
  integer ( kind = 4 ) xlo
  integer ( kind = 4 ) xmid
  integer ( kind = 4 ) xval
  integer ( kind = 4 ) y(n_max)
  integer ( kind = 4 ) yhi
  integer ( kind = 4 ) ylo
  integer ( kind = 4 ) ymid
  integer ( kind = 4 ) yval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))

  compare = iset2_compare ( xval, yval, xlo, ylo )

  if ( compare == -1 ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( compare == 0 ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  compare = iset2_compare ( xval, yval, xhi, yhi )

  if ( compare == +1 ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( compare == 0 ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))

    compare = iset2_compare ( xval, yval, xmid, ymid )

    if ( compare == 0 ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( compare == -1 ) then
      hi = mid
    else if ( compare == +1 ) then
      lo = mid
    end if

  end do

  return
end
function lcm_12n ( n )

!*****************************************************************************80
!
!! LCM_12N computes the least common multiple of the integers 1 through N.
!
!  Example:
!
!    N    LCM_12N
!
!    1          1
!    2          2
!    3          3
!    4         12
!    5         60
!    6         60
!    7        420
!    8        840
!    9       2520
!   10       2520
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the value of N.
!
!    Output, integer ( kind = 4 ) LCM_12N, the least common multiple of the 
!    integers 1 to N.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imult
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lcm_12n
  integer ( kind = 4 ) n

  lcm_12n = 1

  do i = 2, n

    imult = i

    do j = 1, i - 1

      if ( mod ( imult, ( i - j ) ) == 0 ) then
        imult = imult / ( i - j )
      end if

    end do

    lcm_12n = lcm_12n * imult

  end do

  return
end
subroutine lmat_print ( m, n, a, title )

!*****************************************************************************80
!
!! LMAT_PRINT prints an LMAT.
!
!  Discussion:
!
!    An LMAT is an array of L values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2011
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
!    Input, logical A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical a(m,n)
  character ( len = * ) title

  call lmat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine lmat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! LMAT_PRINT_SOME prints some of an LMAT.
!
!  Discussion:
!
!    An LMAT is an array of L values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, logical A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 35
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical a(m,n)
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
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    if ( 100 <= j2hi ) then
      do j = j2lo, j2hi
        j2 = j + 1 - j2lo
        write ( ctemp(j2), '(1x,i1)' ) j / 100
      end do
      write ( *, '(''      '',35a2)' ) ctemp(1:inc)
    end if

    if ( 10 <= j2hi ) then
      do j = j2lo, j2hi
        j2 = j + 1 - j2lo
        write ( ctemp(j2), '(1x,i1)' ) mod ( j / 10, 10 )
      end do
      write ( *, '(''      '',35a2)' ) ctemp(1:inc)
    end if

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(1x,i1)' ) mod ( j, 10 )
    end do
    write ( *, '(''  Col '',35a2)' ) ctemp(1:inc)

    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
      write ( *, '(i5,a1,35(1x,l1))' ) i, ':', a(i,j2lo:j2hi)
    end do

  end do

  return
end
subroutine lmat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! LMAT_TRANSPOSE_PRINT prints an LMAT, transposed.
!
!  Discussion:
!
!    An LMAT is an array of L values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, logical A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical a(m,n)
  character ( len = * ) title

  call lmat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine lmat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! LMAT_TRANSPOSE_PRINT_SOME prints some of an LMAT, transposed.
!
!  Discussion:
!
!    An LMAT is an array of L values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, logical A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 35
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical a(m,n)
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
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    if ( 100 <= i2hi ) then
      do i = i2lo, i2hi
        i2 = i + 1 - i2lo
        write ( ctemp(i2), '(1x,i1)' ) i / 100
      end do
      write ( *, '(''      '',35a2)' ) ctemp(1:inc)
    end if

    if ( 10 <= i2hi ) then
      do i = i2lo, i2hi
        i2 = i + 1 - i2lo
        write ( ctemp(i2), '(1x,i1)' ) mod ( i / 10, 10 )
      end do
      write ( *, '(''      '',35a2)' ) ctemp(1:inc)
    end if

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(1x,i1)' ) mod ( i, 10 )
    end do
    write ( *, '(''  Row '',35a2)' ) ctemp(1:inc)

    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi
      write ( *, '(i5,a1,35(1x,l1))' ) j, ':', a(i2lo:i2hi,j)
    end do

  end do

  return
end
subroutine luhn_check ( digit_num, digit, check_sum )

!*****************************************************************************80
!
!! LUHN_CHECK computes the Luhn checksum for a string of digits.
!
!  Discussion:
!
!    To compute the Luhn checksum, begin at the end of the string, and double
!    every other digit.  If a doubled digit is greater than 9, subtract 9.
!    Then sum the digits to get CHECK_SUM.
!
!    If mod ( CHECK_SUM, 10 ) = 0 the digit sequence is accepted.
!
!    The Luhn check sum will detect any single digit error, as well as
!    most errors involving a transposition of adjacent digits; it cannot
!    detect the transposition of the adjacent digits 0 and 9, however.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT_NUM, the number of digits.
!
!    Input, integer ( kind = 4 ) DIGIT(DIGIT_NUM), the string of digits.
!    Normally, these are true digits, that is, values between 0 and 9.
!
!    Output, integer ( kind = 4 ) CHECK_SUM, the check sum, which
!    should be divisible by 10 if the digit sequence is correct.
!
  implicit none

  integer ( kind = 4 ) digit_num

  integer ( kind = 4 ) check_sum
  integer ( kind = 4 ) digit(digit_num)
  integer ( kind = 4 ) digit_copy(digit_num)

  digit_copy(1:digit_num) = digit(1:digit_num)

  digit_copy(digit_num-1:1:-2) = 2 * digit_copy(digit_num-1:1:-2)

  where ( 9 < digit_copy )
    digit_copy = digit_copy - 9
  end where

  check_sum = sum ( digit_copy(1:digit_num) )

  return
end
subroutine lvec_print ( n, a, title )

!*****************************************************************************80
!
!! LVEC_PRINT prints an LVEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, logical ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 )  n

  logical a(n)
  character ( len = * ) title

  call lvec_print_some ( n, a, 1, n, title )

  return
end
subroutine lvec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! LVEC_PRINT_SOME prints "some" of an LVEC.
!
!  Discussion:
!
!    An LVEC is a vector of logical values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, logical A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices to 
!    print. The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  logical a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,a,1x,l1)' ) i, ':', a(i)
  end do

  return
end
function pause_input ( )

!*****************************************************************************80
!
!! PAUSE_INPUT waits until an input character is entered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character PAUSE_INPUT, the character that was entered.
!
  implicit none

  integer ( kind = 4 ) ios
  character pause_input

  write ( *, '(a)' ) 'Press RETURN to continue.'
  read ( *, '(a)', iostat = ios ) pause_input

  return
end
subroutine perm_check ( n, p, base, ierror )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from BASE to
!    to BASE+N-1 occurs among the N entries of the permutation.
!
!    Set the input quantity BASE to 0, if P is a 0-based permutation,
!    or to 1 if P is a 1-based permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the array to check.
!
!    Input, integer ( kind = 4 ) BASE, the index base.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) find
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seek

  ierror = 0

  do seek = base, base + n - 1

    ierror = 1

    do find = 1, n
      if ( p(find) == seek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a proper permutation.'
      stop
    end if

  end do

  return
end
subroutine perm_cycle ( n, iopt, p, isgn, ncycle )

!*****************************************************************************80
!
!! PERM_CYCLE analyzes a permutation.
!
!  Discussion:
!
!    The routine will count cycles, find the sign of a permutation,
!    and tag a permutation.
!
!  Example:
!
!    Input:
!
!      N = 9
!      IOPT = 1
!      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NCYCLE = 3
!      ISGN = +1
!      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 July 2000
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
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input, integer ( kind = 4 ) IOPT, requests tagging.
!    0, the permutation will not be tagged.
!    1, the permutation will be tagged.
!
!    Input/output, integer ( kind = 4 ) P(N).  On input, P describes a
!    permutation, in the sense that entry I is to be moved to P(I).
!    If IOPT = 0, then P will not be changed by this routine.
!    If IOPT = 1, then on output, P will be "tagged".  That is,
!    one element of every cycle in P will be negated.  In this way,
!    a user can traverse a cycle by starting at any entry I1 of P
!    which is negative, moving to I2 = ABS(P(I1)), then to
!    P(I2), and so on, until returning to I1.
!
!    Output, integer ( kind = 4 ) ISGN, the "sign" of the permutation, which is
!    +1 if the permutation is even, -1 if odd.  Every permutation
!    may be produced by a certain number of pairwise switches.
!    If the number of switches is even, the permutation itself is
!    called even.
!
!    Output, integer ( kind = 4 ) NCYCLE, the number of cycles in the 
!    permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) is
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) ncycle
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_CYCLE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if

  is = 1
  ncycle = n

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      ncycle = ncycle - 1
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    if ( iopt /= 0 ) then
      is = - i4_sign ( p(i) )
    end if

    p(i) = is * abs ( p(i) )

  end do

  isgn = 1 - 2 * mod ( n - ncycle, 2 )

  return
end
subroutine perm_free ( npart, ipart, nfree, ifree )

!*****************************************************************************80
!
!! PERM_FREE reports the number of unused items in a partial permutation.
!
!  Discussion:
!
!    It is assumed that the N objects being permuted are the integers
!    from 1 to N, and that IPART contains a "partial" permutation, that
!    is, the NPART entries of IPART represent the beginning of a
!    permutation of all N items.
!
!    The routine returns in IFREE the items that have not been used yet.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPART, the number of entries in IPART.  
!    NPART may be 0.
!
!    Input, integer ( kind = 4 ) IPART(NPART), the partial permutation, which 
!    should contain, at most once, some of the integers between 1 and
!    NPART+NFREE.
!
!    Input, integer ( kind = 4 ) NFREE, the number of integers that have not 
!    been used in IPART.  This is simply N - NPART.  NFREE may be zero.
!
!    Output, integer ( kind = 4 ) IFREE(NFREE), the integers between 1 and 
!    NPART+NFREE that were not used in IPART.
!
  implicit none

  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nfree)
  integer ( kind = 4 ) ipart(npart)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) match
  integer ( kind = 4 ) n

  n = npart + nfree

  if ( npart < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NPART < 0.'
    write ( *, '(a,i8)' ) '  NPART = ', npart
    stop

  else if ( npart == 0 ) then

    call i4vec_indicator ( n, ifree )

  else if ( nfree < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NFREE < 0.'
    write ( *, '(a,i8)' ) '  NFREE = ', nfree
    stop

  else if ( nfree == 0 ) then

    return

  else

    k = 0

    do i = 1, n

      match = 0

      do j = 1, npart
        if ( ipart(j) == i ) then
          match = j
          exit
        end if
      end do

      if ( match == 0 ) then

        k = k + 1

        if ( nfree < k ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)'    ) 'PERM_FREE - Fatal error!'
          write ( *, '(a)'    ) '  The partial permutation is illegal.'
          write ( *, '(a)'    ) '  It should contain, at most once, some of'
          write ( *, '(a,i8)' ) '  the integers between 1 and N = ', n
          write ( *, '(a)'    ) '  The number of integers that have not'
          write ( *, '(a,i8)' ) '  been used is at least K = ', k
          write ( *, '(a,i8)' ) '  This should be exactly NFREE = ', nfree
          call i4vec_print ( npart, ipart, '  The partial permutation:' )
          stop
        end if

        ifree(k) = i

      end if

    end do

  end if

  return
end
subroutine perm_inverse ( n, p )

!*****************************************************************************80
!
!! PERM_INVERSE inverts a permutation "in place".
!
!  Discussion:
!
!    This algorithm assumes that the entries in the permutation vector are
!    strictly positive.  In particular, the value 0 must not occur.
!
!    When necessary, this function shifts the data temporarily so that
!    this requirement is satisfied.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard 
!    index form.  On output, P describes the inverse permutation
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) p_min

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop
  end if
!
!  Find the least value, and shift data so it begins at 1.
!
  call i4vec_min ( n, p, p_min )
  base = 1
  p(1:n) = p(1:n) - p_min + base
!
!  Check the permutation.
!
  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Invert the permutation.
!
  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    is = - i4_sign ( p(i) )
    p(i) = is * abs ( p(i) )

  end do

  do i = 1, n

    i1 = - p(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = p(i1)
        p(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do
!
!  Reverse the shift.
!
  p(1:n) = p(1:n) + p_min - base

  return
end
subroutine perm_next ( n, p, more, even )

!*****************************************************************************80
!
!! PERM_NEXT computes all of the permutations on N objects, one at a time.
!
!  Discussion:
!
!    If this routine is called with MORE = TRUE, any
!    permutation in P, and EVEN = TRUE, then the successor of the input
!    permutation will be produced, unless P is the last permutation
!    on N letters, in which case P(1) will be set to 0 on return.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt
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
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N).
!    On input, P contains the previous permutation.
!    On output, P contains the next permutation.
!
!    Input/output, logical MORE.
!    On input, MORE = FALSE means this is the first call.
!    On output, MORE = FALSE means there are no more permutations.
!
!    Output, logical EVEN, is TRUE if the output permutation is even.
!
  implicit none

  integer ( kind = 4 ) n

  logical even
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) id
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) p(n)

  if ( .not. more ) then

    call i4vec_indicator ( n, p )

    more = .true.
    even = .true.

    if ( n == 1 ) then
      more = .false.
      return
    end if

    if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
      return
    end if

    do i = 1, n - 3
      if ( p(i+1) /= p(i) + 1 ) then
        return
      end if
    end do

    more = .false.

  else

    if ( n == 1 ) then
      p(1) = 0
      more = .false.
      return
    end if

    if ( even ) then

      ia = p(1)
      p(1) = p(2)
      p(2) = ia
      even = .false.

      if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
        return
      end if

      do i = 1, n - 3
        if ( p(i+1) /= p(i) + 1 ) then
          return
        end if
      end do

      more = .false.
      return

    else

      is = 0
      more = .false.

      do i1 = 2, n

        ia = p(i1)
        i = i1 - 1
        id = 0

        do j = 1, i
          if ( i1 < p(j) ) then
            id = id + 1
          end if
        end do

        is = id + is

        if ( id /= i * mod ( is, 2 ) ) then
          more = .true.
          exit
        end if

      end do

      if ( .not. more ) then
        p(1) = 0
        return
      end if

    end if

    m = mod ( is + 1, 2 ) * ( n + 1 )

    do j = 1, i

      if ( i4_sign ( p(j) - ia ) /= i4_sign ( p(j) - m ) ) then
        m = p(j)
        l = j
      end if

    end do

    p(l) = ia
    p(i1) = m
    even = .true.

  end if

  return
end
subroutine perm_print ( n, p, title )

!*****************************************************************************80
!
!! PERM_PRINT prints a permutation.
!
!  Example:
!
!    Input:
!
!      P = 7 2 4 1 5 3 6
!
!    Printed output:
!
!      1 2 3 4 5 6 7
!      7 2 4 1 5 3 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer which is to be partitioned.
!
!    Input, integer ( kind = 4 ) P(N), the permutation to be converted.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), parameter :: inc = 20
  integer ( kind = 4 ) p(n)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  do ilo = 1, n, inc
    ihi = min ( n, ilo + inc - 1 )
    write ( *, '(a)' ) ' '
    write ( *, '(20i4)' ) ( i,i = ilo, ihi )
    write ( *, '(20i4)' ) p(ilo:ihi)
  end do

  return
end
subroutine perm_uniform ( n, base, seed, p )

!*****************************************************************************80
!
!! PERM_UNIFORM selects a random permutation of N objects.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2008
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
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input, integer ( kind = 4 ) BASE, is 0 for a 0-based permutation and 1 for 
!    a 1-based permutation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) P(N), the permutation.  P(I) is the "new"
!    location of the object originally at I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  do i = 1, n
    p(i) = ( i - 1 ) + base
  end do

  do i = 1, n
    j = i4_uniform ( i, n, seed )
    k    = p(i)
    p(i) = p(j)
    p(j) = k
  end do

  return
end
function pounds_to_kilograms ( lb )

!*****************************************************************************80
!
!! POUNDS_TO_KILOGRAMS converts a measurement in pounds to kilograms.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LB, the weight in pounds.
!
!    Output, real ( kind = 8 ) POUNDS_TO_KILOGRAMS, the corresponding 
!    weight in kilograms.
!
  implicit none

  real ( kind = 8 ) lb
  real ( kind = 8 ) pounds_to_kilograms

  pounds_to_kilograms = 0.4535924D+00 * lb

  return
end
function prime ( n )

!*****************************************************************************80
!
!! PRIME returns any of the first PRIME_MAX prime numbers.
!
!  Discussion:
!
!    PRIME_MAX is 1600, and the largest prime stored is 13499.
!
!    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964, pages 870-873.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the desired prime number.
!    In general, is should be true that 0 <= N <= PRIME_MAX.
!    N = -1 returns PRIME_MAX, the index of the largest prime available.
!    N = 0 is legal, returning PRIME = 1.
!
!    Output, integer ( kind = 4 ) PRIME, the N-th prime.  If N is out of range, 
!    PRIME is returned as -1.
!
  implicit none

  integer ( kind = 4 ), parameter :: prime_max = 1600

  integer ( kind = 4 ), save :: icall = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( prime_max ) :: npvec
  integer ( kind = 4 ) prime

  if ( icall == 0 ) then

    icall = 1

    npvec(1:100) = (/ &
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /)

    npvec(101:200) = (/ &
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /)

    npvec(201:300) = (/ &
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /)

    npvec(301:400) = (/ &
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /)

    npvec(401:500) = (/ &
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /)

    npvec(501:600) = (/ &
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /)

    npvec(601:700) = (/ &
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /)

    npvec(701:800) = (/ &
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /)

    npvec(801:900) = (/ &
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /)

    npvec(901:1000) = (/ &
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

    npvec(1001:1100) = (/ &
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /)

    npvec(1101:1200) = (/ &
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /)

    npvec(1201:1300) = (/ &
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /)

    npvec(1301:1400) = (/ &
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949, &
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /)

    npvec(1401:1500) = (/ &
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /)

   npvec(1501:1600) = (/ &
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, &
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, &
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, &
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, &
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, &
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, &
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, &
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, &
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, &
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /)

  end if

  if ( n == -1 ) then
    prime = prime_max
  else if ( n == 0 ) then
    prime = 1
  else if ( n <= prime_max ) then
    prime = npvec(n)
  else
    prime = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIME - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i8)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop
  end if

  return
end
function prime_ge ( n )

!*****************************************************************************80
!
!! PRIME_GE returns the smallest prime greater than or equal to N.
!
!  Example:
!
!    N     PRIME_GE
!
!    -10    2
!      1    2
!      2    2
!      3    3
!      4    5
!      5    5
!      6    7
!      7    7
!      8   11
!      9   11
!     10   11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number to be bounded.
!
!    Output, integer ( kind = 4 ) PRIME_GE, the smallest prime number that is 
!    greater than or equal to N.  However, if N is larger than the largest
!    prime stored, then PRIME_GE is returned as -1.
!
  implicit none

  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  integer ( kind = 4 ) i_mid
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p_hi
  integer ( kind = 4 ) p_lo
  integer ( kind = 4 ) p_mid
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) prime_ge

  if ( n <= 2 ) then
    prime_ge = 2
    return
  end if

  i_lo = 1
  p_lo = prime(i_lo)
  i_hi = prime(-1)
  p_hi = prime(i_hi)

  if ( p_hi < n ) then
    prime_ge = - p_hi
    return
  end if

  do

    if ( i_lo + 1 == i_hi ) then
      prime_ge = p_hi
      exit
    end if

    i_mid = ( i_lo + i_hi ) / 2
    p_mid = prime(i_mid)

    if ( p_mid < n ) then
      i_lo = i_mid
      p_lo = p_mid
    else if ( n <= p_mid ) then
      i_hi = i_mid
      p_hi = p_mid
    end if

  end do

  return
end
subroutine primer ( n, iprime )

!*****************************************************************************80
!
!! PRIMER computes the prime numbers up to a given limit.
!
!  Discussion:
!
!    PRIMER returns the results of its computations in the vector
!    IPRIME.  IPRIME(I) is -1 if the number I is not prime, and
!    1 if I is prime.
!
!    The algorithm is a simple-minded sieve of Eratosthenes, with
!    no attempt at efficiency.
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
!    Input, integer ( kind = 4 ) N, the dimension of IPRIME, and the maximum
!    value that will be considered.
!
!    Output, integer ( kind = 4 ) IPRIME(N), records the results for each 
!    integer.  IPRIME(I) = -1 if I is not prime, and IPRIME(I) = 1 if I is
!    prime.  By convention, IPRIME(1) will be set to -1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iprime(n)
  integer ( kind = 4 ) next
!
!  IPRIME(I) = 0 means we don't know if I is prime.
!
  iprime(1:n) = 0
!
!  By convention, 1 is not prime.
!
  iprime(1) = -1
  next = 1
!
!  Examine the integers in order.
!
  do next = 2, n

    if ( iprime(next) == 0 ) then
      iprime(next) = 1
      do i = 2 * next, n, next
        iprime(i) = -1
      end do
    end if

  end do

  return
end
function r8_log_10 ( x )

!*****************************************************************************80
!
!! R8_LOG_10 returns the logarithm base 10 of an R8.
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
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose base 2 logarithm is desired.
!    X should not be 0.
!
!    Output, real ( kind = 8 ) R8_LOG_10, the logarithm base 10 of the absolute
!    value of X.  It should be true that |X| = 10**R_LOG_10.
!
  implicit none

  real ( kind = 8 ) r8_log_10
  real ( kind = 8 ) x

  if ( x == 0.0D+00 ) then
    r8_log_10 = -huge ( x )
  else
    r8_log_10 = log10 ( abs ( x ) )
  end if

  return
end
function r8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
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
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

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
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
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
!    05 July 2006
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
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
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
!    12 September 2004
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
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
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
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
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

  real ( kind = 8 ) a(m,n)
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

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
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
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
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

  real ( kind = 8 ) a(m,n)
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
subroutine r8poly_degree ( na, a, degree )

!*****************************************************************************80
!
!! R8POLY_DEGREE returns the degree of a polynomial.
!
!  Discussion:
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:NA), the coefficients of the polynomials.
!
!    Output, integer ( kind = 4 ) DEGREE, the degree of A.
!
  implicit none

  integer ( kind = 4 ) na

  real ( kind = 8 ) a(0:na)
  integer ( kind = 4 ) degree

  degree = na

  do while ( 0 < degree )

    if ( a(degree) /= 0.0D+00 ) then
      return
    end if

    degree = degree - 1

  end do

  return
end
subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! R8POLY_PRINT prints out a polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X**N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mag
  integer ( kind = 4 ) n2
  character plus_minus
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  call r8poly_degree ( n, a, n2 )

  if ( n2 <= 0 ) then
    write ( *, '( ''  p(x) = 0'' )' )
    return
  end if

  if ( a(n2) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( 2 <= n2 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) &
      plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n2-1, 0, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine r8vec_indicator ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR sets an R8VEC to the indicator vector.
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
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine r8vec_mean ( n, a, mean )

!*****************************************************************************80
!
!! R8VEC_MEAN returns the mean of an R8VEC.
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
!    02 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector whose mean is desired.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) mean

  mean = sum ( a(1:n) ) / real ( n, kind = 8 )

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
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
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_variance ( n, a, variance )

!*****************************************************************************80
!
!! R8VEC_VARIANCE returns the variance of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The variance of a vector X of length N is defined as
!
!      mean ( X(1:n) ) = sum ( X(1:n) ) / n
!
!      var ( X(1:n) ) = sum ( ( X(1:n) - mean )^2 ) / ( n - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!    N should be at least 2.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance

  if ( n < 2 ) then

    variance = 0.0D+00

  else

    mean = sum ( a(1:n) ) / real ( n, kind = 8 )

    variance = sum ( ( a(1:n) - mean )**2 )

    variance = variance / real ( n - 1, kind = 8 )

  end if

  return
end
function radians_to_degrees ( radians )

!*****************************************************************************80
!
!! RADIANS_TO_DEGREES converts an angle measure from radians to degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RADIANS, the angle measure in radians.
!
!    Output, real ( kind = 8 ) RADIANS_TO_DEGREES, the angle measure in degrees.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) radians
  real ( kind = 8 ) radians_to_degrees

  radians_to_degrees = ( radians / pi ) * 180.0D+00

  return
end
subroutine random_initialize ( seed_input )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!
!    The G95 compiler is currently generating internal compiler errors
!    when it tries to compile this routine.  Just one more exasperation
!    on the mountain of complications because of the ragged interface with
!    the nonstandard random number generator standard!
!
!    If you don't initialize the FORTRAN90 random number generator
!    routine RANDOM_NUMBER, which is used by calls like
!
!      call random_number ( harvest = x )
!
!    then its behavior is not specified.  That may be OK for you.  But
!    you may want to be able to force the same sequence to be generated
!    each time, or to force a different sequence.  
!
!    To control the sequence of random numbers, you need to set the seed.
!    In FORTRAN90, this is done by calling the RANDOM+SEED routine.
!    You can call it with no arguments, in fact.  But if you call
!    it with no arguments:
!
!      call random_seed ( )
!
!    then its behavior (or more particularly, the behavior of RANDOM_NUMBER)
!    is still not specified.  You might hope that the system will go out
!    and pick a nice random seed for you, but there's no guarantee. 
!
!
!    For example, on the DEC ALPHA, if you compile a program that calls
!    RANDOM_NUMBER, then every time you run it, you get the same sequence
!    of "random" values.  If you compile a program that calls RANDOM_SEED
!    with no arguments, and then calls RANDOM_NUMBER, you still get the
!    same sequence each time you run the program.
!
!    In order to actually try to scramble up the random number generator 
!    a bit, this routine goes through the tedious process of getting the 
!    size of the random number seed, making up values based on the current 
!    time, and setting the random number seed.
!
!    Unfortunately, the RANDOM_SEED routine has a very elastic definition.
!    It does not use a single scalar integer SEED.  Instead, it communicates
!    with the user through an integer vector whose size is not specified.
!    You actually have to "ask" the routine to tell you the size of this
!    vector.  Then you can fill up the vector with values to be used to
!    influence the seeding of the random number routine.  The details of
!    how the seed affects the sequence are also unspecified, but the only
!    thing we are reasonably confident about is that specifying the same
!    seed should result in the same sequence, and specifying different
!    seeds should result in different sequences!
!
!    I assume this is far more than you wanted to know.  (It's certainly
!    more than I wanted to know!)
!
!    The argument SEED is an input quantity only, so it is legal
!    to type
!
!      call random_initialize ( 0 )
!
!    or
!
!      call random_initialize ( 18867 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED_INPUT, the user's "suggestion" for a seed.
!    However, if the input value is 0, the routine will come up with
!    its own "suggestion", based on the system clock.
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) count_max
  integer ( kind = 4 ) count_rate
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_input
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real ( kind = 4 ) t
  integer ( kind = 4 ), parameter :: warm_up = 100

  seed = seed_input
!
!  Initialize the random seed routine.
!
  call random_seed ( )
!
!  Determine the size of the random number seed vector.
!
  call random_seed ( size = seed_size )
!
!  Allocate a vector of the right size to be used as a random seed.
!
  allocate ( seed_vector(seed_size) )
!
!  If the user supplied a SEED value, use that.
!
!  Otherwise, use the system clock value to make up a value that is
!  likely to change based on when this routine is called.
!
  if ( seed /= 0 ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, user SEED = ', seed
    end if

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ', &
        seed
    end if

  end if
!
!  Set the seed vector.  We don't know the significance of the
!  individual entries of the internal seed vector, so we'll just set 
!  all entries to SEED.
!
  seed_vector(1:seed_size) = seed
!
!  Now call RANDOM_SEED, and tell it to use this seed vector.
!
  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times just to "warm it up".
!
  do i = 1, warm_up
    call random_number ( harvest = t )
  end do

  return
end
subroutine rat_factor ( m, n, maxfactor, factor_num, factor, power, &
  mleft, nleft )

!*****************************************************************************80
!
!! RAT_FACTOR factors a rational value into a product of prime factors.
!
!  Discussion:
!
!    ( M / N ) = ( MLEFT / NLEFT ) * Product ( 1 <= I <= FACTOR_NUM )
!      FACTOR(I)**POWER(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the top and bottom of a rational value.
!    The ratio of M and N must be positive.
!
!    Input, integer ( kind = 4 ) MAXFACTOR, the maximum number of factors for
!    which storage has been allocated.
!
!    Output, integer ( kind = 4 ) FACTOR_NUM, the number of prime factors 
!    of M/N.
!
!    Output, integer ( kind = 4 ) FACTOR(MAXFACTOR), the prime factors of M/N.
!
!    Output, integer ( kind = 4 ) POWER(MAXFACTOR).  POWER(I) is the power of
!    the FACTOR(I) in the representation of M/N.
!
!    Output, integer ( kind = 4 ) MLEFT, NLEFT, the top and bottom of 
!    the factor of M / N that remains.  If ABS ( MLEFT / NLEFT ) is not 1, then
!    the rational value was not completely factored.
!
  implicit none

  integer ( kind = 4 ) maxfactor

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mleft
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) p
  integer ( kind = 4 ) power(maxfactor)
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) prime_max

  factor_num = 0

  mleft = m
  nleft = n
!
!  NLEFT should be nonnegative.
!
  if ( nleft < 0 ) then
    mleft = -mleft
    nleft = -nleft
  end if

  if ( m == 0 .or. n == 0 ) then
    return
  end if

  if ( m == n ) then
    factor_num = 1
    factor(1) = 1
    power(1) = 1
    return
  end if
!
!  Find out how many primes we stored.
!
  prime_max = prime ( -1 )

  do i = 1, prime_max

    p = prime ( i )

    if ( mod ( nleft, p ) == 0 .or. &
         mod ( abs ( mleft ), p ) == 0 ) then

      if ( factor_num < maxfactor ) then

        factor_num = factor_num + 1
        factor(factor_num) = p
        power(factor_num) = 0
!
!  Divide MLEFT by PRIME(I) as often as you can.
!
        if ( mod ( abs ( mleft ), p ) == 0  ) then

          do

            power(factor_num) = power(factor_num) + 1
            mleft = mleft / p

            if ( mod ( abs ( mleft ), p ) /= 0 ) then
              exit
            end if

          end do

        end if
!
!  Divide NLEFT by PRIME(I) as often as you can.
!
        if ( mod ( nleft, p ) == 0  ) then

          do

            power(factor_num) = power(factor_num) - 1
            nleft = nleft / p

            if ( mod ( nleft, p ) /= 0 ) then
              exit
            end if

          end do

        end if

        if ( power(factor_num) == 0 ) then
          factor_num = factor_num - 1
        end if

      end if

    end if

  end do

  return
end
subroutine rickey ( ab, bb, er, f, h, hb, hp, r, so, tb, g )

!*****************************************************************************80
!
!! RICKEY evaluates Branch Rickey's baseball index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Schwarz,
!    Looking Beyond the Batting Average,
!    The New York Times, 
!    Sunday, 1 August 2004.
!    
!    Branch Rickey,
!    Goodby to Some Old Baseball Ideas,
!    Life Magazine,
!    2 August 1954.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) AB, number of at-bats.
!
!    Input, integer ( kind = 4 ) BB, base on balls.
!
!    Input, integer ( kind = 4 ) ER, earned runs.
!
!    Input, real ( kind = 8 ) F, the fielding rating.
!
!    Input, integer ( kind = 4 ) H, number of hits.
! 
!    Input, integer ( kind = 4 ) HB, hit batsmen.
!
!    Input, integer ( kind = 4 ) HP, hit by pitcher.
!
!    Input, integer ( kind = 4 ) R, runs.
!
!    Input, integer ( kind = 4 ) SO, strike outs.
!
!    Input, integer ( kind = 4 ) TB, total bases.
!
!    Output, real ( kind = 8 ) G, the Branch Rickey index, an estimate for the
!    expected winning percentage of a team with the given statistics.
!    (0.5 has already been subtracted from this value.)
!
  implicit none

  integer ( kind = 4 ) ab
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) er
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) h
  integer ( kind = 4 ) hb
  real ( kind = 8 ) hitting
  integer ( kind = 4 ) hp
  real ( kind = 8 ) pitching
  integer ( kind = 4 ) r
  integer ( kind = 4 ) so
  integer ( kind = 4 ) tb

  hitting = &
      real (  h + bb + hp,   kind = 8 ) / real ( ab + bb + hp, kind = 8 ) &
    + real ( 3 * ( tb - h ), kind = 8 ) / real ( 4 * ab,       kind = 8 ) &
    + real (              r, kind = 8 ) / real ( h + bb + hp,  kind = 8 )

  pitching = &
      real ( h,       kind = 8 ) / real ( ab,                   kind = 8 ) &
    + real ( bb + hb, kind = 8 ) / real ( ab + bb + hb,         kind = 8 ) &
    + real ( er,      kind = 8 ) / real ( h + bb + hb,          kind = 8 ) &
    - real ( so,      kind = 8 ) / real ( 8 * ( ab + bb + hb ), kind = 8 )

  g = hitting - pitching - f

  return
end
subroutine roots_to_i4poly ( n, x, c )

!*****************************************************************************80
!
!! ROOTS_TO_I4POLY converts polynomial roots to polynomial coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of roots specified.
!
!    Input, integer ( kind = 4 ) X(N), the roots.
!
!    Output, integer ( kind = 4 ) C(0:N), the coefficients of the polynomial.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x(n)
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
  c(0:n-1) = 0
  c(n) = 1
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
  do j = 1, n
    do i = 1, n + 1 - j
      c(n-i) = c(n-i) - x(n+1-i-j+1) * c(n-i+1)
    end do
  end do

  return
end
subroutine roots_to_r8poly ( n, x, c )

!*****************************************************************************80
!
!! ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of roots specified.
!
!    Input, real ( kind = 8 ) X(N), the roots.
!
!    Output, real ( kind = 8 ) C(0:N), the coefficients of the polynomial.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
  c(0:n-1) = 0.0D+00
  c(n) = 1.0D+00
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
  do j = 1, n
    do i = 1, n + 1 - j
      c(n-i) = c(n-i) - x(n+1-i-j+1) * c(n-i+1)
    end do
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
subroutine tuple_next2 ( n, xmin, xmax, x, rank )

!*****************************************************************************80
!
!! TUPLE_NEXT2 computes the next element of an integer tuple space.
!
!  Discussion:
!
!    The elements X are N vectors.
!
!    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
!
!    The elements are produced one at a time.
!
!    The first element is
!      (XMIN(1), XMIN(2), ..., XMIN(N)),
!    the second is (probably)
!      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
!    and the last element is
!      (XMAX(1), XMAX(2), ..., XMAX(N))
!
!    Intermediate elements are produced in a lexicographic order, with
!    the first index more important than the last, and the ordering of
!    values at a fixed index implicitly defined by the sign of
!    XMAX(I) - XMIN(I).
!
!  Example:
!
!    N = 2,
!    XMIN = (/ 1, 10 /)
!    XMAX = (/ 3,  8 /)
!
!    RANK    X
!    ----  -----
!      1   1 10
!      2   1  9
!      3   1  8
!      4   2 10
!      5   2  9
!      6   2  8
!      7   3 10
!      8   3  9
!      9   3  8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, integer ( kind = 4 ) XMIN(N), XMAX(N), the "minimum" and "maximum"
!    entry values.  These values are minimum and maximum only in the sense of 
!    the lexicographic ordering.  In fact, XMIN(I) may be less than, equal to,
!    or greater than XMAX(I).
!
!    Input/output, integer ( kind = 4 ) X(N), on input the previous tuple.
!    On output, the next tuple.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank of the item.  On first
!    call, set RANK to 0 to start up the sequence.  On return, if RANK is zero,
!    there are no more items in the sequence.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) xmin(n)
  integer ( kind = 4 ) xmax(n)

  if ( rank < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( product ( 1 + abs ( xmax(1:n) - xmin(1:n) ) ) < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( rank == 0 ) then
    x(1:n) = xmin(1:n)
    rank = 1
    return
  end if

  rank = rank + 1
  i = n

  do

    if ( x(i) /= xmax(i) ) then
      x(i) = x(i) + sign ( 1, xmax(i) - xmin(i) )
      exit
    end if

    x(i) = xmin(i)

    if ( i == 1 ) then
      rank = 0
      exit
    end if

    i = i - 1

  end do

  return
end
subroutine tvec_even ( nt, t )

!*****************************************************************************80
!
!! TVEC_EVEN computes evenly spaced angles between 0 and 2*PI.
!
!  Discussion:
!
!    The computation realizes that 0 = 2 * PI.
!
!  Example:
!
!    NT = 4
!
!    T = ( 0, PI/2, PI, 3*PI/2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles, in radians.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t(nt)

  do i = 1, nt
    t(i) = real ( 2 * ( i - 1 ), kind = 8 ) * pi &
         / real ( nt, kind = 8 )
  end do

  return
end
subroutine tvec_even2 ( nt, t )

!*****************************************************************************80
!
!! TVEC_EVEN2 computes evenly spaced angles between 0 and 2*PI.
!
!  Discussion:
!
!    The computation realizes that 0 = 2 * PI.  The values are equally
!    spaced in the circle, do not include 0, and are symmetric about 0.
!
!  Example:
!
!    NT = 4
!
!    T = ( PI/4, 3*PI/4, 5*PI/4, 7*PI/4 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles, in radians.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t(nt)

  do i = 1, nt
    t(i) = real ( 2 * i - 1, kind = 8 ) * pi &
         / real ( nt, kind = 8 )
  end do

  return
end
subroutine tvec_even3 ( nt, t )

!*****************************************************************************80
!
!! TVEC_EVEN3 computes evenly spaced angles between 0 and 2*PI.
!
!  Discussion:
!
!    The angles begin with 0 and end with 2*PI.
!
!  Example:
!
!    NT = 4
!
!    T = ( 0, 2*PI/3, 4*PI/3 2*PI )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles, in radians.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t(nt)

  if ( nt == 1 ) then
    t(1) = pi
  else
    do i = 1, nt
      t(i) = real ( 2 * ( i - 1 ), kind = 8 ) * pi &
           / real ( nt - 1, kind = 8 )
    end do
  end if

  return
end
subroutine tvec_even_bracket ( nt, theta1, theta2, t )

!*****************************************************************************80
!
!! TVEC_EVEN_BRACKET computes evenly spaced angles between THETA1 and THETA2.
!
!  Discussion:
!
!    The interval between THETA1 and THETA2 is divided into NT-1 subintervals.
!
!    The angles returned are the breakpoints of these subintervals,
!    including THETA1 and THETA2.
!
!  Example:
!
!    NT = 4
!    THETA1 = 30
!    THETA2 = 90
!
!    T = ( 30, 50, 70, 90 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the limiting angles.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  if ( nt == 1 ) then

    t(1) = ( theta1 + theta2 ) / 2.0D+00

  else

    do i = 1, nt
      t(i) = ( real ( nt - i,     kind = 8 ) * theta1   &
             + real (      i - 1, kind = 8 ) * theta2 ) &
             / real ( nt     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine tvec_even_bracket2 ( nt, theta1, theta2, t )

!*****************************************************************************80
!
!! TVEC_EVEN_BRACKET2 computes evenly spaced angles from THETA1 to THETA2.
!
!  Discussion:
!
!    The interval between THETA1 and THETA2 is divided into NT+1 subintervals.
!
!    The angles returned are the internal NT breakpoints of the subintervals.
!
!  Example:
!
!    NT = 5
!    THETA1 = 30
!    THETA2 = 90
!
!    T = ( 40, 50, 60, 70, 80 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the limiting angles.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  do i = 1, nt
    t(i) = ( real ( nt + 1 - i, kind = 8 ) * theta1   &
           + real (          i, kind = 8 ) * theta2 ) &
           / real ( nt + 1,     kind = 8 )
  end do

  return
end
subroutine tvec_even_bracket3 ( nt, theta1, theta2, t )

!*****************************************************************************80
!
!! TVEC_EVEN_BRACKET3 computes evenly spaced angles from THETA1 to THETA2.
!
!  Discussion:
!
!    The interval between THETA1 and THETA2 is divided into NT subintervals.
!
!    The angles returned are the midpoints of each subinterval.
!
!  Example:
!
!    NT = 3
!    THETA1 = 30
!    THETA2 = 90
!
!    T = ( 40, 60, 80 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the limiting angles.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  do i = 1, nt
    t(i) = ( real ( 2 * nt - 2 * i + 1, kind = 8 ) * theta1   &
           + real (          2 * i - 1, kind = 8 ) * theta2 ) &
           / real ( 2 * nt,             kind = 8 )
  end do

  return
end
subroutine upc_check_digit ( p, l, r, c )

!*****************************************************************************80
!
!! UPC_CHECK_DIGIT returns the check digit of a UPC.
!
!  Discussion:
!
!    UPC stands for Universal Price Code.
!
!    A full UPC is a string of 12 digits, in groups of size 1, 5, 5, and 1,
!    of the form P-LLLLL-RRRRR-C, where:
!
!      P is the one-digit product type code.
!      L is the five-digit manufacturer code.
!      R is the five_digit product code
!      C is the check digit.
!
!  Example:
!
!    0-72890-00011-8
!    0-12345-67890-5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the one-digit product type code.
!
!    Input, integer ( kind = 4 ) L, the five-digit manufacturer code.
!
!    Input, integer ( kind = 4 ) R, the five-digit product code.
!
!    Output, integer ( kind = 4 ) C, the check digit.
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc(5)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) r
  integer ( kind = 4 ) rc(5)

  if ( p < 0 .or. 9 < p ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UPC_CHECK_DIGIT - Fatal error!'
    write ( *, '(a)' ) '  P < 0 or 9 < P!'
    stop
  end if

  if ( l < 0 .or. 99999 < l ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UPC_CHECK_DIGIT - Fatal error!'
    write ( *, '(a)' ) '  L < 0 or 99999 < L!'
    stop
  end if

  if ( r < 0 .or. 99999 < r ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UPC_CHECK_DIGIT - Fatal error!'
    write ( *, '(a)' ) '  R < 0 or 99999 < R!'
    stop
  end if

  call i4_to_digits_decimal ( l, 5, lc )
  call i4_to_digits_decimal ( r, 5, rc )

  c = ( p + lc(2) + lc(4) + rc(1) + rc(3) + rc(5) ) * 3 &
          + lc(1) + lc(3) + lc(5) + rc(2) + rc(4)

  c = mod ( c, 10 )

  c = mod ( 10 - c, 10 )

  return
end
function versine_pulse ( t, ta, tb, v1, amp )

!*****************************************************************************80
!
!! VERSINE_PULSE adds a versine pulse to a constant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) TA, the time at which the pulse begins.
!
!    Input, real ( kind = 8 ) TB, the time at which the pulse finishes.
!
!    Input, real ( kind = 8 ) V1, the constant value.
!
!    Input, real ( kind = 8 ) AMP, the amplitude of the pulse.
!
!    Output, real ( kind = 8 ) VERSINE_PULSE, the value of the signal at time T.
!
  implicit none

  real ( kind = 8 ) amp
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) ta
  real ( kind = 8 ) tb
  real ( kind = 8 ) v
  real ( kind = 8 ) v1
  real ( kind = 8 ) versine_pulse

  v = v1

  if ( ta <= t .and. t <= tb ) then
        v = v + ( 0.5D+00 * amp &
          * ( 1.0D+00 - cos ( 2.0D+00 * pi * ( t - ta ) / ( tb - ta ) ) ) )
  end if

  versine_pulse = v

  return
end
