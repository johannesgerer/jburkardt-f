program main

!*****************************************************************************80
!
!! MAIN is the main program for LEGENDRE_RULE_FAST.
!
!  Discussion:
!
!    This program computes a standard Gauss-Legendre quadrature rule
!    and writes it to a file.
!
!  Usage:
!
!    legendre_rule_fast ( n, a, b )
!
!    where
!
!    * n is the number of points in the rule;
!    * a is the left endpoint;
!    * b is the right endpoint.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  integer ( kind = 4 ) arg_num
  real ( kind = 8 ) b
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  integer ( kind = 4 ) n
  character ( len = 255 ) string

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_RULE_FAST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute a Gauss-Legendre rule for approximating'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Integral ( b <= x <= b ) f(x) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  of order N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The computed rule is written to 3 files:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    leg_oN_w.txt - the weight file'
  write ( *, '(a)' ) '    leg_oN_x.txt - the abscissa file.'
  write ( *, '(a)' ) '    leg_oN_r.txt - the region file.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get N.
!
  if ( 1 <= arg_num ) then
  
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, n, ierror, last )
    
  else
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the rule order N:'
    read ( *, * ) n
    
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The requested order N ', n
!
!  Get A.
!
  if ( 2 <= arg_num ) then
  
    iarg = 2
    call getarg ( iarg, string )
    call s_to_r8 ( string, a, ierror, last )
    
  else
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the left endpoint A:'
    read ( *, * ) a
    
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The left endpoint A is ', a
!
!  Get B.
!
  if ( 3 <= arg_num ) then
  
    iarg = 3
    call getarg ( iarg, string )
    call s_to_r8 ( string, b, ierror, last )
    
  else
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the right endpoint B:'
    read ( *, * ) b
    
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The right endpoint B is ', b
!
!  Construct the rule and output it.
!
  call legendre_handle ( n, a, b )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_RULE_FAST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
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
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.  
!
!  Discussion:
!
!    CH_EQI ( 'A', 'a' ) is TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!     CH  DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal. 
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
!    If CH was 'illegal', then DIGIT is -1.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then
 
    digit = iachar ( ch ) - 48
 
  else if ( ch == ' ' ) then
 
    digit = 0
 
  else

    digit = -1

  end if
 
  return
end
subroutine digit_to_ch ( digit, ch )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Discussion:
!
!    Instead of CHAR, we now use the ACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!    DIGIT   CH 
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character CH, the corresponding character.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    ch = achar ( digit + 48 )

  else

    ch = '*'

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
!    15 January 2008
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
subroutine i4_to_s_left ( i4, s )

!*****************************************************************************80
!
!! I4_TO_S_LEFT converts an I4 to a left-justified string.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ).
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!        I4  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!   1234567  ******  <-- Not enough room!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be left-justified.  If there is not enough space,
!    the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * )  s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = i4
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = -ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit of IVAL, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

    if ( ival == 0 ) then
      exit
    end if

  end do
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '
 
  return
end
subroutine legendre_compute_glr ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_GLR: Legendre quadrature by the Glaser-Liu-Rokhlin method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 October 2009
!
!  Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
!    A fast algorithm for the calculation of the roots of special functions, 
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) p
  real ( kind = 8 ) pp
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
!
!  Get the value and derivative of the N-th Legendre polynomial at 0.0.
!
  call legendre_compute_glr0 ( n, p, pp )
!
!  If N is odd, then zero is a root.
!
  if ( mod ( n, 2 ) == 1 ) then

    x((n+1)/2) = 0.0D+00
    w((n+1)/2) = pp
!
!  If N is even, we have to compute a root.
!
  else

    call legendre_compute_glr2 ( p, n, x((n/2)+1), w((n/2)+1) )

  end if
!
!  Get the complete set of roots and derivatives.
!
  call legendre_compute_glr1 ( n, x, w )
!
!  Compute W.
!
  w(1:n) = 2.0D+00 / &
    ( 1.0D+00 - x(1:n) ) / ( 1.0D+00 + x(1:n) ) / w(1:n) / w(1:n)

  w(1:n) = 2.0D+00 * w(1:n) / sum ( w(1:n) )

  return
end
subroutine legendre_compute_glr0 ( n, p, pp )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_GLR0 gets a starting value for the fast algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2009
!
!  Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
!    A fast algorithm for the calculation of the roots of special functions, 
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Legendre polynomial.
!
!    Output, real ( kind = 8 ) P, PP, the value of the N-th Legendre polynomial
!    and its derivative at 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) k
  real ( kind = 8 ) p
  real ( kind = 8 ) pm1
  real ( kind = 8 ) pm2
  real ( kind = 8 ) pp
  real ( kind = 8 ) ppm1
  real ( kind = 8 ) ppm2
  real ( kind = 8 ) rk
!
!  Compute coefficients of P_m(0), Pm'(0), m = 0,..,N
!
  pm2 = 0.0D+00
  pm1 = 1.0D+00
  ppm2 = 0.0D+00
  ppm1 = 0.0D+00

  do k = 0, n - 1
    rk = real ( k, kind = 8 )
    p = - rk * pm2 / ( rk + 1.0D+00 )
    pp = ( ( 2.0D+00 * rk + 1.0D+00 ) * pm1 &
                     - rk             * ppm2 ) &
         / (           rk + 1.0D+00 )
    pm2 = pm1
    pm1 = p
    ppm2 = ppm1
    ppm1 = pp
  end do

  return
end
subroutine legendre_compute_glr1 ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_GLR1 gets the complete set of Legendre points and weights.
!
!  Discussion:
!
!    This routine requires that a starting estimate be provided for one
!    root and its derivative.  This information will be stored in entry
!    (N+1)/2 if N is odd, or N/2 if N is even, of X and W.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2009
!
!  Author:
!
!    Original C++ version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
!    A fast algorithm for the calculation of the roots of special functions, 
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Legendre polynomial.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, a starting value
!    has been set in one entry.  On output, the roots of the Legendre 
!    polynomial.
!
!    Input/output, real ( kind = 8 ) W(N).  On input, a starting value
!    has been set in one entry.  On output, the derivatives of the Legendre 
!    polynomial at the zeros.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) M, the number of terms in the Taylor expansion.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 30
  integer ( kind = 4 ) n

  real ( kind = 8 ) dk
  real ( kind = 8 ) dn
  real ( kind = 8 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n2
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) rk2_leg
  integer ( kind = 4 ) s
  real ( kind = 8 ) ts_mult
  real ( kind = 8 ) u(m+2)
  real ( kind = 8 ) up(m+1)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xp

  if ( mod ( n, 2 ) == 1 ) then
    n2 = ( n - 1 ) / 2 - 1
    s = 1
  else
    n2 = n / 2 - 1
    s = 0
  end if

  dn = real ( n, kind = 8 )

  do j = n2 + 1, n - 2

    xp = x(j+1)

    h = rk2_leg ( pi/2.0D+00, -pi/2.0D+00, xp, n ) - xp

    u(1) = 0.0D+00
    u(2) = 0.0D+00
    u(3) = w(j+1)

    up(1) = 0.0D+00
    up(2) = u(3)

    do k = 0, m - 2

      dk = real ( k, kind = 8 )

      u(k+4) = &
      ( &
        2.0D+00 * xp * ( dk + 1.0D+00 ) * u(k+3) &
        + ( dk * ( dk + 1.0D+00 ) - dn * ( dn + 1.0D+00 ) ) * u(k+2) / ( dk + 1.0D+00 ) &
      ) / ( 1.0D+00 - xp ) / ( 1.0D+00 + xp ) / ( dk + 2.0D+00 )

      up(k+3) = ( dk + 2.0D+00 ) * u(k+4)

    end do

    do l = 0, 4
      h = h - ts_mult ( u, h, m ) / ts_mult ( up, h, m-1 )
    end do

    x(j+2) = xp + h
    w(j+2) = ts_mult ( up, h, m - 1 )   

  end do

  do k = 0, n2 + s
    x(k+1) = - x(n-1-k+1)
    w(k+1) = w(n-1-k+1)
  end do

  return
end
subroutine legendre_compute_glr2 ( pn0, n, x1, d1 )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_GLR2 finds the first real root.
!
!  Discussion:
!
!    This routine is only called if N is even.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 October 2009
!
!  Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
!    A fast algorithm for the calculation of the roots of special functions, 
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PN0, the value of the N-th Legendre polynomial
!    at 0.
!
!    Input, integer ( kind = 4 ) N, the order of the Legendre polynomial.
!
!    Output, real ( kind = 8 ) X1, the first real root.
!
!    Output, real ( kind = 8 ) D1, the derivative at X1.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) M, the number of terms in the Taylor expansion.
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 30

  real ( kind = 8 ) d1
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) pn0
  real ( kind = 8 ) rk
  real ( kind = 8 ) rn
  real ( kind = 8 ) scale
  real ( kind = 8 ) step
  real ( kind = 8 ) theta
  real ( kind = 8 ) u(m+1)
  real ( kind = 8 ) up(m+1)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x1k(m+1)

  k = ( n + 1 ) / 2

  theta = pi * real ( 4 * k - 1, kind = 8 ) / real ( 4 * n + 2, kind = 8 )

  x1 = ( 1.0D+00 - real ( n - 1, kind = 8 ) &
    / real ( 8 * n * n * n, kind = 8 ) &
    - 1.0D+00 / real ( 384 * n * n * n * n, kind = 4 ) &
    * ( 39.0D+00 - 28.0D+00 / ( sin ( theta ) * sin ( theta ) ) ) ) &
    * cos ( theta )
!
!  Scaling.
!
  scale = 1.0D+00 / x1
!
!  Recurrence relation for Legendre polynomials.
!
  u(1:m+1) = 0.0D+00
  up(1:m+1) = 0.0D+00

  rn = real ( n, kind = 8 )

  u(1) = pn0

  do k = 0, m - 2, 2

    rk = real ( k, kind = 8 )

    u(k+3) = ( rk * ( rk + 1.0D+00 ) - rn * ( rn + 1.0D+00 ) ) * u(k+1) &
      / ( rk + 1.0D+00 ) / ( rk + 2.0D+00 ) / scale / scale

    up(k+2) = ( rk + 2.0D+00 ) * u(k+3) * scale

  end do
!
!  Flip for more accuracy in inner product calculation
!
  u = u(m+1:1:-1)
  up = up(m+1:1:-1)

  x1k(1:m+1) = 1.0D+00

  step = huge ( step )
  l = 0
!
!  Newton iteration.
!
  eps = epsilon ( eps )

  do while ( eps < abs ( step ) .and. l < 10 )
    l = l + 1
    step = dot_product ( u(1:m+1),  x1k(1:m+1) ) &
         / dot_product ( up(1:m+1), x1k(1:m+1) )
    x1 = x1 - step
    x1k(1) = 1.0D+00
    x1k(2) = scale * x1
    do kk = 3, m + 1
      x1k(kk) = x1k(kk-1) * scale * x1
    end do
    x1k(1:m+1) = x1k(m+1:1:-1)
  end do

  d1 = dot_product ( up(1:m+1), x1k(1:m+1) )

  return
end
subroutine legendre_handle ( n, a, b )

!*****************************************************************************80
!
!! LEGENDRE_HANDLE computes the requested Gauss-Legendre rule and outputs it.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Input, real ( kind = 8 ) A, B, the left and right endpoints.
! 
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  character ( len = 80 ) output_r
  character ( len = 80 ) output_w
  character ( len = 80 ) output_x
  real ( kind = 8 ) r(2)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  character ( len = 10 ) tag
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  r(1) = a
  r(2) = b
!
!  Compute the rule.
!
  allocate ( w(n) )
  allocate ( x(n) )

  call cpu_time ( t1 )
  call legendre_compute_glr ( n, x, w )
  call cpu_time ( t2 )
  
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a)' ) '  Computation required ', t2 - t1, ' seconds.'
!
!  Rescale the data.
!
  call rescale ( n, a, b, x, w )
!
!  Write the data to files.
!
  call i4_to_s_left ( n, tag )

  output_w = 'leg_o' // trim ( tag ) // '_w.txt'
  output_x = 'leg_o' // trim ( tag ) // '_x.txt'
  output_r = 'leg_o' // trim ( tag ) // '_r.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weight file will be   "' // trim ( output_w ) // '".'
  write ( *, '(a)' ) '  Abscissa file will be "' // trim ( output_x ) // '".'
  write ( *, '(a)' ) '  Region file will be   "' // trim ( output_r ) // '".'
            
  call r8mat_write ( output_w, 1, n, w )
  call r8mat_write ( output_x, 1, n, x )
  call r8mat_write ( output_r, 1, 2, r )
      
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * )  output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine rescale ( n, a, b, x, w )

!*****************************************************************************80
!
!! RESCALE rescales a Legendre quadrature rule from [-1,+1] to [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original MATLAB version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
!    A fast algorithm for the calculation of the roots of special functions, 
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the new interval.
!
!    Input/output, real ( kind = 8 ) X(N), on input, the abscissas for [-1,+1].
!    On output, the abscissas for [A,B].
!
!    Input/output, real ( kind = 8 ) W(N), on input, the weights for [-1,+1].
!    On output, the weights for [A,B].
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  x(1:n) = ( ( a + b ) + ( b - a ) * x(1:n) ) / 2.0D+00
  w(1:n) = ( b - a ) * w(1:n) / 2.0D+00

  return
end
function rk2_leg ( t1, t2, x, n )

!*****************************************************************************80
!
!! RK2_LEG advances the value of X(T) using a Runge-Kutta method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2009
!
!  Author:
!
!    Original C++ version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T1, T2, the range of the integration interval.
!
!    Input, real ( kind = 8 ) X, the value of X at T1.
!
!    Input, integer ( kind = 4 ) N, the number of steps to take.
!
!    Output, real ( kind = 8 ) RK2_LEG, the value of X at T2.
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) h
  integer ( kind = 4 ) j
  real ( kind = 8 ) k1
  real ( kind = 8 ) k2
  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ) n
  real ( kind = 8 ) rk2_leg
  real ( kind = 8 ) snn1
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  x2 = x

  h = ( t2 - t1 ) / real ( m, kind = 8 )
  snn1 = sqrt ( real ( n * ( n + 1 ), kind = 8 ) )
  t = t1

  do j = 0, m - 1

    f = ( 1.0D+00 - x2 ) * ( 1.0D+00 + x2 )
    k1 = - h * f / ( snn1 * sqrt ( f ) - 0.5D+00 * x2 * sin ( 2.0D+00 * t ) )
    x2 = x2 + k1

    t = t + h

    f = ( 1.0D+00 - x2 ) * ( 1.0D+00 + x2 )
    k2 = - h * f / ( snn1 * sqrt ( f ) - 0.5D+00 * x2 * sin ( 2.0D+00 * t ) )
    x2 = x2 + 0.5D+00 * ( k2 - k1 )

  end do

  rk2_leg = x2

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S 
!    used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) length
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = 8 )".
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  integer ( kind = 4 ) s_length
  character :: TAB = achar ( 9 )

  s_length = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( s_length < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length+1 == s_length ) then
    length = s_length
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

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

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
function ts_mult ( u, h, n )

!*****************************************************************************80
!
!! TS_MULT...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2009
!
!  Author:
!
!    Original C++ version by Nick Hale.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U(N+1), ...
!
!    Input, real ( kind = 8 ) H, ...
!
!    Input, integer ( kind = 4 ) N, ...
!
!    Output, real ( kind = 8 ) TS_MULT, ...
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) h
  real ( kind = 8 ) hk
  integer ( kind = 4 ) k
  real ( kind = 8 ) ts
  real ( kind = 8 ) ts_mult
  real ( kind = 8 ) u(n+1)
  
  ts = 0.0D+00
  hk = 1.0D+00
  do k = 1, n
    ts = ts + u(k+1) * hk
    hk = hk * h
  end do

  ts_mult = ts

  return
end
