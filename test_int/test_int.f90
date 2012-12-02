function besj0 ( x )

!*****************************************************************************80
!
!! BESJ0 computes J Bessel functions of order 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    William Cody, Argonne National Laboratory,
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argment of the Bessel function.
!
!    Output, real ( kind = 8 ) BESJ0, the value of the Bessel function.
!
  implicit none

  real ( kind = 8 ) besj0
  integer ( kind = 4 ) jint
  real ( kind = 8 ) result
  real ( kind = 8 ) x

  jint = 0
  call caljy0 ( x, result, jint )
  besj0 = result

  return
end
function besy0 ( x )

!*****************************************************************************80
!
!! BESY0 computes Y Bessel functions of order 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    William Cody, Argonne National Laboratory,
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argment of the Bessel function.
!
!    Output, real ( kind = 8 ) BESY0, the value of the Bessel function.
!
  implicit none

  real ( kind = 8 ) besy0
  integer ( kind = 4 ) jint
  real ( kind = 8 ) result
  real ( kind = 8 ) x

  jint = 1
  call caljy0 ( x, result, jint )
  besy0 = result

  return
end
subroutine caljy0 ( arg, result, jint )

!*****************************************************************************80
!
!! CALJY0 calculates J and Y Bessel functions of zero order.
!
!  Discussion:
!
!    This routine computes zero-order Bessel functions of the first and
!    second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
!    for Y0, and |X| <= XMAX for J0.
!
!    The routine CALJY0 is intended for internal use only, all
!    computations within the packet being concentrated in
!    this one routine.  The function subprograms invoke CALJY0 with
!    the statement
!      CALL CALJY0(ARG,RESULT,JINT),
!    where the parameter usage is as follows:
!
!      Function                  Parameters for CALJY0
!       call              ARG             RESULT          JINT
!
!     BESJ0(ARG)     |ARG| .LE. XMAX       J0(ARG)          0
!     BESY0(ARG)   0 .LT. ARG .LE. XMAX    Y0(ARG)          1
!
!    The main computation uses unpublished minimax rational
!    approximations for X .LE. 8.0, and an approximation from Hart,
!    et. al., for arguments larger than 8.0   Part of this
!    transportable packet is patterned after the machine dependent
!    FUNPACK program BESJ0(X), but cannot match that version for
!    efficiency or accuracy.  This version uses rational functions
!    that are theoretically accurate to at least 18 significant decimal
!    digits for X <= 8, and at least 18 decimal places for 8 < X.  The
!    accuracy achieved depends on the arithmetic system, the compiler,
!    the intrinsic functions, and proper selection of the machine
!    dependent constants.
!
!  Machine dependent constants:
!
!    XMAX   = largest acceptable argument.  The functions AINT, SIN
!             and COS must perform properly for  ABS(X) .LE. XMAX.
!             We recommend that XMAX be a small integer multiple of
!             sqrt(1/eps), where eps is the smallest positive number
!             such that 1 < 1+eps.
!
!    XSMALL = positive argument such that  1.0-(X/2)^2 = 1.0
!             to machine precision for all  ABS(X) .LE. XSMALL.
!             We recommend that  XSMALL < sqrt(eps)/beta, where beta
!             is the floating-point radix (usually 2 or 16).
!
!    Approximate values for some important machines are
!
!                            eps      XMAX     XSMALL
!
!    CDC 7600      (S.P.)  7.11D-15  1.34D+08  2.98D-08
!    CRAY-1        (S.P.)  7.11D-15  1.34D+08  2.98D-08
!    IBM PC (8087) (S.P.)  5.96D-08  8.19D+03  1.22D-04
!    IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09
!    IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13
!    UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10
!    VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Argonne National Laboratory.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
!    John Rice, Henry Thatcher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument of the Bessel function.
!
!    Output, real ( kind = 8 ) RESULT, the value of the Bessel function.
!
!    Input, integer ( kind = 4 ) JINT, is 0 for the J0 Bessel function, and
!    1 for the Y0 Bessel function.
!
  implicit none

  real ( kind = 8 ) arg
  real ( kind = 8 ) ax
  real ( kind = 8 ) cons
  real ( kind = 8 ) down
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jint
  real ( kind = 8 ) pi2
  real ( kind = 8 ) pj0(7)
  real ( kind = 8 ) pj1(8)
  real ( kind = 8 ) plg(4)
  real ( kind = 8 ) prod
  real ( kind = 8 ) py0(6)
  real ( kind = 8 ) py1(7)
  real ( kind = 8 ) py2(8)
  real ( kind = 8 ) p0(6)
  real ( kind = 8 ) p1(6)
  real ( kind = 8 ) p17
  real ( kind = 8 ) qj0(5)
  real ( kind = 8 ) qj1(7)
  real ( kind = 8 ) qlg(4)
  real ( kind = 8 ) qy0(5)
  real ( kind = 8 ) qy1(6)
  real ( kind = 8 ) qy2(7)
  real ( kind = 8 ) q0(5)
  real ( kind = 8 ) q1(5)
  real ( kind = 8 ) resj
  real ( kind = 8 ) result
  real ( kind = 8 ) r0
  real ( kind = 8 ) r1
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) twopi
  real ( kind = 8 ) twopi1
  real ( kind = 8 ) twopi2
  real ( kind = 8 ) two56
  real ( kind = 8 ) up
  real ( kind = 8 ) w
  real ( kind = 8 ) wsq
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xmax = 8.19D+03
  real ( kind = 8 ) xnum
  real ( kind = 8 ), parameter :: xsmall = 1.22D-09
  real ( kind = 8 ) xj0
  real ( kind = 8 ) xj1
  real ( kind = 8 ) xj01
  real ( kind = 8 ) xj02
  real ( kind = 8 ) xj11
  real ( kind = 8 ) xj12
  real ( kind = 8 ) xy
  real ( kind = 8 ) xy0
  real ( kind = 8 ) xy01
  real ( kind = 8 ) xy02
  real ( kind = 8 ) xy1
  real ( kind = 8 ) xy11
  real ( kind = 8 ) xy12
  real ( kind = 8 ) xy2
  real ( kind = 8 ) xy21
  real ( kind = 8 ) xy22
  real ( kind = 8 ) z
  real ( kind = 8 ) zsq
!-------------------------------------------------------------------
!  Mathematical constants
!    CONS = ln(.5) + Euler's gamma
!-------------------------------------------------------------------
  data p17/1.716D-1/
  data two56,cons/256.0D+0,-1.1593151565841244881D-1/
  data pi2,twopi/6.3661977236758134308D-1,6.2831853071795864769D+0/
  data twopi1,twopi2/6.28125D+0,1.9353071795864769253D-3/
!-------------------------------------------------------------------
!  Zeroes of Bessel functions
!-------------------------------------------------------------------
  data xj0/2.4048255576957727686D+0/,xj1/5.5200781102863106496D+0/
  data xy0/8.9357696627916752158D-1/,xy1/3.9576784193148578684D+0/
  data xy2/7.0860510603017726976D+0/
  data xj01/ 616.0D+0/, xj02/-1.4244423042272313784D-03/
  data xj11/1413.0D+0/, xj12/ 5.4686028631064959660D-04/
  data xy01/ 228.0D+0/, xy02/ 2.9519662791675215849D-03/
  data xy11/1013.0D+0/, xy12/ 6.4716931485786837568D-04/
  data xy21/1814.0D+0/, xy22/ 1.1356030177269762362D-04/
!-------------------------------------------------------------------
!  Coefficients for rational approximation to ln(x/a)
!--------------------------------------------------------------------
  data plg/-2.4562334077563243311D+01,2.3642701335621505212D+02, &
           -5.4989956895857911039D+02,3.5687548468071500413D+02/
  data qlg/-3.5553900764052419184D+01,1.9400230218539473193D+02, &
           -3.3442903192607538956D+02,1.7843774234035750207D+02/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J0(X) / (X^2 - XJ0^2),  XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
  data pj0/6.6302997904833794242D+06,-6.2140700423540120665D+08, &
          2.7282507878605942706D+10,-4.1298668500990866786D+11, &
         -1.2117036164593528341D-01, 1.0344222815443188943D+02, &
         -3.6629814655107086448D+04/

  data qj0/4.5612696224219938200D+05, 1.3985097372263433271D+08, &
         2.6328198300859648632D+10, 2.3883787996332290397D+12, &
          9.3614022392337710626D+02/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J0(X) / (X^2 - XJ1^2),  4.0  <  |X|  <=  8.0
!-------------------------------------------------------------------
  data pj1/4.4176707025325087628D+03, 1.1725046279757103576D+04, &
        1.0341910641583726701D+04,-7.2879702464464618998D+03, &
        -1.2254078161378989535D+04,-1.8319397969392084011D+03, &
         4.8591703355916499363D+01, 7.4321196680624245801D+02/

  data qj1/3.3307310774649071172D+02,-2.9458766545509337327D+03, &
          1.8680990008359188352D+04,-8.4055062591169562211D+04, &
          2.4599102262586308984D+05,-3.5783478026152301072D+05, &
         -2.5258076240801555057D+01/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY0) J0(X)) / (X^2 - XY0^2),
!        XSMALL  <  |X|  <=  3.0
!--------------------------------------------------------------------
  data py0/1.0102532948020907590D+04,-2.1287548474401797963D+06, &
         2.0422274357376619816D+08,-8.3716255451260504098D+09, &
         1.0723538782003176831D+11,-1.8402381979244993524D+01/

  data qy0/6.6475986689240190091D+02, 2.3889393209447253406D+05, &
         5.5662956624278251596D+07, 8.1617187777290363573D+09, &
         5.8873865738997033405D+11/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY1) J0(X)) / (X^2 - XY1^2),
!        3.0  <  |X|  <=  5.5
!--------------------------------------------------------------------
  data py1/-1.4566865832663635920D+04, 4.6905288611678631510D+06, &
         -6.9590439394619619534D+08, 4.3600098638603061642D+10, &
         -5.5107435206722644429D+11,-2.2213976967566192242D+13, &
          1.7427031242901594547D+01/

  data qy1/ 8.3030857612070288823D+02, 4.0669982352539552018D+05, &
          1.3960202770986831075D+08, 3.4015103849971240096D+10, &
          5.4266824419412347550D+12, 4.3386146580707264428D+14/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY2) J0(X)) / (X^2 - XY2^2),
!        5.5  <  |X|  <=  8.0
!--------------------------------------------------------------------
  data py2/ 2.1363534169313901632D+04,-1.0085539923498211426D+07, &
          2.1958827170518100757D+09,-1.9363051266772083678D+11, &
         -1.2829912364088687306D+11, 6.7016641869173237784D+14, &
         -8.0728726905150210443D+15,-1.7439661319197499338D+01/

  data qy2/ 8.7903362168128450017D+02, 5.3924739209768057030D+05, &
          2.4727219475672302327D+08, 8.6926121104209825246D+10, &
          2.2598377924042897629D+13, 3.9272425569640309819D+15, &
          3.4563724628846457519D+17/
!-------------------------------------------------------------------
!  Coefficients for Hart's approximation,  8.0 < |X|.
!-------------------------------------------------------------------
  data p0/3.4806486443249270347D+03, 2.1170523380864944322D+04, &
        4.1345386639580765797D+04, 2.2779090197304684302D+04, &
        8.8961548424210455236D-01, 1.5376201909008354296D+02/

  data q0/3.5028735138235608207D+03, 2.1215350561880115730D+04, &
        4.1370412495510416640D+04, 2.2779090197304684318D+04, &
        1.5711159858080893649D+02/

  data p1/-2.2300261666214198472D+01,-1.1183429920482737611D+02, &
        -1.8591953644342993800D+02,-8.9226600200800094098D+01, &
        -8.8033303048680751817D-03,-1.2441026745835638459D+00/

  data q1/1.4887231232283756582D+03, 7.2642780169211018836D+03, &
        1.1951131543434613647D+04, 5.7105024128512061905D+03, &
        9.0593769594993125859D+01/
!
!  Check for error conditions.
!
  ax = abs ( arg )

  if ( jint == 1 .and. arg <= 0.0D+00 ) then
    result = - r8_huge ( )
    return
  end if

  if ( xmax < ax ) then
    result = 0.0D+00
    return
  end if

  if ( ax <= xsmall ) then
    if ( jint == 0 ) then
      result = 1.0D+00
    else
      result = pi2 * ( log ( ax ) + cons )
    end if
    return
  end if
!
!  Calculate J0 for appropriate interval, preserving
!  accuracy near the zero of J0.
!
  if ( ax <= 8.0D+00 ) then

    zsq = ax * ax

    if ( ax <= 4.0D+00 ) then

      xnum = ( pj0(5) * zsq + pj0(6) ) * zsq + pj0(7)
      xden = zsq + qj0(5)

      do i = 1, 4
        xnum = xnum * zsq + pj0(i)
        xden = xden * zsq + qj0(i)
      end do

      prod = ( ( ax - xj01 / two56 ) - xj02 ) * ( ax + xj0 )

    else

      wsq = 1.0D+00 - zsq / 64.0D+00
      xnum = pj1(7) * wsq + pj1(8)
      xden = wsq + qj1(7)

      do i = 1, 6
        xnum = xnum * wsq + pj1(i)
        xden = xden * wsq + qj1(i)
      end do

      prod = ( ax + xj1 ) * ( ( ax - xj11 / two56 ) - xj12 )

    end if

    result = prod * xnum / xden

    if ( jint == 0 ) then
      return
    end if
!
!  Calculate Y0.
!  First find  RESJ = pi/2 ln(x/xn) J0(x), where XN is a zero of Y0.
!
    if ( ax <= 3.0D+00 ) then
      up = ( ax - xy01 / two56 ) - xy02
      xy = xy0
    else if ( ax <= 5.5D+00 ) then
      up = ( ax - xy11 / two56 ) - xy12
      xy = xy1
    else
      up = ( ax - xy21 / two56 ) - xy22
      xy = xy2
    end if

    down = ax + xy

    if ( abs ( up ) < p17 * down ) then

      w = up / down
      wsq = w * w
      xnum = plg(1)
      xden = wsq + qlg(1)

      do i = 2, 4
        xnum = xnum * wsq + plg(i)
        xden = xden * wsq + qlg(i)
      end do

      resj = pi2 * result * w * xnum / xden

    else

      resj = pi2 * result * log ( ax / xy )

    end if
!
!  Now calculate Y0 for appropriate interval, preserving
!  accuracy near the zero of Y0.
!
    if ( ax <= 3.0D+00 ) then

      xnum = py0(6) * zsq + py0(1)
      xden = zsq + qy0(1)

      do i = 2, 5
        xnum = xnum * zsq + py0(i)
        xden = xden * zsq + qy0(i)
      end do

    else if (ax <= 5.5D+00 ) then

      xnum = py1(7) * zsq + py1(1)
      xden = zsq + qy1(1)

      do i = 2, 6
        xnum = xnum * zsq + py1(i)
        xden = xden * zsq + qy1(i)
      end do

    else

      xnum = py2(8) * zsq + py2(1)
      xden = zsq + qy2(1)

      do i = 2, 7
        xnum = xnum * zsq + py2(i)
        xden = xden * zsq + qy2(i)
      end do

    end if

    result = resj + up * down * xnum / xden
!
!  Calculate J0 or Y0 for 8.0 < |ARG|.
!
  else

    z = 8.0D+00 / ax
    w = ax / twopi
    w = aint ( w ) + 0.125D+00
    w = ( ax - w * twopi1 ) - w * twopi2
    zsq = z * z
    xnum = p0(5) * zsq + p0(6)
    xden = zsq + q0(5)
    up = p1(5) * zsq + p1(6)
    down = zsq + q1(5)

    do i = 1, 4
      xnum = xnum * zsq + p0(i)
      xden = xden * zsq + q0(i)
      up = up * zsq + p1(i)
      down = down * zsq + q1(i)
    end do

    r0 = xnum / xden
    r1 = up / down

    if ( jint == 0 ) then
      result = sqrt ( pi2 / ax ) * ( r0 * cos ( w ) - z * r1 * sin ( w ) )
    else
      result = sqrt ( pi2 / ax ) * ( r0 * sin ( w ) + z * r1 * cos ( w ) )
    end if

  end if

  return
end
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
subroutine csevl ( x, cs, n, value )

!*****************************************************************************80
!
!! CSEVL evaluates a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 2000
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    R. Broucke,
!    Algorithm 446,
!    Communications of the ACM,
!    Volume 16, 254 (1973).
!
!    Leslie Fox, Ian Parker,
!    Chebyshev Polynomials in Numerical Analysis,
!    Oxford Press, page 56.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value at which the series
!    is to be evaluated.
!    X should be in the range -1.0 <= X <= 1.0.  The routine will refuse
!    to proceed if X <= -1.1 or 1.1 <= X.
!
!    Input, real ( kind = 8 ) CS(N), the coefficients of the Chebyshev series.
!    In evaluating CS, only half the first coef is summed.
!
!    Output, real ( kind = 8 ) VALUE, the value of the Chebyshev series at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) cs(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  Number of terms N <= 0.'
    stop
  end if

  if ( 1000 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  Number of terms 1000 < N.'
    stop
  end if

  if ( x < -1.1D+00 .or. 1.1D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSEVL - Fatal error!'
    write ( *, '(a)' ) '  X <= -1.1 or 1.1 <= X.'
    stop
  end if

  b1 = 0.0D+00
  b0 = 0.0D+00
  do i = 1, n
    b2 = b1
    b1 = b0
    b0 = 2.0D+00 * x * b1 - b2 + cs(n+1-i)
  end do

  value = 0.5D+00 * ( b0 - b2 )

  return
end
function error_f ( x )

!*****************************************************************************80
!
!! ERROR_F computes the error function.
!
!  Definition:
!
!    ERF(X) = ( 2 / SQRT ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( -T^2 ) dT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the error function.
!
!    Output, real ( kind = 8 ) ERROR_F, the value of the error function at X.
!
  implicit none

  real ( kind = 8 ) error_f
  real ( kind = 8 ) error_fc
  real ( kind = 8 ), parameter, dimension ( 13 ) :: erfcs = (/ &
    -0.049046121234691808D+00, -0.14226120510371364D+00, &
     0.010035582187599796D+00, -0.000576876469976748D+00, &
     0.000027419931252196D+00, -0.000001104317550734D+00, &
     0.000000038488755420D+00, -0.000000001180858253D+00, &
     0.000000000032334215D+00, -0.000000000000799101D+00, &
     0.000000000000017990D+00, -0.000000000000000371D+00, &
     0.000000000000000007D+00 /)
  integer ( kind = 4 ), save :: nterf = 0
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ), parameter :: sqrtpi = 1.7724538509055160D+00
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xbig = 0.0D+00
  real ( kind = 8 ) y
!
!  Initialize the Chebyshev series.
!
  if ( nterf == 0 ) then
    call inits ( erfcs, 13, 0.1D+00 * epsilon ( erfcs ), nterf )
    xbig = sqrt ( - log ( sqrtpi * epsilon ( xbig ) ) )
    sqeps = sqrt ( 2.0D+00 * epsilon ( sqeps ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then
    error_f = 2.0D+00 * x / sqrtpi
  else if ( y <= 1.0D+00 ) then
    call csevl ( 2.0D+00 * x**2 - 1.0D+00, erfcs, nterf, value )
    error_f = x * ( 1.0D+00 + value )
  else if ( y <= xbig ) then
    error_f = sign ( 1.0D+00 - error_fc ( y ), x )
  else
    error_f = sign ( 1.0D+00, x )
  end if

  return
end
function error_fc ( x )

!*****************************************************************************80
!
!! ERROR_FC computes the complementary error function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the error function.
!
!    Output, real ( kind = 8 ) ERROR_FC, the value of the complementary
!    error function at X.
!
  implicit none

  real ( kind = 8 ) error_fc
  real ( kind = 8 ), parameter, dimension ( 13 ) :: erfcs = (/ &
   -0.049046121234691808D+00, -0.14226120510371364D+00, &
    0.010035582187599796D+00, -0.000576876469976748D+00, &
    0.000027419931252196D+00, -0.000001104317550734D+00, &
    0.000000038488755420D+00, -0.000000001180858253D+00, &
    0.000000000032334215D+00, -0.000000000000799101D+00, &
    0.000000000000017990D+00, -0.000000000000000371D+00, &
    0.000000000000000007D+00 /)
  real ( kind = 8 ), parameter, dimension ( 24 ) :: erfccs = (/ &
    0.0715179310202925D+00, &
   -0.026532434337606719D+00, &
    0.001711153977920853D+00, &
   -0.000163751663458512D+00, &
    0.000019871293500549D+00, &
   -0.000002843712412769D+00, &
    0.000000460616130901D+00, &
   -0.000000082277530261D+00, &
    0.000000015921418724D+00, &
   -0.000000003295071356D+00, &
    0.000000000722343973D+00, &
   -0.000000000166485584D+00, &
    0.000000000040103931D+00, &
   -0.000000000010048164D+00, &
    0.000000000002608272D+00, &
   -0.000000000000699105D+00, &
    0.000000000000192946D+00, &
   -0.000000000000054704D+00, &
    0.000000000000015901D+00, &
   -0.000000000000004729D+00, &
    0.000000000000001432D+00, &
   -0.000000000000000439D+00, &
    0.000000000000000138D+00, &
   -0.000000000000000048D+00 /)
  real ( kind = 8 ), parameter, dimension ( 23 ) :: erc2cs = (/ &
   -0.069601346602309501D+00, &
   -0.041101339362620893D+00, &
    0.003914495866689626D+00, &
   -0.000490639565054897D+00, &
    0.000071574790013770D+00, &
   -0.000011530716341312D+00, &
    0.000001994670590201D+00, &
   -0.000000364266647159D+00, &
    0.000000069443726100D+00, &
   -0.000000013712209021D+00, &
    0.000000002788389661D+00, &
   -0.000000000581416472D+00, &
    0.000000000123892049D+00, &
   -0.000000000026906391D+00, &
    0.000000000005942614D+00, &
   -0.000000000001332386D+00, &
    0.000000000000302804D+00, &
   -0.000000000000069666D+00, &
    0.000000000000016208D+00, &
   -0.000000000000003809D+00, &
    0.000000000000000904D+00, &
   -0.000000000000000216D+00, &
    0.000000000000000052D+00 /)
  real ( kind = 8 ) eta
  integer ( kind = 4 ), save :: nterc2 = 0
  integer ( kind = 4 ), save :: nterf = 0
  integer ( kind = 4 ), save :: nterfc = 0
  real ( kind = 8 ), save :: sqeps = 0.0D+00
  real ( kind = 8 ), parameter :: sqrtpi = 1.7724538509055160D+00
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xmax = 0.0D+00
  real ( kind = 8 ), save :: xsml = 0.0D+00
  real ( kind = 8 ) y

  if ( nterf == 0 ) then

    eta = 0.1D+00 * epsilon ( eta )
    call inits ( erfcs, 13, eta, nterf )
    call inits ( erfccs, 24, eta, nterfc )
    call inits ( erc2cs, 23, eta, nterc2 )

    xsml = -sqrt ( - log ( sqrtpi * epsilon ( xsml ) ) )
    xmax = sqrt ( - log ( sqrtpi * tiny ( xmax ) ) )
    xmax = xmax - 0.5D+00 * log ( xmax ) / xmax - 0.01D+00
    sqeps = sqrt ( 2.0D+00 * epsilon ( sqeps ) )

  end if

  if ( x <= xsml ) then
    error_fc = 2.0D+00
    return
  end if
!
!  X so big that ERFC will underflow.
!
  if ( xmax < x ) then
    error_fc = 0.0D+00
    return
  end if

  y = abs ( x )
!
!  erfc(x) = 1.0D+00 - erf(x) for -1 <= x <= 1.
!
  if ( y <= 1.0D+00 ) then

    if ( y < sqeps ) then
      error_fc = 1.0D+00 - 2.0D+00 * x / sqrtpi
    else if ( sqeps <= y ) then
      call csevl ( 2.0D+00 * x * x - 1.0D+00, erfcs, nterf, value )
      error_fc = 1.0D+00 - x * ( 1.0D+00 + value )
    end if

    return

  end if
!
!  For 1 < |x| <= xmax, erfc(x) = 1.0D+00 - erf(x)
!
  y = y * y

  if ( y <= 4.0D+00 ) then
    call csevl ( ( 8.0D+00 / y - 5.0D+00 ) &
      / 3.0D+00, erc2cs, nterc2, value )
    error_fc = exp ( -y ) / abs ( x ) * ( 0.5D+00 + value )
  else
    call csevl ( 8.0D+00 / y - 1.0D+00, erfccs, nterfc, value )
    error_fc = exp ( -y ) / abs ( x ) * ( 0.5D+00 + value )
  end if

  if ( x < 0.0D+00 ) then
    error_fc = 2.0D+00 - error_fc
  end if

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
!      Gamma = limit ( M -> +oo ) ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 April 1999
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

  euler_constant = 0.577215664901532860606512090082402431042D+00

  return
end
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    However, when using the G95, the values returned by HUGE were
!    not equal to 2147483647, apparently, and were causing severe
!    and obscure errors in my random number generator, which needs to
!    add I4_HUGE to the seed whenever the seed is negative.  So I
!    am backing away from invoking HUGE, whereas I4_HUGE is under
!    my control.
!
!    Explanation: because under G95 the default integer type is 64 bits!
!    So HUGE ( 1 ) = a very very huge integer indeed, whereas
!    I4_HUGE ( ) = the same old 32 bit big value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end
subroutine i4_to_halton_number_sequence ( seed, base, n, r )

!*****************************************************************************80
!
!! I4_TO_HALTON_NUMBER_SEQUENCE: next N elements of a scalar Halton sequence.
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
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, pages 84-90, 1960.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED, the index of the desired element.
!    Only the absolute value of SEED is considered.
!    SEED = 0 is allowed, and returns R = 0.
!
!    Input, integer ( kind = 4 ) BASE, the Halton base, which should
!    be a prime number.  This routine only checks that BASE is greater
!    than 1.
!
!    Input, integer ( kind = 4 ) N, the number of elements desired.
!
!    Output, real ( kind = 8 ) R(N), the SEED-th through (SEED+N-1)-th
!    elements of the Halton sequence for base BASE.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  real ( kind = 8 ) base_inv
  integer ( kind = 4 ) digit(n)
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed2(n)
!
!  Set SEED2 = ( SEED, SEED+1, SEED+2, ..., SEED+N-1 )
!
  call i4vec_indicator ( n, seed2 )

  seed2(1:n) = seed2(1:n) + abs ( seed ) - 1

  if ( base <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON_NUMBER_SEQUENCE - Fatal error!'
    write ( *, '(a)' ) '  The input base BASE is <= 1!'
    write ( *, '(a,i8)' ) '  BASE = ', base
    stop
  end if

  base_inv = 1.0D+00 / real ( base )

  r(1:n) = 0.0D+00

  do while ( any ( seed2(1:n) /= 0 ) )
    digit(1:n) = mod ( seed2(1:n), base )
    r(1:n) = r(1:n) + real ( digit(1:n) ) * base_inv
    base_inv = base_inv / real ( base )
    seed2(1:n) = seed2(1:n) / base
  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector A(I)=I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
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
subroutine inits ( os, nos, eta, need )

!*****************************************************************************80
!
!! INITS finds the number of Chebyshev terms needed to achieve a given accuracy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2004
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) OS(NOS), the array of coefficients.
!
!    Input, integer ( kind = 4 ) NOS, the number of coefficients.
!
!    Input, real ( kind = 8 ) ETA, the requested accuracy.
!    A typical value of ETA is ( EPSILON ( 1.0 ) ) / 10.0.
!
!    Output, integer ( kind = 4 ) NEED, the number of terms needed for
!    the given accuracy.
!
  implicit none

  integer ( kind = 4 ) nos

  real ( kind = 8 ) err
  real ( kind = 8 ) eta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) need
  real ( kind = 8 ) os(nos)

  if ( nos < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INITS - Fatal error!'
    write ( *, '(a)' ) '  The number of coefficients is less than 1.'
    stop
  end if

  need = nos

  err = 0.0D+00
  do i = nos, 1, -1
    err = err + abs ( os(i) )
    if ( eta < err ) then
      need = i
      return
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INITS - Warning!'
  write ( *, '(a,g14.6)' ) '  The requested accuracy, ETA = ', eta
  write ( *, '(a)' ) '  is smaller than the Chebyshev coefficients'
  write ( *, '(a)' ) '  can guarantee.'

  return
end
subroutine p00_even ( prob, int_num, result )

!*****************************************************************************80
!
!! P00_EVEN uses evenly spaced points to integrate a function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) INT_NUM, the number of sample points.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  integer ( kind = 4 ) int_num

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fx(int_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) prob
  real ( kind = 8 ) result
  real ( kind = 8 ) x(int_num)

  call p00_lim ( prob, a, b )

  if ( int_num == 1 ) then
    x(1) = ( a + b ) / 2.0D+00
  else
    do i = 1, int_num
      x(i) = ( real ( int_num - i,     kind = 8 ) * a    &
             + real (           i - 1, kind = 8 ) * b  ) &
             / real ( int_num     - 1, kind = 8 )
    end do
  end if

  call p00_fun ( prob, int_num, x, fx )

  result = ( b - a ) * sum ( fx(1:int_num) ) &
    / real ( int_num, kind = 8 )

  return
end
subroutine p00_exact ( prob, exact )

!*****************************************************************************80
!
!! P00_EXACT returns the exact integral for any problem.
!
!  Discussion:
!
!    This routine provides a "generic" interface to the exact integral
!    routines for the various problems, and allows a problem to be called
!    by number (PROB) rather than by name.
!
!    In some cases, the "exact" value of the integral is in fact
!    merely a respectable approximation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the desired test problem.
!
!    Output, real ( kind = 8 ) EXACT, the exact value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_exact ( exact )
  else if ( prob == 2 ) then
    call p02_exact ( exact )
  else if ( prob == 3 ) then
    call p03_exact ( exact )
  else if ( prob == 4 ) then
    call p04_exact ( exact )
  else if ( prob == 5 ) then
    call p05_exact ( exact )
  else if ( prob == 6 ) then
    call p06_exact ( exact )
  else if ( prob == 7 ) then
    call p07_exact ( exact )
  else if ( prob == 8 ) then
    call p08_exact ( exact )
  else if ( prob == 9 ) then
    call p09_exact ( exact )
  else if ( prob == 10 ) then
    call p10_exact ( exact )
  else if ( prob == 11 ) then
    call p11_exact ( exact )
  else if ( prob == 12 ) then
    call p12_exact ( exact )
  else if ( prob == 13 ) then
    call p13_exact ( exact )
  else if ( prob == 14 ) then
    call p14_exact ( exact )
  else if ( prob == 15 ) then
    call p15_exact ( exact )
  else if ( prob == 16 ) then
    call p16_exact ( exact )
  else if ( prob == 17 ) then
    call p17_exact ( exact )
  else if ( prob == 18 ) then
    call p18_exact ( exact )
  else if ( prob == 19 ) then
    call p19_exact ( exact )
  else if ( prob == 20 ) then
    call p20_exact ( exact )
  else if ( prob == 21 ) then
    call p21_exact ( exact )
  else if ( prob == 22 ) then
    call p22_exact ( exact )
  else if ( prob == 23 ) then
    call p23_exact ( exact )
  else if ( prob == 24 ) then
    call p24_exact ( exact )
  else if ( prob == 25 ) then
    call p25_exact ( exact )
  else if ( prob == 26 ) then
    call p26_exact ( exact )
  else if ( prob == 27 ) then
    call p27_exact ( exact )
  else if ( prob == 28 ) then
    call p28_exact ( exact )
  else if ( prob == 29 ) then
    call p29_exact ( exact )
  else if ( prob == 30 ) then
    call p30_exact ( exact )
  else if ( prob == 31 ) then
    call p31_exact ( exact )
  else if ( prob == 32 ) then
    call p32_exact ( exact )
  else if ( prob == 33 ) then
    call p33_exact ( exact )
  else if ( prob == 34 ) then
    call p34_exact ( exact )
  else if ( prob == 35 ) then
    call p35_exact ( exact )
  else if ( prob == 36 ) then
    call p36_exact ( exact )
  else if ( prob == 37 ) then
    call p37_exact ( exact )
  else if ( prob == 38 ) then
    call p38_exact ( exact )
  else if ( prob == 39 ) then
    call p39_exact ( exact )
  else if ( prob == 40 ) then
    call p40_exact ( exact )
  else if ( prob == 41 ) then
    call p41_exact ( exact )
  else if ( prob == 42 ) then
    call p42_exact ( exact )
  else if ( prob == 43 ) then
    call p43_exact ( exact )
  else if ( prob == 44 ) then
    call p44_exact ( exact )
  else if ( prob == 45 ) then
    call p45_exact ( exact )
  else if ( prob == 46 ) then
    call p46_exact ( exact )
  else if ( prob == 47 ) then
    call p47_exact ( exact )
  else if ( prob == 48 ) then
    call p48_exact ( exact )
  else if ( prob == 49 ) then
    call p49_exact ( exact )
  else if ( prob == 50 ) then
    call p50_exact ( exact )
  else if ( prob == 51 ) then
    call p51_exact ( exact )
  else if ( prob == 52 ) then
    call p52_exact ( exact )
  else if ( prob == 53 ) then
    call p53_exact ( exact )
  else if ( prob == 54 ) then
    call p54_exact ( exact )
  else if ( prob == 55 ) then
    call p55_exact ( exact )
  else if ( prob == 56 ) then
    call p56_exact ( exact )
  else if ( prob == 57 ) then
    call p57_exact ( exact )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_EXACT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_fun ( prob, n, x, fx )

!*****************************************************************************80
!
!! P00_FUN evaluates the integrand for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the desired test problem.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) prob
  real ( kind = 8 ) x(n)

  if ( prob == 1 ) then
    call p01_fun ( n, x, fx )
  else if ( prob == 2 ) then
    call p02_fun ( n, x, fx )
  else if ( prob == 3 ) then
    call p03_fun ( n, x, fx )
  else if ( prob == 4 ) then
    call p04_fun ( n, x, fx )
  else if ( prob == 5 ) then
    call p05_fun ( n, x, fx )
  else if ( prob == 6 ) then
    call p06_fun ( n, x, fx )
  else if ( prob == 7 ) then
    call p07_fun ( n, x, fx )
  else if ( prob == 8 ) then
    call p08_fun ( n, x, fx )
  else if ( prob == 9 ) then
    call p09_fun ( n, x, fx )
  else if ( prob == 10 ) then
    call p10_fun ( n, x, fx )
  else if ( prob == 11 ) then
    call p11_fun ( n, x, fx )
  else if ( prob == 12 ) then
    call p12_fun ( n, x, fx )
  else if ( prob == 13 ) then
    call p13_fun ( n, x, fx )
  else if ( prob == 14 ) then
    call p14_fun ( n, x, fx )
  else if ( prob == 15 ) then
    call p15_fun ( n, x, fx )
  else if ( prob == 16 ) then
    call p16_fun ( n, x, fx )
  else if ( prob == 17 ) then
    call p17_fun ( n, x, fx )
  else if ( prob == 18 ) then
    call p18_fun ( n, x, fx )
  else if ( prob == 19 ) then
    call p19_fun ( n, x, fx )
  else if ( prob == 20 ) then
    call p20_fun ( n, x, fx )
  else if ( prob == 21 ) then
    call p21_fun ( n, x, fx )
  else if ( prob == 22 ) then
    call p22_fun ( n, x, fx )
  else if ( prob == 23 ) then
    call p23_fun ( n, x, fx )
  else if ( prob == 24 ) then
    call p24_fun ( n, x, fx )
  else if ( prob == 25 ) then
    call p25_fun ( n, x, fx )
  else if ( prob == 26 ) then
    call p26_fun ( n, x, fx )
  else if ( prob == 27 ) then
    call p27_fun ( n, x, fx )
  else if ( prob == 28 ) then
    call p28_fun ( n, x, fx )
  else if ( prob == 29 ) then
    call p29_fun ( n, x, fx )
  else if ( prob == 30 ) then
    call p30_fun ( n, x, fx )
  else if ( prob == 31 ) then
    call p31_fun ( n, x, fx )
  else if ( prob == 32 ) then
    call p32_fun ( n, x, fx )
  else if ( prob == 33 ) then
    call p33_fun ( n, x, fx )
  else if ( prob == 34 ) then
    call p34_fun ( n, x, fx )
  else if ( prob == 35 ) then
    call p35_fun ( n, x, fx )
  else if ( prob == 36 ) then
    call p36_fun ( n, x, fx )
  else if ( prob == 37 ) then
    call p37_fun ( n, x, fx )
  else if ( prob == 38 ) then
    call p38_fun ( n, x, fx )
  else if ( prob == 39 ) then
    call p39_fun ( n, x, fx )
  else if ( prob == 40 ) then
    call p40_fun ( n, x, fx )
  else if ( prob == 41 ) then
    call p41_fun ( n, x, fx )
  else if ( prob == 42 ) then
    call p42_fun ( n, x, fx )
  else if ( prob == 43 ) then
    call p43_fun ( n, x, fx )
  else if ( prob == 44 ) then
    call p44_fun ( n, x, fx )
  else if ( prob == 45 ) then
    call p45_fun ( n, x, fx )
  else if ( prob == 46 ) then
    call p46_fun ( n, x, fx )
  else if ( prob == 47 ) then
    call p47_fun ( n, x, fx )
  else if ( prob == 48 ) then
    call p48_fun ( n, x, fx )
  else if ( prob == 49 ) then
    call p49_fun ( n, x, fx )
  else if ( prob == 50 ) then
    call p50_fun ( n, x, fx )
  else if ( prob == 51 ) then
    call p51_fun ( n, x, fx )
  else if ( prob == 52 ) then
    call p52_fun ( n, x, fx )
  else if ( prob == 53 ) then
    call p53_fun ( n, x, fx )
  else if ( prob == 54 ) then
    call p54_fun ( n, x, fx )
  else if ( prob == 55 ) then
    call p55_fun ( n, x, fx )
  else if ( prob == 56 ) then
    call p56_fun ( n, x, fx )
  else if ( prob == 57 ) then
    call p57_fun ( n, x, fx )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FUN - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_gauss_legendre ( prob, int_num, result )

!*****************************************************************************80
!
!! P00_GAUSS_LEGENDRE applies a composite Gauss-Legendre rule.
!
!  Discussion:
!
!    A 4 point rule is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) INT_NUM, the number of subintervals.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  integer ( kind = 4 ), parameter :: gauss_num = 4

  real ( kind = 8 ) a
  real ( kind = 8 ) a_sub
  real ( kind = 8 ) b
  real ( kind = 8 ) b_sub
  real ( kind = 8 ) fx(gauss_num)
  real ( kind = 8 ), parameter :: gauss_abs(gauss_num) = (/ &
    -0.861136311594052575223946488893D+00, &
    -0.339981043584856264802665759103D+00, &
     0.339981043584856264802665759103D+00, &
     0.861136311594052575223946488893D+00 /)
  real ( kind = 8 ), parameter :: gauss_weight(gauss_num) = (/ &
    0.347854845137453857373063949222D+00, &
    0.652145154862546142626936050778D+00, &
    0.652145154862546142626936050778D+00, &
    0.347854845137453857373063949222D+00 /)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) int_i
  integer ( kind = 4 ) int_num
  integer ( kind = 4 ) prob
  real ( kind = 8 ) result
  real ( kind = 8 ) x(gauss_num)

  call p00_lim ( prob, a, b )

  h = ( b - a ) / real ( int_num, kind = 8 )

  result = 0.0D+00

  do int_i = 1, int_num

    a_sub = ( real ( int_num - int_i + 1, kind = 8 ) * a   &
            + real (           int_i - 1, kind = 8 ) * b ) &
            / real ( int_num,             kind = 8 )

    b_sub = ( real ( int_num - int_i, kind = 8 ) * a   &
            + real (           int_i, kind = 8 ) * b ) &
            / real ( int_num,         kind = 8 )

    x(1:gauss_num) = 0.5D+00 * ( ( 1.0D+00 - gauss_abs(1:gauss_num) ) * a_sub &
      + ( 1.0D+00 + gauss_abs(1:gauss_num) ) * b_sub )

    call p00_fun ( prob, gauss_num, x, fx )

    result = result + 0.5D+00 * h &
      * dot_product ( gauss_weight(1:gauss_num), fx(1:gauss_num) )

  end do

  return
end
subroutine p00_halton ( prob, int_num, result )

!*****************************************************************************80
!
!! P00_HALTON applies a Halton sequence rule to integrate a function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, pages 84-90, 1960.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) INT_NUM, the number of sample points.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  integer ( kind = 4 ) int_num

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) base
  real ( kind = 8 ) fx(int_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) prob
  real ( kind = 8 ) result
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(int_num)

  call p00_lim ( prob, a, b )

  seed = 1
  base = 2
  call i4_to_halton_number_sequence ( seed, base, int_num, x )
  x(1:int_num) = a + ( b - a ) * x(1:int_num)

  call p00_fun ( prob, int_num, x, fx )

  result = ( b - a ) * sum ( fx(1:int_num) ) &
    / real ( int_num, kind = 8 )

  return
end
subroutine p00_lim ( prob, a, b )

!*****************************************************************************80
!
!! P00_LIM returns the integration limits for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the desired test problem.
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_lim ( a, b )
  else if ( prob == 2 ) then
    call p02_lim ( a, b )
  else if ( prob == 3 ) then
    call p03_lim ( a, b )
  else if ( prob == 4 ) then
    call p04_lim ( a, b )
  else if ( prob == 5 ) then
    call p05_lim ( a, b )
  else if ( prob == 6 ) then
    call p06_lim ( a, b )
  else if ( prob == 7 ) then
    call p07_lim ( a, b )
  else if ( prob == 8 ) then
    call p08_lim ( a, b )
  else if ( prob == 9 ) then
    call p09_lim ( a, b )
  else if ( prob == 10 ) then
    call p10_lim ( a, b )
  else if ( prob == 11 ) then
    call p11_lim ( a, b )
  else if ( prob == 12 ) then
    call p12_lim ( a, b )
  else if ( prob == 13 ) then
    call p13_lim ( a, b )
  else if ( prob == 14 ) then
    call p14_lim ( a, b )
  else if ( prob == 15 ) then
    call p15_lim ( a, b )
  else if ( prob == 16 ) then
    call p16_lim ( a, b )
  else if ( prob == 17 ) then
    call p17_lim ( a, b )
  else if ( prob == 18 ) then
    call p18_lim ( a, b )
  else if ( prob == 19 ) then
    call p19_lim ( a, b )
  else if ( prob == 20 ) then
    call p20_lim ( a, b )
  else if ( prob == 21 ) then
    call p21_lim ( a, b )
  else if ( prob == 22 ) then
    call p22_lim ( a, b )
  else if ( prob == 23 ) then
    call p23_lim ( a, b )
  else if ( prob == 24 ) then
    call p24_lim ( a, b )
  else if ( prob == 25 ) then
    call p25_lim ( a, b )
  else if ( prob == 26 ) then
    call p26_lim ( a, b )
  else if ( prob == 27 ) then
    call p27_lim ( a, b )
  else if ( prob == 28 ) then
    call p28_lim ( a, b )
  else if ( prob == 29 ) then
    call p29_lim ( a, b )
  else if ( prob == 30 ) then
    call p30_lim ( a, b )
  else if ( prob == 31 ) then
    call p31_lim ( a, b )
  else if ( prob == 32 ) then
    call p32_lim ( a, b )
  else if ( prob == 33 ) then
    call p33_lim ( a, b )
  else if ( prob == 34 ) then
    call p34_lim ( a, b )
  else if ( prob == 35 ) then
    call p35_lim ( a, b )
  else if ( prob == 36 ) then
    call p36_lim ( a, b )
  else if ( prob == 37 ) then
    call p37_lim ( a, b )
  else if ( prob == 38 ) then
    call p38_lim ( a, b )
  else if ( prob == 39 ) then
    call p39_lim ( a, b )
  else if ( prob == 40 ) then
    call p40_lim ( a, b )
  else if ( prob == 41 ) then
    call p41_lim ( a, b )
  else if ( prob == 42 ) then
    call p42_lim ( a, b )
  else if ( prob == 43 ) then
    call p43_lim ( a, b )
  else if ( prob == 44 ) then
    call p44_lim ( a, b )
  else if ( prob == 45 ) then
    call p45_lim ( a, b )
  else if ( prob == 46 ) then
    call p46_lim ( a, b )
  else if ( prob == 47 ) then
    call p47_lim ( a, b )
  else if ( prob == 48 ) then
    call p48_lim ( a, b )
  else if ( prob == 49 ) then
    call p49_lim ( a, b )
  else if ( prob == 50 ) then
    call p50_lim ( a, b )
  else if ( prob == 51 ) then
    call p51_lim ( a, b )
  else if ( prob == 52 ) then
    call p52_lim ( a, b )
  else if ( prob == 53 ) then
    call p53_lim ( a, b )
  else if ( prob == 54 ) then
    call p54_lim ( a, b )
  else if ( prob == 55 ) then
    call p55_lim ( a, b )
  else if ( prob == 56 ) then
    call p56_lim ( a, b )
  else if ( prob == 57 ) then
    call p57_lim ( a, b )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_LIM - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_midpoint ( prob, int_num, result )

!*****************************************************************************80
!
!! P00_MIDPOINT applies the composite midpoint rule to integrate a function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) INT_NUM, the number of subintervals.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  integer ( kind = 4 ) int_num

  real ( kind = 8 ) a
  real ( kind = 8 ) a_sub
  real ( kind = 8 ) b
  real ( kind = 8 ) b_sub
  real ( kind = 8 ) fx(int_num)
  integer ( kind = 4 ) int_i
  integer ( kind = 4 ) prob
  real ( kind = 8 ) result
  real ( kind = 8 ) x(int_num)

  call p00_lim ( prob, a, b )

  do int_i = 1, int_num

    a_sub = ( real ( int_num - int_i + 1, kind = 8 ) * a   &
            + real (           int_i - 1, kind = 8 ) * b ) &
            / real ( int_num,             kind = 8 )

    b_sub = ( real ( int_num - int_i, kind = 8 ) * a   &
            + real (           int_i, kind = 8 ) * b ) &
            / real ( int_num,         kind = 8 )

    x(int_i) = 0.5D+00 * ( a_sub + b_sub )

  end do

  call p00_fun ( prob, int_num, x, fx )

  result =  ( b - a ) * sum ( fx(1:int_num) ) / real ( int_num, kind = 8 )

  return
end
subroutine p00_montecarlo ( prob, int_num, result )

!*****************************************************************************80
!
!! P00_MONTECARLO applies the Monte Carlo rule to integrate a function.
!
!  Discussion:
!
!    This routine originally used an automatic array for X.  However,
!    under the G95 compiler, this was causing bizarre errors.  Replacing
!    the automatic array by an allocatable array made the problems
!    disappear.  Not an entirely satisfactory conclusion!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) INT_NUM, the number of sample points.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  integer ( kind = 4 ) int_num

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: fx(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) prob
  real ( kind = 8 ) result
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:)

  seed = 123456789

  call p00_lim ( prob, a, b )

  allocate ( x(1:int_num) )
  allocate ( fx(1:int_num) )

  call r8vec_uniform ( int_num, a, b, seed, x )

  call p00_fun ( prob, int_num, x, fx )

  result = ( b - a ) * sum ( fx(1:int_num) ) &
    / real ( int_num, kind = 8 )

  deallocate ( fx )
  deallocate ( x )

  return
end
subroutine p00_prob_num ( prob_num )

!*****************************************************************************80
!
!! P00_PROB_NUM returns the number of test integration problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) PROB_NUM, the number of test integration
!    problems.
!
  implicit none

  integer ( kind = 4 ) prob_num

  prob_num = 57

  return
end
subroutine p00_simpson ( prob, int_num, result )

!*****************************************************************************80
!
!! P00_SIMPSON applies the composite Simpson rule to integrate a function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) INT_NUM, the number of subintervals.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a_sub
  real ( kind = 8 ) b
  real ( kind = 8 ) b_sub
  real ( kind = 8 ) fx1(1)
  real ( kind = 8 ) fx2(1)
  real ( kind = 8 ) fx3(1)
  real ( kind = 8 ) h
  integer ( kind = 4 ) int_i
  integer ( kind = 4 ) int_num
  integer ( kind = 4 ) prob
  real ( kind = 8 ) result
  real ( kind = 8 ) x1(1)
  real ( kind = 8 ) x2(1)
  real ( kind = 8 ) x3(1)

  call p00_lim ( prob, a, b )

  h = ( b - a ) / real ( int_num, kind = 8 )

  result = 0.0D+00

  do int_i = 1, int_num

    a_sub = ( real ( int_num - int_i + 1, kind = 8 ) * a   &
            + real (           int_i - 1, kind = 8 ) * b ) &
            / real ( int_num,             kind = 8 )

    b_sub = ( real ( int_num - int_i ) * a &
            + real (           int_i ) * b ) / real ( int_num, kind = 8 )

    x1(1) = a_sub
    call p00_fun ( prob, 1, x1, fx1 )
    x2(1) = 0.5D+00 * ( a_sub + b_sub )
    call p00_fun ( prob, 1, x2, fx2 )
    x3(1) = b_sub
    call p00_fun ( prob, 1, x3, fx3 )

    result = result + h * ( &
                     fx1(1) &
         + 4.0D+00 * fx2(1) &
         +           fx3(1) ) / 6.0D+00

  end do

  return
end
subroutine p00_trapezoid ( prob, int_num, result )

!*****************************************************************************80
!
!! P00_TRAPEZOID applies the composite trapezoid rule to integrate a function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the problem index.
!
!    Input, integer ( kind = 4 ) INT_NUM, the number of subintervals.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a_sub
  real ( kind = 8 ) b
  real ( kind = 8 ) b_sub
  real ( kind = 8 ) fx1(1)
  real ( kind = 8 ) fx2(1)
  real ( kind = 8 ) h
  integer ( kind = 4 ) int_i
  integer ( kind = 4 ) int_num
  integer ( kind = 4 ) prob
  real ( kind = 8 ) result
  real ( kind = 8 ) x1(1)
  real ( kind = 8 ) x2(1)

  call p00_lim ( prob, a, b )

  h = ( b - a ) / real ( int_num, kind = 8 )

  result = 0.0D+00

  do int_i = 1, int_num

    a_sub = ( real ( int_num - int_i + 1, kind = 8 ) * a   &
            + real (           int_i - 1, kind = 8 ) * b ) &
            / real ( int_num,             kind = 8 )

    b_sub = ( real ( int_num - int_i, kind = 8 ) * a   &
            + real (           int_i, kind = 8 ) * b ) &
            / real ( int_num,         kind = 8 )

    x1(1) = a_sub
    x2(1) = b_sub

    call p00_fun ( prob, 1, x1, fx1 )
    call p00_fun ( prob, 1, x2, fx2 )

    result = result + 0.5D+00 * h * ( fx1(1) + fx2(1) )

  end do

  return
end
subroutine p01_exact ( exact )

!*****************************************************************************80
!
!! P01_EXACT returns the exact integral for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = exp ( 1.0D+00 ) - 1.0D+00

  return
end
subroutine p01_fun ( n, x, fx )

!*****************************************************************************80
!
!! P01_FUN evaluates the integrand for problem 1.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    exp ( x )
!
!  Antiderivative:
!
!    exp ( x )
!
!  Exact Integral:
!
!    exp ( 1 ) - 1
!
!  Approximate Integral (25 digits):
!
!    1.718281828459045235360287...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) x(n)

  fx(1:n) = exp ( x(1:n) )

  return
end
subroutine p01_lim ( a, b )

!*****************************************************************************80
!
!! P01_LIM returns the integration limits for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p02_exact ( exact )

!*****************************************************************************80
!
!! P02_EXACT returns the exact integral for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.7D+00

  return
end
subroutine p02_fun ( n, x, fx )

!*****************************************************************************80
!
!! P02_FUN evaluates the integrand for problem 2.
!
!  Discussion:
!
!    The integrand is discontinuous at X = 0.3.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    if ( x < 0.3 )
!      f(x) = 0
!    else
!      f(x) = 1
!
!  Antiderivative:
!
!    if ( x < 0.3 )
!      g(x) = 0
!    else
!      g(x) = X - 0.3
!
!  Exact Integral:
!
!    0.7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) < 0.3D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = 1.0D+00
    end if

  end do

  return
end
subroutine p02_lim ( a, b )

!*****************************************************************************80
!
!! P02_LIM returns the integration limits for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p03_exact ( exact )

!*****************************************************************************80
!
!! P03_EXACT returns the exact integral for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 2.0D+00 / 3.0D+00

  return
end
subroutine p03_fun ( n, x, fx )

!*****************************************************************************80
!
!! P03_FUN evaluates the integrand for problem 3.
!
!  Discussion:
!
!    The integrand is not differentiable at X = 0.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    sqrt ( x )
!
!  Antiderivative:
!
!    ( 2 / 3 ) * x^(3/2)
!
!  Exact Integral:
!
!    2 / 3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = sqrt ( x(1:n) )

  return
end
subroutine p03_lim ( a, b )

!*****************************************************************************80
!
!! P03_LIM returns the integration limits for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p04_exact ( exact )

!*****************************************************************************80
!
!! P04_EXACT returns the estimated integral for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.47942822668880166736D+00

  return
end
subroutine p04_fun ( n, x, fx )

!*****************************************************************************80
!
!! P04_FUN evaluates the integrand for problem 4.
!
!  Interval:
!
!    -1 <= x <= 1
!
!  Integrand:
!
!    0.92 * cosh ( x ) - cos ( x )
!
!  Antiderivative:
!
!    0.92 * sinh ( x ) - sin ( x )
!
!  Exact Integral:
!
!    1.84 * sinh ( 1 ) - 2 * sin ( 1 )
!
!  Approximate Integral (20 digits):
!
!    0.47942822668880166736...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Clenshaw, Alan Curtis,
!    A Method for Numerical Integration on an Automatic Computer,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 197-205.
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = 0.92D+00 * cosh ( x(1:n) ) - cos ( x(1:n) )

  return
end
subroutine p04_lim ( a, b )

!*****************************************************************************80
!
!! P04_LIM returns the integration limits for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = -1.0D+00
  b = 1.0D+00

  return
end
subroutine p05_exact ( exact )

!*****************************************************************************80
!
!! P05_EXACT returns the estimated integral for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 1.5822329637296729331D+00

  return
end
subroutine p05_fun ( n, x, fx )

!*****************************************************************************80
!
!! P05_FUN evaluates the integrand for problem 5.
!
!  Interval:
!
!    -1 <= x <= 1
!
!  Integrand:
!
!    1 / ( x^4 + x^2 + 0.9 )
!
!  Approximate Integral (20 digits):
!
!    1.5822329637296729331...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Clenshaw, Alan Curtis,
!    A Method for Numerical Integration on an Automatic Computer,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 197-205.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = 1.0D+00 / ( x(1:n)**4 + x(1:n)**2 + 0.9D+00 )

  return
end
subroutine p05_lim ( a, b )

!*****************************************************************************80
!
!! P05_LIM returns the integration limits for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = -1.0D+00
  b = 1.0D+00

  return
end
subroutine p06_exact ( exact )

!*****************************************************************************80
!
!! P06_EXACT returns the exact integral for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 1.460447131787105D+00

  return
end
subroutine p06_fun ( n, x, fx )

!*****************************************************************************80
!
!! P06_FUN evaluates the integrand for problem 6.
!
!  Interval:
!
!    -1 <= x <= 1
!
!  Integrand:
!
!    sqrt ( abs ( x + 0.5 ) )
!
!  Exact Integral:
!
!    ( sqrt ( 2 ) + 3 * sqrt ( 6 ) ) / 6 = 1.460447131787105
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Clenshaw, Alan Curtis,
!    A Method for Numerical Integration on an Automatic Computer,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 197-205.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = sqrt ( abs ( x(1:n) + 0.5D+00 ) )

  return
end
subroutine p06_lim ( a, b )

!*****************************************************************************80
!
!! P06_LIM returns the integration limits for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = -1.0D+00
  b = 1.0D+00

  return
end
subroutine p07_exact ( exact )

!*****************************************************************************80
!
!! P07_EXACT returns the exact integral for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 2.0D+00

  return
end
subroutine p07_fun ( n, x, fx )

!*****************************************************************************80
!
!! P07_FUN evaluates the integrand for problem 7.
!
!  Discussion:
!
!    The integrand is singular at X = 0.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    1 / sqrt ( X )
!
!  Antiderivative:
!
!    2 * sqrt ( X )
!
!  Exact Integral:
!
!    2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( 0.0D+00 < x(i) ) then
      fx(i) = 1.0D+00 / sqrt ( x(i) )
    else
      fx(i) = 0.0D+00
    end if

  end do

  return
end
subroutine p07_lim ( a, b )

!*****************************************************************************80
!
!! P07_LIM returns the integration limits for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p08_exact ( exact )

!*****************************************************************************80
!
!! P08_EXACT returns the estimated integral for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.86697298733991103757D+00

  return
end
subroutine p08_fun ( n, x, fx )

!*****************************************************************************80
!
!! P08_FUN evaluates the integrand for problem 8.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    1 / ( 1 + X^4 )
!
!  Antiderivative:
!
!    (1/8) * sqrt ( 2 )
!    * ln ( ( X^2 + sqrt ( 2 ) * X + 1 ) / ( X^2 - sqrt ( 2 ) * X + 1 ) )
!    + (1/4) * sqrt ( 2 ) * arctan ( sqrt ( 2 ) * X / ( 1 - X^2 ) )
!
!  Approximate Integral (20 digits):
!
!    0.86697298733991103757...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = 1.0D+00 / ( 1.0D+00 + x(1:n)**4 )

  return
end
subroutine p08_lim ( a, b )

!*****************************************************************************80
!
!! P08_LIM returns the integration limits for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p09_exact ( exact )

!*****************************************************************************80
!
!! P09_EXACT returns the estimated integral for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 1.1547005383792515290D+00

  return
end
subroutine p09_fun ( n, x, fx )

!*****************************************************************************80
!
!! P09_FUN evaluates the integrand for problem 9.
!
!  Discussion:
!
!    The integrand is oscillatory, going through 5 periods in [0,1].
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    2 / ( 2 + sin ( 10 * pi * X ) )
!
!  Antiderivative:
!
!    1 / ( 5 * pi * sqrt ( 3 ) ) *
!    arctan ( ( 1 + 2 * tan ( 5 * pi * X ) ) / sqrt ( 3 ) )
!
!  Exact Integral:
!
!    2 / sqrt ( 3 )
!
!  Approximate Integral (20 digits):
!
!    1.1547005383792515290...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  fx(1:n) = 2.0D+00 / ( 2.0D+00 + sin ( 10.0D+00 * pi * x(1:n) ) )

  return
end
subroutine p09_lim ( a, b )

!*****************************************************************************80
!
!! P09_LIM returns the integration limits for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p10_exact ( exact )

!*****************************************************************************80
!
!! P10_EXACT returns the estimated integral for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.6931471805599453094172321D+00

  return
end
subroutine p10_fun ( n, x, fx )

!*****************************************************************************80
!
!! P10_FUN evaluates the integrand for problem 10.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    1 / ( 1 + X )
!
!  Antiderivative:
!
!    ln ( 1 + X )
!
!  Exact Integral:
!
!    ln ( 2 )
!
!  Approximate Integral (25 digits):
!
!    0.6931471805599453094172321...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = 1.0D+00 / ( 1.0D+00 + x(1:n) )

  return
end
subroutine p10_lim ( a, b )

!*****************************************************************************80
!
!! P10_LIM returns the integration limits for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p11_exact ( exact )

!*****************************************************************************80
!
!! P11_EXACT returns the estimated integral for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.37988549304172247537D+00

  return
end
subroutine p11_fun ( n, x, fx )

!*****************************************************************************80
!
!! P11_FUN evaluates the integrand for problem 11.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    1 / ( 1 + exp ( X ) )
!
!  Antiderivative:
!
!    ln ( exp ( X ) / ( 1 + exp ( X ) ) )
!
!  Exact Integral:
!
!    ln ( 2 * E / ( 1 + E ) )
!
!  Approximate Integral (20 digits):
!
!    0.37988549304172247537...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = 1.0D+00 / ( 1.0D+00 + exp ( x(1:n) ) )

  return
end
subroutine p11_lim ( a, b )

!*****************************************************************************80
!
!! P11_LIM returns the integration limits for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p12_exact ( exact )

!*****************************************************************************80
!
!! P12_EXACT returns the estimated integral for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.77750463411224827642D+00

  return
end
subroutine p12_fun ( n, x, fx )

!*****************************************************************************80
!
!! P12_FUN evaluates the integrand for problem 12.
!
!  Discussion:
!
!    The integrand has a removable singularity at X = 0.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    X / ( exp ( X ) - 1 )
!
!  Antiderivative:
!
!    The Debye function.
!
!  Approximate Integral (20 digits):
!
!    0.77750463411224827642...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 1.0D+00
    else
      fx(i) = x(i) / ( exp ( x(i) ) - 1.0D+00 )
    end if

  end do

  return
end
subroutine p12_lim ( a, b )

!*****************************************************************************80
!
!! P12_LIM returns the integration limits for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p13_exact ( exact )

!*****************************************************************************80
!
!! P13_EXACT returns the estimated integral for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) exact
  real ( kind = 8 ) r8_si

  call p13_lim ( a, b )

  exact = r8_si ( b ) - r8_si ( a )

  return
end
subroutine p13_fun ( n, x, fx )

!*****************************************************************************80
!
!! P13_FUN evaluates the integrand for problem 13.
!
!  Interval:
!
!    0 <= x <= 10
!
!  Integrand:
!
!    sin ( X ) / X
!
!  Approximate Integral (20 digits):
!
!    1.6583475942188740493...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 1.0D+00
    else
      fx(i) = sin ( x(i) ) / x(i)
    end if

  end do

  return
end
subroutine p13_lim ( a, b )

!*****************************************************************************80
!
!! P13_LIM returns the integration limits for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 10.0D+00

  return
end
subroutine p14_exact ( exact )

!*****************************************************************************80
!
!! P14_EXACT returns the estimated integral for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.500000211166D+00

  return
end
subroutine p14_fun ( n, x, fx )

!*****************************************************************************80
!
!! P14_FUN evaluates the integrand for problem 14.
!
!  Discussion:
!
!    For X's that aren't actually very big, the function becomes very
!    small.  Some compilers may product code that fails in these cases.
!    An attempt has been made to return a value of 0 when the computed
!    value of F(X) would be extremely small.
!
!  Interval:
!
!    0 <= x <= 10
!
!  Integrand:
!
!    sqrt ( 50 ) * exp ( - 50 * pi * x * x )
!
!  Exact Integral:
!
!    0.5 * erf ( 50 * sqrt ( 2 * pi ) )
!
!  Approximate Integral:
!
!    0.500000211166...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), save :: x_max = 0.0D+00

  if ( x_max == 0.0D+00 ) then
    x_max = sqrt ( log ( max ( epsilon ( x_max ), 1.0D-10 ) ) &
     / ( - 50.0D+00 * pi ) )
  end if

  do i = 1, n

    if ( x_max < abs ( x(i) ) ) then
      fx(i) = 0.0D+00
    else
      fx(i) = sqrt ( 50.0D+00 ) * exp ( - 50.0D+00 * pi * x(i) * x(i) )
    end if

  end do

  return
end
subroutine p14_lim ( a, b )

!*****************************************************************************80
!
!! P14_LIM returns the integration limits for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 10.0D+00

  return
end
subroutine p15_exact ( exact )

!*****************************************************************************80
!
!! P15_EXACT returns the exact integral for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 1.0D+00

  return
end
subroutine p15_fun ( n, x, fx )

!*****************************************************************************80
!
!! P15_FUN evaluates the integrand for problem 15.
!
!  Interval:
!
!    0 <= x <= 10
!
!  Integrand:
!
!    25 * exp ( - 25 * X )
!
!  Antiderivative:
!
!    - exp ( - 25 * X )
!
!  Exact Integral:
!
!    1 - exp ( - 250 )
!
!  Approximate Integral:
!
!    1.00000000...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) x(n)

  fx(1:n) = 25.0D+00 * exp ( - 25.0D+00 * x(1:n) )

  return
end
subroutine p15_lim ( a, b )

!*****************************************************************************80
!
!! P15_LIM returns the integration limits for problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 10.0D+00

  return
end
subroutine p16_exact ( exact )

!*****************************************************************************80
!
!! P16_EXACT returns the exact integral for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.49936338107645674464D+00

  return
end
subroutine p16_fun ( n, x, fx )

!*****************************************************************************80
!
!! P16_FUN evaluates the integrand for problem 16.
!
!  Interval:
!
!    0 <= x <= 10
!
!  Integrand:
!
!    50.0 / ( pi * ( 2500.0 * X * X + 1.0 ) )
!
!  Antiderivative:
!
!    ( 1 / pi ) * arctan ( 50 * X )
!
!  Approximate Integral (20 digits):
!
!    0.49936338107645674464...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)

  fx(1:n) = 50.0D+00 / pi / ( 2500.0D+00 * x(1:n) * x(1:n) + 1.0D+00 )

  return
end
subroutine p16_lim ( a, b )

!*****************************************************************************80
!
!! P16_LIM returns the integration limits for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p17_exact ( exact )

!*****************************************************************************80
!
!! P17_EXACT returns the estimated integral for problem 17.
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
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.5D+00

  return
end
subroutine p17_fun ( n, x, fx )

!*****************************************************************************80
!
!! P17_FUN evaluates the integrand for problem 17.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    ( sin ( 50 * pi * X ) )^2
!
!  Antiderivative:
!
!    1/2 X - sin ( 100 * pi * X ) / ( 200 * pi )
!
!  Approximate Integral:
!
!    0.5
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
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  fx(1:n) = ( sin ( 50.0D+00 * pi * x(1:n) ) )**2

  return
end
subroutine p17_lim ( a, b )

!*****************************************************************************80
!
!! P17_LIM returns the integration limits for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p18_exact ( exact )

!*****************************************************************************80
!
!! P18_EXACT returns the estimated integral for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.17055734950243820437D+00

  return
end
subroutine p18_fun ( n, x, fx )

!*****************************************************************************80
!
!! P18_FUN evaluates the integrand for problem 18.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    X / ( exp ( X ) + 1 )
!
!  Approximate Integral (20 digits):
!
!    0.17055734950243820437...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hermann Engels,
!    Numerical Quadrature and Cubature,
!    Academic Press, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = x(1:n) / ( exp ( x(1:n) ) + 1.0D+00 )

  return
end
subroutine p18_lim ( a, b )

!*****************************************************************************80
!
!! P18_LIM returns the integration limits for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p19_exact ( exact )

!*****************************************************************************80
!
!! P19_EXACT returns the exact integral for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = - 1.0D+00

  return
end
subroutine p19_fun ( n, x, fx )

!*****************************************************************************80
!
!! P19_FUN evaluates the integrand for problem 19.
!
!  Discussion:
!
!    The integrand is singular at X = 0.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    ln ( X )
!
!  Antiderivative:
!
!    X * ln ( X ) - X
!
!  Exact Integral:
!
!    -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) <= 1.0D-15 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = log ( x(i) )
    end if

  end do

  return
end
subroutine p19_lim ( a, b )

!*****************************************************************************80
!
!! P19_LIM returns the integration limits for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p20_exact ( exact )

!*****************************************************************************80
!
!! P20_EXACT returns the estimated integral for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  real ( kind = 8 ) exact

  exact = 1.5643964440690497731D+00

  return
end
subroutine p20_fun ( n, x, fx )

!*****************************************************************************80
!
!! P20_FUN evaluates the integrand for problem 20.
!
!  Interval:
!
!    -1 <= x <= 1
!
!  Integrand:
!
!    1 / ( X^2 + 1.005 )
!
!  Antiderivative:
!
!    ( 1 / sqrt ( 1.005 ) ) * arctan ( X / sqrt ( 1.005 ) )
!
!  Approximate Integral (20 digits):
!
!    1.5643964440690497731...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = 1.0D+00 / ( x(1:n)**2 + 1.005D+00 )

  return
end
subroutine p20_lim ( a, b )

!*****************************************************************************80
!
!! P20_LIM returns the integration limits for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = -1.0D+00
  b = 1.0D+00

  return
end
subroutine p21_exact ( exact )

!*****************************************************************************80
!
!! P21_EXACT returns the estimated integral for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.21080273631018169851D+00

  return
end
subroutine p21_fun ( n, x, fx )

!*****************************************************************************80
!
!! P21_FUN evaluates the integrand for problem 21.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!       ( sech (   10.0 * ( x - 0.2 ) ) )^2
!     + ( sech (  100.0 * ( x - 0.4 ) ) )^4
!     + ( sech ( 1000.0 * ( x - 0.6 ) ) )^6
!
!  Exact Integral:
!
!    ( 1 + tanh ( 8 ) * tanh ( 2 ) ) / 10.0 + 2 / 150 + 2 / 1875
!
!  Approximate Integral (20 digits):
!
!    0.21080273631018169851...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) x(n)

  fx(1:n) = &
      ( 1.0D+00 / cosh (   10.0D+00 * ( x(1:n) - 0.2D+00 ) ) )**2 &
    + ( 1.0D+00 / cosh (  100.0D+00 * ( x(1:n) - 0.4D+00 ) ) )**4 &
    + ( 1.0D+00 / cosh ( 1000.0D+00 * ( x(1:n) - 0.6D+00 ) ) )**6

  return
end
subroutine p21_lim ( a, b )

!*****************************************************************************80
!
!! P21_LIM returns the integration limits for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p22_exact ( exact )

!*****************************************************************************80
!
!! P22_EXACT returns the estimated integral for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = 0.125D+00 * log ( 9.0D+00 ) + pi / sqrt ( 48.0D+00 )

  return
end
subroutine p22_fun ( n, x, fx )

!*****************************************************************************80
!
!! P22_FUN evaluates the integrand for problem 22.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    1 / ( X^4 + X^2 + 1 )
!
!  Exact integral:
!
!    ln ( 9 ) / 8 + pi / sqrt ( 48 )
!
!  Approximate Integral (20 digits):
!
!    0.72810291322558188550...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = 1.0D+00 / ( x(1:n)**4 + x(1:n)**2 + 1.0D+00 )

  return
end
subroutine p22_lim ( a, b )

!*****************************************************************************80
!
!! P22_LIM returns the integration limits for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p23_exact ( exact )

!*****************************************************************************80
!
!! P23_EXACT returns the estimated integral for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.62471325642771360429D+00

  return
end
subroutine p23_fun ( n, x, fx )

!*****************************************************************************80
!
!! P23_FUN evaluates the integrand for problem 23.
!
!  Discussion:
!
!    The integrand has a singularity at X = 0.
!    The integrand is discontinuous at X = 0.
!    The integrand is arbitrarily oscillatory as X decreases to 0.
!    The integrand becomes unbounded as X decreases to 0.
!
!    Integral ( 0 < X < 1 ) ( 1 / X ) sin ( 1 / X ) dX
!    = Integral ( 1 < X < +oo ) ( 1 / X ) * sin ( X ) dX.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    ( 1 / x ) sin ( 1 / x )
!
!  Approximate Integral (20 digits):
!
!    0.62471325642771360429...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = ( 1.0D+00 / x(i) ) * sin ( 1.0D+00 / x(i) )
    end if

  end do

  return
end
subroutine p23_lim ( a, b )

!*****************************************************************************80
!
!! P23_LIM returns the integration limits for problem 23.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p24_exact ( exact )

!*****************************************************************************80
!
!! P24_EXACT returns the estimated integral for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = - 0.0067547455D+00

  return
end
subroutine p24_fun ( n, x, fx )

!*****************************************************************************80
!
!! P24_FUN evaluates the integrand for problem 24.
!
!  Discussion:
!
!    The integrand is continuous, but nowhere differentiable.
!
!  Interval:
!
!    0 <= X <= 0.5
!
!  Integrand:
!
!    ( 1 / pi ) * sum ( 1 <= I < +oo ) 2^(-I) * cos ( 7^I * pi * X )
!
!  Approximate Integral:
!
!    - 0.0067547455
!
!  Antiderivative:
!
!    ( 1 / pi^2 ) * sum ( 1 <= I < +oo ) 14^(-I) * sin ( 7^I * pi * X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Herbert Salzer, Norman Levine,
!    Table of a Weierstrass Continuous Nondifferentiable Function,
!    Mathematics of Computation,
!    Volume 15, pages 120 - 130, 1961.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ), parameter :: n_term = 40
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  do i = 1, n

    fx(i) = 0.0D+00
    do j = 1, n_term
      fx(i) = fx(i) + cos ( 7.0D+00**j * pi * x(i) ) / 2.0D+00**j
    end do

    fx(i) = fx(i) / pi

  end do

  return
end
subroutine p24_lim ( a, b )

!*****************************************************************************80
!
!! P24_LIM returns the integration limits for problem 24.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 0.5D+00

  return
end
subroutine p25_exact ( exact )

!*****************************************************************************80
!
!! P25_EXACT returns the estimated integral for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.3D+00 * log ( 0.3D+00 ) + 0.7D+00 * log ( 0.7D+00 ) - 1.0D+00

  return
end
subroutine p25_fun ( n, x, fx )

!*****************************************************************************80
!
!! P25_FUN evaluates the integrand for problem 25.
!
!  Interval:
!
!    0 <= X <= 1.
!
!  Integrand:
!
!    ln ( abs ( x - 0.7 ) )
!
!  Exact Integral:
!
!    0.3 * ln ( 0.3 ) + 0.7 * ln ( 0.7 ) - 1
!
!  Approximate Integral (20 digits):
!
!    -1.6108643020548934630
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kendall Atkinson,
!    An Introduction to Numerical Analysis,
!    Prentice Hall, 1984, page 303.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 0.7D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = log ( abs ( x(i) - 0.7D+00 ) )
    end if

  end do

  return
end
subroutine p25_lim ( a, b )

!*****************************************************************************80
!
!! P25_LIM returns the integration limits for problem 25.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p26_exact ( exact )

!*****************************************************************************80
!
!! P26_EXACT returns the exact integral for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 7.9549265210128452745D+00

  return
end
subroutine p26_fun ( n, x, fx )

!*****************************************************************************80
!
!! P26_FUN evaluates the integrand for problem 26.
!
!  Interval:
!
!    0 <= x <= 2 pi
!
!  Integrand:
!
!    exp ( cos ( x ) )
!
!  Exact Integral:
!
!    2 * Pi * I0(1)
!
!  Approximate Integral (20 digits):
!
!    7.9549265210128452745...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kendall Atkinson,
!    An Introduction to Numerical Analysis,
!    Prentice Hall, 1984, page 262.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = exp ( cos ( x(1:n) ) )

  return
end
subroutine p26_lim ( a, b )

!*****************************************************************************80
!
!! P26_LIM returns the integration limits for problem 26.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  a = 0.0D+00
  b = 2.0D+00 * pi

  return
end
subroutine p27_exact ( exact )

!*****************************************************************************80
!
!! P27_EXACT returns the exact integral for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 5.0D+00 - 6.0D+00 * log ( 2.0D+00 )

  return
end
subroutine p27_fun ( n, x, fx )

!*****************************************************************************80
!
!! P27_FUN evaluates the integrand for problem 27.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    1 / ( X^(1/2) + X^(1/3) )
!
!  Exact Integral:
!
!    5 - 6 * ln ( 2 )
!
!  Approximate Integral (20 digits):
!
!    0.84111691664032814350...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = 1.0D+00 / ( sqrt ( x(i) ) + x(i)**(1.0D+00/3.0D+00) )
    end if

  end do

  return
end
subroutine p27_lim ( a, b )

!*****************************************************************************80
!
!! P27_LIM returns the integration limits for problem 27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p28_exact ( exact )

!*****************************************************************************80
!
!! P28_EXACT returns the exact integral for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = ( 50.0D+00 / 2501.0D+00 ) * ( 1.0D+00 - exp ( - 2.0D+00 * pi ) )

  return
end
subroutine p28_fun ( n, x, fx )

!*****************************************************************************80
!
!! P28_FUN evaluates the integrand for problem 28.
!
!  Interval:
!
!    0 <= X <= 2 PI
!
!  Integrand:
!
!    exp ( - X ) * sin ( 50 * X )
!
!  Exact Integral:
!
!    50 / ( 2501 ) * ( 1 - exp ( - 2 * PI ) )
!
!  Approximate Integral (20 digits):
!
!    0.019954669277654778312...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kendall Atkinson,
!    An Introduction to Numerical Analysis,
!    Prentice Hall, 1984, page 303.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = exp ( - x(1:n) ) * sin ( 50.0D+00 * x(1:n) )

  return
end
subroutine p28_lim ( a, b )

!*****************************************************************************80
!
!! P28_LIM returns the integration limits for problem 28.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  a = 0.0D+00
  b = 2.0D+00 * pi

  return
end
subroutine p29_exact ( exact )

!*****************************************************************************80
!
!! P29_EXACT returns the exact integral for problem 29.
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
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 1.0D+00 - log ( 2.0D+00 )

  return
end
subroutine p29_fun ( n, x, fx )

!*****************************************************************************80
!
!! P29_FUN evaluates the integrand for problem 29.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    F ( X ) = 1 / ( X + 2 )   for 0 < X < E - 2
!            = 0               otherwise
!
!  Exact Integral:
!
!    1 - ln ( 2 )
!
!  Approximate Integral (20 digits):
!
!    0.30685281944005469058...
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
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( 0.0D+00 <= x(i) .and. x(i) <= exp ( 1.0D+00 ) - 2.0D+00 ) then
      fx(i) = 1.0D+00 / ( x(i) + 2.0D+00 )
    else
      fx(i) = 0.0D+00
    end if

  end do

  return
end
subroutine p29_lim ( a, b )

!*****************************************************************************80
!
!! P29_LIM returns the integration limits for problem 29.
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
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p30_exact ( exact )

!*****************************************************************************80
!
!! P30_EXACT returns the exact integral for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = -4.5275696251606720278D+00

  return
end
subroutine p30_fun ( n, x, fx )

!*****************************************************************************80
!
!! P30_FUN evaluates the integrand for problem 30.
!
!  Interval:
!
!    2 <= x <= 7
!
!  Integrand:
!
!          cos (       x )
!    + 5 * cos ( 1.6 * x )
!    - 2 * cos ( 2.0 * x )
!    + 5 * cos ( 4.5 * x )
!    + 7 * cos ( 9.0 * x )
!
!  Antiderivative:
!
!          sin (       x )
!    + 5 * sin ( 1.6 * x ) / 1.6
!    - 2 * sin ( 2.0 * x ) / 2.0
!    + 5 * sin ( 4.5 * x ) / 4.5
!    + 7 * sin ( 9.0 * x ) / 9.0
!
!  Exact Integral:
!
!    -4.5275696251606720278
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dianne OLeary,
!    Scientific Computing with Case Studies,
!    SIAM, 2008,
!    ISBN13: 978-0-898716-66-5,
!    LC: QA401.O44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  real ( kind = 8 ) x(n)

  fx(1:n) = &
                cos (           x(1:n) ) &
    + 5.0D+00 * cos ( 1.6D+00 * x(1:n) ) &
    - 2.0D+00 * cos ( 2.0D+00 * x(1:n) ) &
    + 5.0D+00 * cos ( 4.5D+00 * x(1:n) ) &
    + 7.0D+00 * cos ( 9.0D+00 * x(1:n) )

  return
end
subroutine p30_lim ( a, b )

!*****************************************************************************80
!
!! P30_LIM returns the integration limits for problem 30.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 2.0D+00
  b = 7.0D+00

  return
end
subroutine p31_exact ( exact )

!*****************************************************************************80
!
!! P31_EXACT returns the exact integral for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 2.0D+00 * atan ( 4.0D+00 )

  return
end
subroutine p31_fun ( n, x, fx )

!*****************************************************************************80
!
!! P31_FUN evaluates the integrand for problem 31.
!
!  Discussion:
!
!    A simple Newton-Cotes quadrature rule, in which the order of the
!    rule is increased, but the interval is not subdivided, diverges
!    for this integrand.
!
!    This is Runge's function, a standard example of the perils of
!    using high order polynomial interpolation at equally spaced nodes.
!    Since this is exactly what is really going on in a Newton Cotes
!    rule, it is little wonder that the result is so poor.
!
!  Interval:
!
!    -4 <= x <= 4
!
!  Integrand:
!
!    1 / ( 1 + x^2 )
!
!  Antiderivative:
!
!    arctan ( x )
!
!  Exact Integral:
!
!    2 * arctan ( 4 )
!
!  Approximate Integral (20 digits):
!
!    2.6516353273360649301...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kendall Atkinson,
!    An Introduction to Numerical Analysis,
!    Prentice Hall, 1984, page 266.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = 1.0D+00 / ( 1.0D+00 + x(1:n)**2 )

  return
end
subroutine p31_lim ( a, b )

!*****************************************************************************80
!
!! P31_LIM returns the integration limits for problem 31.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = - 4.0D+00
  b =   4.0D+00

  return
end
subroutine p32_exact ( exact )

!*****************************************************************************80
!
!! P32_EXACT returns the exact integral for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = - 0.5D+00 * ( exp ( pi ) + 1.0D+00 )

  return
end
subroutine p32_fun ( n, x, fx )

!*****************************************************************************80
!
!! P32_FUN evaluates the integrand for problem 32.
!
!  Interval:
!
!    0 <= X <= PI
!
!  Integrand:
!
!    exp ( X ) * cos ( X )
!
!  Antiderivative:
!
!    0.5 * exp ( X ) * ( sin ( X ) + cos ( X ) )
!
!  Exact Integral:
!
!    - 0.5 * ( exp ( PI ) + 1 )
!
!  Approximate Integral (20 digits):
!
!    -12.070346316389634503...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kendall Atkinson,
!    An Introduction to Numerical Analysis,
!    Prentice Hall, 1984, page 254, 277, 297.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = exp ( x(1:n) ) * cos ( x(1:n) )

  return
end
subroutine p32_lim ( a, b )

!*****************************************************************************80
!
!! P32_LIM returns the integration limits for problem 32.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  a = 0.0D+00
  b = pi

  return
end
subroutine p33_exact ( exact )

!*****************************************************************************80
!
!! P33_EXACT returns the exact integral for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = 0.5D+00 * sqrt ( pi )

  return
end
subroutine p33_fun ( n, x, fx )

!*****************************************************************************80
!
!! P33_FUN evaluates the integrand for problem 33.
!
!  Discussion:
!
!    The integrand is singular at both endpoints of the interval.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    sqrt ( - ln ( X ) )
!
!  Exact Integral:
!
!    sqrt ( pi ) / 2
!
!  Approximate Integral (20 digits):
!
!    0.88622692545275801365...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kendall Atkinson,
!    An Introduction to Numerical Analysis,
!    Prentice Hall, 1984, page 307.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) <= 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = sqrt ( - log ( x(i) ) )
    end if

  end do

  return
end
subroutine p33_lim ( a, b )

!*****************************************************************************80
!
!! P33_LIM returns the integration limits for problem 33.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p34_exact ( exact )

!*****************************************************************************80
!
!! P34_EXACT returns the exact integral for problem 34.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 1627879.0D+00 / 1500.0D+00

  return
end
subroutine p34_fun ( n, x, fx )

!*****************************************************************************80
!
!! P34_FUN evaluates the integrand for problem 34.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    ( 10 * X - 1 ) * ( 10 * X - 1.1 ) * ( 10 * X - 1.2 ) * ( 10 * X - 1.3 )
!
!  Exact Integral:
!
!    1627879 / 1500
!
!  Approximate Integral (20 digits):
!
!    1085.2526666666666666...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hermann Engels,
!    Numerical Quadrature and Cubature,
!    Academic Press, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = ( 10.0D+00 * x(1:n) - 1.0D+00 ) * ( 10.0D+00 * x(1:n) - 1.1D+00 ) &
    * ( 10.0D+00 * x(1:n) - 1.2D+00 ) * ( 10.0D+00 * x(1:n) - 1.3D+00 )

  return
end
subroutine p34_lim ( a, b )

!*****************************************************************************80
!
!! P34_LIM returns the integration limits for problem 34.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p35_exact ( exact )

!*****************************************************************************80
!
!! P35_EXACT returns the exact integral for problem 35.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 26.0D+00

  return
end
subroutine p35_fun ( n, x, fx )

!*****************************************************************************80
!
!! P35_FUN evaluates the integrand for problem 35.
!
!  Interval:
!
!    -9 <= X <= 100
!
!  Integrand:
!
!    1 / sqrt ( abs ( X ) )
!
!  Exact Integral:
!
!    26
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hermann Engels,
!    Numerical Quadrature and Cubature,
!    Academic Press, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = 1.0D+00 / sqrt ( abs ( x(i) ) )
    end if

  end do

  return
end
subroutine p35_lim ( a, b )

!*****************************************************************************80
!
!! P35_LIM returns the integration limits for problem 35.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = -9.0D+00
  b = 100.0D+00

  return
end
subroutine p36_exact ( exact )

!*****************************************************************************80
!
!! P36_EXACT returns the exact integral for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) exact

  call p36_param_get ( alpha )

  exact = 1.0D+00 / ( alpha + 1.0D+00 )**2

  return
end
subroutine p36_fun ( n, x, fx )

!*****************************************************************************80
!
!! P36_FUN evaluates the integrand for problem 36.
!
!  Discussion:
!
!    The problem has a parameter ALPHA that can be set by calling
!    P36_PARAM_SET.  It had a default value of -0.9.
!
!    The integrand has an endpoint singularity at X=0.
!
!    Suggested values of ALPHA include -0.9 through 2.6.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    X^alpha * ln ( 1 / X )
!
!  Exact Integral:
!
!    1 / ( alpha + 1 )^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) value

  call p36_param_get ( alpha )

  do i = 1, n

    if ( x(i) <= 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = x(i)**(alpha) * log ( 1.0D+00 / x(i) )
    end if

  end do

  return
end
subroutine p36_lim ( a, b )

!*****************************************************************************80
!
!! P36_LIM returns the integration limits for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p36_param ( action, name, value )

!*****************************************************************************80
!
!! P36_PARAM gets or sets the parameter values for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.
!    'get' to get the value.
!    'set' to set the value.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'alpha' is the only option.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If the action is 'get', then VALUE returns the current parameter value.
!    If ACTION is 'set', then the parameter value is set to VALUE.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = -0.9D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      value = alpha
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P36_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      alpha = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P36_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P36_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action.'
    stop
  end if

  return
end
subroutine p36_param_get ( alpha )

!*****************************************************************************80
!
!! P36_PARAM_GET returns the parameter values for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the current value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p36_param ( 'get', 'alpha', alpha )

  return
end
subroutine p36_param_set ( alpha )

!*****************************************************************************80
!
!! P36_PARAM_SET sets the parameter values for problem 36.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the new value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p36_param ( 'set', 'alpha', alpha )

  return
end
subroutine p37_exact ( exact )

!*****************************************************************************80
!
!! P37_EXACT returns the exact integral for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  call p37_param_get ( alpha )

  exact = atan ( ( 4.0D+00 - pi ) * 4.0D+00**( alpha - 1.0D+00 ) ) &
        + atan (             pi   * 4.0D+00**( alpha - 1.0D+00 ) )

  return
end
subroutine p37_fun ( n, x, fx )

!*****************************************************************************80
!
!! P37_FUN evaluates the integrand for problem 37.
!
!  Discussion:
!
!    The problem has a parameter ALPHA that can be set by calling
!    P37_PARAM_SET.  It had a default value of 5.0.
!
!    The integrand has a peak of height 4^ALPHA at X = PI/4.
!
!    Suggested values of ALPHA include 0 through 20.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    4^(-ALPHA) / ( (X-PI/4)^2 + 16^(-ALPHA) )
!
!  Exact Integral:
!
!    atan ( ( 4 - PI ) * 4^(ALPHA-1) ) + atan ( PI * 4^(ALPHA-1) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  call p37_param_get ( alpha )

  do i = 1, n

    fx(i) = 4.0D+00**( -alpha ) &
      / ( ( x(i) - 0.25D+00 * pi )**2 + 16.0D+00**(-alpha) )

  end do

  return
end
subroutine p37_lim ( a, b )

!*****************************************************************************80
!
!! P37_LIM returns the integration limits for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p37_param ( action, name, value )

!*****************************************************************************80
!
!! P37_PARAM gets or sets the parameter values for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.
!    'get' to get the value.
!    'set' to set the value.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'alpha' is the only option.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If the action is 'get', then VALUE returns the current parameter value.
!    If ACTION is 'set', then the parameter value is set to VALUE.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 5.0D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      value = alpha
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P37_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      alpha = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P37_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P37_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action.'
    stop
  end if

  return
end
subroutine p37_param_get ( alpha )

!*****************************************************************************80
!
!! P37_PARAM_GET returns the parameter values for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the current value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p37_param ( 'get', 'alpha', alpha )

  return
end
subroutine p37_param_set ( alpha )

!*****************************************************************************80
!
!! P37_PARAM_SET sets the parameter values for problem 37.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the new value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p37_param ( 'set', 'alpha', alpha )

  return
end
subroutine p38_exact ( exact )

!*****************************************************************************80
!
!! P38_EXACT returns the exact integral for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) besj0
  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  call p38_param_get ( alpha )

  x = 2.0D+00**alpha

  exact = pi * besj0 ( x )

  return
end
subroutine p38_fun ( n, x, fx )

!*****************************************************************************80
!
!! P38_FUN evaluates the integrand for problem 38.
!
!  Discussion:
!
!    The problem has a parameter ALPHA that can be set by calling
!    P38_PARAM_SET.
!
!    The integrand oscillates more strongly as ALPHA is increased.
!
!    The suggested range for ALPHA is 0 to 10.
!
!  Interval:
!
!    0 <= X <= PI
!
!  Integrand:
!
!    cos ( 2^ALPHA * sin ( x ) )
!
!  Exact Integral:
!
!    pi * J0 ( 2^ALPHA )
!
!    where J0 ( x ) is the J Bessel function of order 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) value

  call p38_param_get ( alpha )

  fx(1:n) = cos ( 2.0D+00**alpha * sin ( x(1:n) ) )

  return
end
subroutine p38_lim ( a, b )

!*****************************************************************************80
!
!! P38_LIM returns the integration limits for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  a = 0.0D+00
  b = pi

  return
end
subroutine p38_param ( action, name, value )

!*****************************************************************************80
!
!! P38_PARAM gets or sets the parameter values for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.
!    'get' to get the value.
!    'set' to set the value.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'alpha' is the only option.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If the action is 'get', then VALUE returns the current parameter value.
!    If ACTION is 'set', then the parameter value is set to VALUE.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 3.0D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      value = alpha
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P38_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      alpha = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P38_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P38_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action.'
    stop
  end if

  return
end
subroutine p38_param_get ( alpha )

!*****************************************************************************80
!
!! P38_PARAM_GET returns the parameter values for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the current value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p38_param ( 'get', 'alpha', alpha )

  return
end
subroutine p38_param_set ( alpha )

!*****************************************************************************80
!
!! P38_PARAM_SET sets the parameter values for problem 38.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the new value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p38_param ( 'set', 'alpha', alpha )

  return
end
subroutine p39_exact ( exact )

!*****************************************************************************80
!
!! P39_EXACT returns the exact integral for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) exact

  call p39_param_get ( alpha )

  exact = ( ( 2.0D+00 / 3.0D+00 )**( alpha + 1.0D+00 ) &
          + ( 1.0D+00 / 3.0D+00 )**( alpha + 1.0D+00 ) ) / ( alpha + 1.0D+00 )

  return
end
subroutine p39_fun ( n, x, fx )

!*****************************************************************************80
!
!! P39_FUN evaluates the integrand for problem 39.
!
!  Discussion:
!
!    The problem has a parameter ALPHA that can be set by calling
!    P39_PARAM_SET.
!
!    The integrand has a singularity at an internal point ( x = 1/3 )
!    with small binary period.
!
!    The suggested range for ALPHA is -0.8 through 2.1.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    ( abs ( x - 1/3 ) )^alpha
!
!  Exact Integral:
!
!    ( (2/3)^(alpha+1) + (1/3)^(alpha+1) ) / ( alpha + 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) alpha

  call p39_param_get ( alpha )

  do i = 1, n

    if ( x(i) - 1.0D+00 / 3.0D+00 == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = ( abs ( x(i) - 1.0D+00 / 3.0D+00 ) )**alpha
    end if

  end do

  return
end
subroutine p39_lim ( a, b )

!*****************************************************************************80
!
!! P39_LIM returns the integration limits for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p39_param ( action, name, value )

!*****************************************************************************80
!
!! P39_PARAM gets or sets the parameter values for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.
!    'get' to get the value.
!    'set' to set the value.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'alpha' is the only option.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If the action is 'get', then VALUE returns the current parameter value.
!    If ACTION is 'set', then the parameter value is set to VALUE.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = -0.5D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      value = alpha
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P39_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      alpha = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P39_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P39_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action.'
    stop
  end if

  return
end
subroutine p39_param_get ( alpha )

!*****************************************************************************80
!
!! P39_PARAM_GET returns the parameter values for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the current value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p39_param ( 'get', 'alpha', alpha )

  return
end
subroutine p39_param_set ( alpha )

!*****************************************************************************80
!
!! P39_PARAM_SET sets the parameter values for problem 39.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the new value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p39_param ( 'set', 'alpha', alpha )

  return
end
subroutine p40_exact ( exact )

!*****************************************************************************80
!
!! P40_EXACT returns the exact integral for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  call p40_param_get ( alpha )

  exact = ( ( 1.0 - 0.25D+00 * pi )**( alpha + 1.0D+00 ) &
          + (     + 0.25D+00 * pi )**( alpha + 1.0D+00 ) ) &
          / ( alpha + 1.0D+00 )

  return
end
subroutine p40_fun ( n, x, fx )

!*****************************************************************************80
!
!! P40_FUN evaluates the integrand for problem 40.
!
!  Discussion:
!
!    The problem has a parameter ALPHA that can be set by calling
!    P40_PARAM_SET.
!
!    The integrand has a singularity at an internal point ( x = pi/4 ).
!
!    The suggested range for ALPHA is -0.8 through 2.1.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    ( abs ( x - pi/4 ) )^alpha
!
!  Exact Integral:
!
!    ( (1-pi/4)^(alpha+1) + (pi/4)^(alpha+1) ) / ( alpha + 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  call p40_param_get ( alpha )

  fx(1:n) = ( abs ( x(1:n) - 0.25D+00 * pi ) )**alpha

  return
end
subroutine p40_lim ( a, b )

!*****************************************************************************80
!
!! P40_LIM returns the integration limits for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p40_param ( action, name, value )

!*****************************************************************************80
!
!! P40_PARAM gets or sets the parameter values for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.
!    'get' to get the value.
!    'set' to set the value.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'alpha' is the only option.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If the action is 'get', then VALUE returns the current parameter value.
!    If ACTION is 'set', then the parameter value is set to VALUE.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = -0.5D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      value = alpha
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P40_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      alpha = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P40_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P40_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action.'
    stop
  end if

  return
end
subroutine p40_param_get ( alpha )

!*****************************************************************************80
!
!! P40_PARAM_GET returns the parameter values for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the current value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p40_param ( 'get', 'alpha', alpha )

  return
end
subroutine p40_param_set ( alpha )

!*****************************************************************************80
!
!! P40_PARAM_SET sets the parameter values for problem 40.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the new value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p40_param ( 'set', 'alpha', alpha )

  return
end
subroutine p41_exact ( exact )

!*****************************************************************************80
!
!! P41_EXACT returns the exact integral for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  call p41_param_get ( alpha )

  exact = pi / sqrt ( ( 1.0D+00 + 2.0D+00**(-alpha) )**2 - 1.0D+00 )

  return
end
subroutine p41_fun ( n, x, fx )

!*****************************************************************************80
!
!! P41_FUN evaluates the integrand for problem 41.
!
!  Discussion:
!
!    The problem has a parameter ALPHA that can be set by calling
!    P41_PARAM_SET.
!
!    The integrand has a singularity at both endpoints, whose
!    severity increases with ALPHA.
!
!    The suggested range for ALPHA is 1 through 20.
!
!  Interval:
!
!    -1 <= X <= 1
!
!  Integrand:
!
!    1 / ( sqrt ( 1 - x^2 ) * ( x + 1 + 2^(-alpha) ) )
!
!  Exact Integral:
!
!    pi / sqrt ( ( 1 + 2^(-alpha) ) - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) alpha

  call p41_param_get ( alpha )

  do i = 1, n

    if ( 1.0D+00 - x(i)**2 == 0.0D+00 .or. &
         x(i) + 1.0D+00 + 0.5D+00**alpha == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = 1.0D+00 / &
        ( sqrt ( 1.0D+00 - x(i)**2 ) * ( x(i) + 1.0D+00 + 0.5D+00**alpha ) )
    end if

  end do

  return
end
subroutine p41_lim ( a, b )

!*****************************************************************************80
!
!! P41_LIM returns the integration limits for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = -1.0D+00
  b = 1.0D+00

  return
end
subroutine p41_param ( action, name, value )

!*****************************************************************************80
!
!! P41_PARAM gets or sets the parameter values for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.
!    'get' to get the value.
!    'set' to set the value.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'alpha' is the only option.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If the action is 'get', then VALUE returns the current parameter value.
!    If ACTION is 'set', then the parameter value is set to VALUE.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 3.0D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      value = alpha
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P41_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      alpha = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P41_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P41_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action.'
    stop
  end if

  return
end
subroutine p41_param_get ( alpha )

!*****************************************************************************80
!
!! P41_PARAM_GET returns the parameter values for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the current value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p41_param ( 'get', 'alpha', alpha )

  return
end
subroutine p41_param_set ( alpha )

!*****************************************************************************80
!
!! P41_PARAM_SET sets the parameter values for problem 41.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the new value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p41_param ( 'set', 'alpha', alpha )

  return
end
subroutine p42_exact ( exact )

!*****************************************************************************80
!
!! P42_EXACT returns the exact integral for problem 42.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) exact
  real ( kind = 8 ) r8_gamma

  call p42_param_get ( alpha )

  exact = 2.0D+00**( alpha - 2.0D+00 ) * ( r8_gamma ( alpha / 2.0D+00 ) )**2 &
    / r8_gamma ( alpha )

  return
end
subroutine p42_fun ( n, x, fx )

!*****************************************************************************80
!
!! P42_FUN evaluates the integrand for problem 42.
!
!  Discussion:
!
!    The problem has a parameter ALPHA that can be set by calling
!    P42_PARAM_SET.
!
!    The integrand has a singularity at X = 0 if ALPHA < 1.
!
!    The suggested range for ALPHA is 0.1 through 2.
!
!  Interval:
!
!    0 <= X <= pi/2
!
!  Integrand:
!
!    ( sin(x) )^( alpha - 1 )
!
!  Exact Integral:
!
!    2^( alpha - 2 ) * ( Gamma(alpha/2) )^2 / Gamma(alpha)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) base

  call p42_param_get ( alpha )

  do i = 1, n

    base = sin ( x(i) )

    if ( base == 0.0D+00 ) then

      if ( 1.0D+00 < alpha ) then
        fx(i) = 0.0D+00
      else if ( alpha == 1.0D+00 ) then
        fx(i) = 1.0D+00
      else
        fx(i) = 0.0D+00
      end if

    else

      fx(i) = base**( alpha - 1.0D+00 )

    end if

  end do

  return
end
subroutine p42_lim ( a, b )

!*****************************************************************************80
!
!! P42_LIM returns the integration limits for problem 42.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  a = 0.0D+00
  b = pi / 2.0D+00

  return
end
subroutine p42_param ( action, name, value )

!*****************************************************************************80
!
!! P42_PARAM gets or sets the parameter values for problem 42.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.
!    'get' to get the value.
!    'set' to set the value.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'alpha' is the only option.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If the action is 'get', then VALUE returns the current parameter value.
!    If ACTION is 'set', then the parameter value is set to VALUE.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 0.3D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      value = alpha
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P42_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      alpha = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P42_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P42_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action.'
    stop
  end if

  return
end
subroutine p42_param_get ( alpha )

!*****************************************************************************80
!
!! P42_PARAM_GET returns the parameter values for problem 42.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the current value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p42_param ( 'get', 'alpha', alpha )

  return
end
subroutine p42_param_set ( alpha )

!*****************************************************************************80
!
!! P42_PARAM_SET sets the parameter values for problem 42.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the new value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p42_param ( 'set', 'alpha', alpha )

  return
end
subroutine p43_exact ( exact )

!*****************************************************************************80
!
!! P43_EXACT returns the exact integral for problem 43.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) exact
  real ( kind = 8 ) r8_gamma

  call p43_param_get ( alpha )

  exact = r8_gamma ( alpha )

  return
end
subroutine p43_fun ( n, x, fx )

!*****************************************************************************80
!
!! P43_FUN evaluates the integrand for problem 43.
!
!  Discussion:
!
!    The problem has a parameter ALPHA that can be set by calling
!    P43_PARAM_SET.
!
!    The suggested parameter range is 0.1 <= ALPHA <= 2.0.
!
!    The integrand has an algebraic endpoint singularity at X = 1
!    times a singular factor.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    ( ln ( 1 / x ) )^( alpha - 1 )
!
!  Exact Integral:
!
!    Gamma(alpha)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 84.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) alpha

  call p43_param_get ( alpha )

  do i = 1, n

    if ( x(i) <= 0.0D+00 ) then

      fx(i) = 0.0D+00

    else if ( x(i) == 0.0D+00 ) then

      if ( alpha - 1.0D+00 < 0.0D+00 ) then
        fx(i) = 0.0D+00
      else if ( alpha - 1.0D+00 == 0.0D+00 ) then
        fx(i) = 1.0D+00
      else
        fx(i) = 0.0D+00
      end if

    else if ( x(i) == 1.0D+00 ) then

      if ( alpha - 1.0D+00 < 0.0D+00 ) then
        fx(i) = 0.0D+00
      else if ( alpha - 1.0D+00 == 0.0D+00 ) then
        fx(i) = 1.0D+00
      else
        fx(i) = 0.0D+00
      end if

    else

      fx(i) = ( log ( 1.0D+00 / x(i) ) )**( alpha - 1.0D+00 )

    end if

  end do

  return
end
subroutine p43_lim ( a, b )

!*****************************************************************************80
!
!! P43_LIM returns the integration limits for problem 43.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p43_param ( action, name, value )

!*****************************************************************************80
!
!! P43_PARAM gets or sets the parameter values for problem 43.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.
!    'get' to get the value.
!    'set' to set the value.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'alpha' is the only option.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If the action is 'get', then VALUE returns the current parameter value.
!    If ACTION is 'set', then the parameter value is set to VALUE.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 0.3D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      value = alpha
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P43_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      alpha = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P43_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P43_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action.'
    stop
  end if

  return
end
subroutine p43_param_get ( alpha )

!*****************************************************************************80
!
!! P43_PARAM_GET returns the parameter values for problem 43.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the current value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p43_param ( 'get', 'alpha', alpha )

  return
end
subroutine p43_param_set ( alpha )

!*****************************************************************************80
!
!! P43_PARAM_SET sets the parameter values for problem 43.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the new value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p43_param ( 'set', 'alpha', alpha )

  return
end
subroutine p44_exact ( exact )

!*****************************************************************************80
!
!! P44_EXACT returns the exact integral for problem 44.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) exact

  call p44_param_get ( alpha )

  exact = ( 20.0D+00 * sin ( 2.0D+00**alpha ) &
    - 2.0D+00**alpha * cos ( 2.0D+00**alpha ) &
    + 2.0D+00**alpha * exp ( -20.0D+00 ) ) / ( 400.0D+00 + 4.0D+00**alpha )

  return
end
subroutine p44_fun ( n, x, fx )

!*****************************************************************************80
!
!! P44_FUN evaluates the integrand for problem 44.
!
!  Discussion:
!
!    The problem has a parameter ALPHA that can be set by calling
!    P44_PARAM_SET.
!
!    The suggested parameter range is 0.0 <= ALPHA <= 9.0.
!
!    As ALPHA increases, the integrand becomes more oscillatory.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    exp ( 20 * ( x - 1 ) ) * sin ( 2^alpha * x )
!
!  Exact Integral:
!
!    ( 20 sin ( 2^alpha ) - 2^alpha cos ( 2^alpha )
!    + 2^alpha exp ( -20 ) ) / ( 400 + 4^alpha )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 84.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) alpha

  call p44_param_get ( alpha )

  fx(1:n) = exp ( 20.0D+00 * ( x(1:n) - 1.0D+00 ) ) &
    * sin ( 2.0D+00**alpha * x(1:n) )

  return
end
subroutine p44_lim ( a, b )

!*****************************************************************************80
!
!! P44_LIM returns the integration limits for problem 44.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p44_param ( action, name, value )

!*****************************************************************************80
!
!! P44_PARAM gets or sets the parameter values for problem 44.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.
!    'get' to get the value.
!    'set' to set the value.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'alpha' is the only option.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If the action is 'get', then VALUE returns the current parameter value.
!    If ACTION is 'set', then the parameter value is set to VALUE.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 2.0D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      value = alpha
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P44_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      alpha = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P44_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P44_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action.'
    stop
  end if

  return
end
subroutine p44_param_get ( alpha )

!*****************************************************************************80
!
!! P44_PARAM_GET returns the parameter values for problem 44.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the current value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p44_param ( 'get', 'alpha', alpha )

  return
end
subroutine p44_param_set ( alpha )

!*****************************************************************************80
!
!! P44_PARAM_SET sets the parameter values for problem 44.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the new value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p44_param ( 'set', 'alpha', alpha )

  return
end
subroutine p45_exact ( exact )

!*****************************************************************************80
!
!! P45_EXACT returns the exact integral for problem 45.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) besj0
  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  call p45_param_get ( alpha )

  exact = pi * cos ( 2.0D+00**( alpha - 1.0D+00 ) ) &
    * besj0 ( 2.0D+00**( alpha - 1.0D+00 ) )

  return
end
subroutine p45_fun ( n, x, fx )

!*****************************************************************************80
!
!! P45_FUN evaluates the integrand for problem 45.
!
!  Discussion:
!
!    The problem has a parameter ALPHA that can be set by calling
!    P45_PARAM_SET.
!
!    The suggested parameter range is 0.0 <= ALPHA <= 8.0.
!
!    The function is singular at 0 and 1.
!
!  Interval:
!
!    0 <= X <= 1
!
!  Integrand:
!
!    cos ( 2^alpha * x ) / sqrt ( x * ( 1 - x ) )
!
!  Exact Integral:
!
!    pi * cos ( 2^(alpha-1) ) * J0 ( 2^(alpha-1) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 84.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) alpha

  call p45_param_get ( alpha )

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else if ( x(i) == 1.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = cos ( 2.0D+00**alpha * x(i) ) / sqrt ( x(i) * ( 1.0D+00 - x(i) ) )
    end if

  end do

  return
end
subroutine p45_lim ( a, b )

!*****************************************************************************80
!
!! P45_LIM returns the integration limits for problem 45.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p45_param ( action, name, value )

!*****************************************************************************80
!
!! P45_PARAM gets or sets the parameter values for problem 45.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.
!    'get' to get the value.
!    'set' to set the value.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'alpha' is the only option.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If the action is 'get', then VALUE returns the current parameter value.
!    If ACTION is 'set', then the parameter value is set to VALUE.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 2.0D+00
  character ( len = * ) name
  logical s_eqi
  real ( kind = 8 ) value

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      value = alpha
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P45_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      alpha = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P45_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P45_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action.'
    stop
  end if

  return
end
subroutine p45_param_get ( alpha )

!*****************************************************************************80
!
!! P45_PARAM_GET returns the parameter values for problem 45.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the current value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p45_param ( 'get', 'alpha', alpha )

  return
end
subroutine p45_param_set ( alpha )

!*****************************************************************************80
!
!! P45_PARAM_SET sets the parameter values for problem 45.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the new value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p45_param ( 'set', 'alpha', alpha )

  return
end
subroutine p46_exact ( exact )

!*****************************************************************************80
!
!! P46_EXACT returns the exact integral for problem 46.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 6.0690909595647754101D+00

  return
end
subroutine p46_fun ( n, x, fx )

!*****************************************************************************80
!
!! P46_FUN evaluates the integrand for problem 46.
!
!  Discussion:
!
!    The problem has a parameter ALPHA that can be set by calling
!    P63_PARAM_SET.
!
!    The integrand is the radius of an ellipse as a function of angle.
!
!    The integral represents the arc length of the ellipse.
!
!    The suggested parameter range is 0.0 <= ALPHA < 1.0.  ALPHA is
!    the eccentricity of the ellipse.
!
!  Interval:
!
!    0 <= theta <= 2 pi
!
!  Integrand:
!
!    r(theta) = ( 1 - alpha^2 ) / ( 1 - alpha * cos ( theta ) )
!
!  Exact Integral:
!
!    When alpha = sin ( pi / 12 ), then
!
!      6.0690909595647754101
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Crandall,
!    Projects in Scientific Computing,
!    Springer, 2000, page 47.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) alpha

  call p46_param_get ( alpha )

  fx(1:n) = ( 1.0D+00 - alpha**2 ) / ( 1.0D+00 - alpha * cos ( x(1:n) ) )

  return
end
subroutine p46_lim ( a, b )

!*****************************************************************************80
!
!! P46_LIM returns the integration limits for problem 46.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  a = 0.0D+00
  b = 2.0D+00 * pi

  return
end
subroutine p46_param ( action, name, value )

!*****************************************************************************80
!
!! P46_PARAM gets or sets the parameter values for problem 46.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION.
!    'get' to get the value.
!    'set' to set the value.
!
!    Input, character ( len = * ) NAME, the name of the parameter.
!    'alpha' is the only option.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If the action is 'get', then VALUE returns the current parameter value.
!    If ACTION is 'set', then the parameter value is set to VALUE.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 0.0D+00
  character ( len = * ) name
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  logical s_eqi
  logical, save :: set = .false.
  real ( kind = 8 ) value

  if ( .not. set ) then
    alpha = sin ( pi / 12.0D+00 )
    set = .true.
  end if

  if ( s_eqi ( action, 'get' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      value = alpha
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P46_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else if ( s_eqi ( action, 'set' ) ) then
    if ( s_eqi ( name, 'alpha' ) ) then
      alpha = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P46_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name.'
      stop
    end if
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P46_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action.'
    stop
  end if

  return
end
subroutine p46_param_get ( alpha )

!*****************************************************************************80
!
!! P46_PARAM_GET returns the parameter values for problem 46.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the current value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p46_param ( 'get', 'alpha', alpha )

  return
end
subroutine p46_param_set ( alpha )

!*****************************************************************************80
!
!! P46_PARAM_SET sets the parameter values for problem 46.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the new value of the parameter.
!
  implicit none

  real ( kind = 8 ) alpha

  call p46_param ( 'set', 'alpha', alpha )

  return
end
subroutine p47_exact ( exact )

!*****************************************************************************80
!
!! P47_EXACT returns the exact integral for problem 47.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = - 4.0D+00 / 9.0D+00

  return
end
subroutine p47_fun ( n, x, fx )

!*****************************************************************************80
!
!! P47_FUN evaluates the integrand for problem 47.
!
!  Discussion:
!
!    The function is singular at the left endpoint.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    sqrt ( x ) * ln ( x )
!
!  Exact Integral:
!
!    -4/9 = -0.4444...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 101.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = sqrt ( x(i) ) * log ( x(i) )
    end if

  end do

  return
end
subroutine p47_lim ( a, b )

!*****************************************************************************80
!
!! P47_LIM returns the integration limits for problem 47.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p48_exact ( exact )

!*****************************************************************************80
!
!! P48_EXACT returns the exact integral for problem 48.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = -4.0D+00

  return
end
subroutine p48_fun ( n, x, fx )

!*****************************************************************************80
!
!! P48_FUN evaluates the integrand for problem 48.
!
!  Discussion:
!
!    The function is singular at the left endpoint.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    ln ( x ) / sqrt ( x )
!
!  Exact Integral:
!
!    -4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 103.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = log ( x(i) ) / sqrt ( x(i) )
    end if

  end do

  return
end
subroutine p48_lim ( a, b )

!*****************************************************************************80
!
!! P48_LIM returns the integration limits for problem 48.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p49_exact ( exact )

!*****************************************************************************80
!
!! P49_EXACT returns the exact integral for problem 49.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 61.0D+00 * log ( 2.0D+00 ) &
    + 77.0D+00 * log ( 7.0D+00 ) / 4.0D+00 - 27.0D+00

  return
end
subroutine p49_fun ( n, x, fx )

!*****************************************************************************80
!
!! P49_FUN evaluates the integrand for problem 49.
!
!  Discussion:
!
!    The function is singular at two internal points, 1 and sqrt(2).
!
!  Interval:
!
!    0 <= x <= 3
!
!  Integrand:
!
!    x^3 * log ( abs ( ( x^2 - 1 ) * ( x^2 - 2 ) ) )
!
!  Exact Integral:
!
!    61 log ( 2 ) + (77/4) log ( 7 ) - 27 = 52.7408...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 104.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( ( x(i)**2 - 1.0D+00 ) * ( x(i)**2 - 2.0D+00 ) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = x(i)**3 * log ( abs ( ( x(i)**2 - 1.0D+00 ) &
        * ( x(i)**2 - 2.0D+00 ) ) )
    end if

  end do

  return
end
subroutine p49_lim ( a, b )

!*****************************************************************************80
!
!! P49_LIM returns the integration limits for problem 49.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 3.0D+00

  return
end
subroutine p50_exact ( exact )

!*****************************************************************************80
!
!! P50_EXACT returns the exact integral for problem 50.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  real ( kind = 8 ) euler_constant
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_ci
  real ( kind = 8 ) t

  t = 10.0D+00 * pi

  exact = ( - euler_constant ( ) - log ( t ) + r8_ci ( t ) ) / t

  return
end
subroutine p50_fun ( n, x, fx )

!*****************************************************************************80
!
!! P50_FUN evaluates the integrand for problem 50.
!
!  Discussion:
!
!    The function has a removable singularity at x = 0.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    log ( x ) * sin ( 10 * pi * x )
!
!  Exact Integral:
!
!    ( - gamma - log ( 10 * pi ) + Ci ( 10 * pi ) ) / 10 * pi = -0.1281316...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 106.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = log ( x(i) ) * sin ( 10.0D+00 * pi * x(i) )
    end if

  end do

  return
end
subroutine p50_lim ( a, b )

!*****************************************************************************80
!
!! P50_LIM returns the integration limits for problem 50.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p51_exact ( exact )

!*****************************************************************************80
!
!! P51_EXACT returns the exact integral for problem 51.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_ci
  real ( kind = 8 ) r8_si

  exact = - ( r8_ci ( 1.0D+00 ) * sin ( 1.0D+00 ) + &
    ( 0.5D+00 * pi - r8_si ( 1.0D+00 ) ) * cos ( 1.0D+00 ) ) / pi

  return
end
subroutine p51_fun ( n, x, fx )

!*****************************************************************************80
!
!! P51_FUN evaluates the integrand for problem 51.
!
!  Discussion:
!
!    The function has a removable singularity at x = 0.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    ln ( x ) / ( 1 + ( ln(x) )^2 )^2
!
!  Exact Integral:
!
!    - ( ci(1) * sin(1) + ( pi/2 - si(1) ) * cos(1) ) / pi = - 0.1892752...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 108.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = log ( x(i) ) / ( 1.0D+00 + ( log ( x(i) ) )**2 )**2
    end if

  end do

  return
end
subroutine p51_lim ( a, b )

!*****************************************************************************80
!
!! P51_LIM returns the integration limits for problem 51.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p52_exact ( exact )

!*****************************************************************************80
!
!! P52_EXACT returns the exact integral for problem 52.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = log ( 125.0D+00 / 631.0D+00 ) / 18.0D+00

  return
end
subroutine p52_fun ( n, x, fx )

!*****************************************************************************80
!
!! P52_FUN evaluates the integrand for problem 52.
!
!  Discussion:
!
!    The function has a singularity at x = 0.
!
!  Interval:
!
!    -1 <= x <= 5
!
!  Integrand:
!
!    1 / ( x * ( 5 * x^3 + 6 ) )
!
!  Exact Integral:
!
!    ln ( 125 / 631 ) / 18 = -0.08994401...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 109.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    if ( x(i) == 0.0D+00 ) then
      fx(i) = 0.0D+00
    else
      fx(i) = 1.0D+00 / ( x(i) * ( 5.0D+00 * x(i)**3 + 6.0D+00 ) )
    end if

  end do

  return
end
subroutine p52_lim ( a, b )

!*****************************************************************************80
!
!! P52_LIM returns the integration limits for problem 52.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = -1.0D+00
  b = 5.0D+00

  return
end
subroutine p53_exact ( exact )

!*****************************************************************************80
!
!! P53_EXACT returns the exact integral for problem 53.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  exact = 0.5D+00 * pi - atan ( 1.0D+00 / sqrt ( 2.0D+00 ) ) &
    + log ( 3.0D+00 ) / 2.0D+00

  return
end
subroutine p53_fun ( n, x, fx )

!*****************************************************************************80
!
!! P53_FUN evaluates the integrand for problem 53.
!
!  Discussion:
!
!    The integrand is singular at x = -1 + sqrt ( 3 ) = 0.732...
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    1 / sqrt ( abs ( x^2 + 2 * x - 2 ) )
!
!  Exact Integral:
!
!    pi / 2 - arctan ( 1 / sqrt ( 2 ) ) + ln ( 3 ) / 2 = 1.504622...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenga,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK: A Subroutine Package for Automatic Integration,
!    Springer, 1983, page 110.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = 1.0D+00 / sqrt ( abs ( x(1:n)**2 + 2.0D+00 * x(1:n) - 2.0D+00 ) )

  return
end
subroutine p53_lim ( a, b )

!*****************************************************************************80
!
!! P53_LIM returns the integration limits for problem 53.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p54_exact ( exact )

!*****************************************************************************80
!
!! P54_EXACT returns the exact integral for problem 54.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 2.0D+00 / sqrt ( 3.0D+00 )

  return
end
subroutine p54_fun ( n, x, fx )

!*****************************************************************************80
!
!! P54_FUN evaluates the integrand for problem 54.
!
!  Discussion:
!
!    The reference claims that this integrand is more closely approximated
!    by the trapezoid rule than by Gauss-Legendre quadrature.
!
!    Points  Trapezoid  Gauss-Legendre
!     4      1.91667    2.53883
!    12      2.1594     2.25809
!
!    However, the stated results hardly give one confidence in
!    the convergence of the trapezoid results, and I am unable to
!    confirm them, because my results for 4 points give good results
!    (about 1.14) for BOTH Trapezoid and Gauss-Legendre!
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    2 / ( 2 + sin ( 10 * PI * X ) )
!
!  Exact Integral:
!
!    2 / sqrt ( 3 ) = 1.1547...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Prem Kythe, Pratap Puri,
!    Computational Methods for Linear Integral Equations,
!    Birkhaeuser, 2002.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  fx(1:n) = 2.0D+00 / ( 2.0D+00 + sin ( 10.0D+00 * pi * x(1:n) ) )

  return
end
subroutine p54_lim ( a, b )

!*****************************************************************************80
!
!! P54_LIM returns the integration limits for problem 54.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p55_exact ( exact )

!*****************************************************************************80
!
!! P55_EXACT returns the exact integral for problem 55.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) error_f
  real ( kind = 8 ) exact
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x0

  call p55_lim ( a, b )
  call p55_param ( 'get', 'c', c )
  call p55_param ( 'get', 'x0', x0 )

  exact = sqrt ( pi ) * &
    ( error_f ( c * ( b - x0 ) ) - error_f ( c * ( a - x0 ) ) ) &
    / ( 2.0D+00 * c )

  return
end
subroutine p55_fun ( n, x, fx )

!*****************************************************************************80
!
!! P55_FUN evaluates the integrand for problem 55.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    exp ( - c^2 * ( x - x0 )^2 )
!
!  Exact Integral:
!
!    sqrt ( pi )
!    * ( erf ( c * ( b - x0 ) ) - erf ( c * ( a - x0 ) ) )
!    / ( 2 * c )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) c
  real ( kind = 8 ) x0

  call p55_param ( 'get', 'c', c )
  call p55_param ( 'get', 'x0', x0 )

  fx(1:n) = exp ( - c**2 * ( x(1:n) - x0 )**2 )

  return
end
subroutine p55_lim ( a, b )

!*****************************************************************************80
!
!! P55_LIM returns the integration limits for problem 55.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine p55_param ( action, name, value )

!*****************************************************************************80
!
!! P55_PARAM sets or gets real scalar parameters for problem 55.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION,
!    'get' to get a parameter.
!    'set' to set a parameter.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    'C' is the coefficient.
!    'X0' is the base point.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    * If ACTION = 'set', then VALUE is an input quantity, and is the
!      new value to be assigned to NAME.
!    * If ACTION = 'get', then VALUE is an output quantity, and is the
!      current value of NAME.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: c = 3.00D+00
  character ( len = * ) name
  real ( kind = 8 ) value
  real ( kind = 8 ), save :: x0 = 0.75D+00

  if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

    if ( name == 'C' .or. name == 'c' ) then
      value = c
    else if ( name == 'X0' .or. name == 'x0' ) then
      value = x0
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P55_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

    if ( name == 'C' .or. name == 'c' ) then
      c = value
    else if ( name == 'X0' .or. name == 'x0' ) then
      x0 = value
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P55_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Unrecognized name = "' // trim ( name ) // '".'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P55_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop

  end if

  return
end
subroutine p56_exact ( exact )

!*****************************************************************************80
!
!! P56_EXACT returns the estimated integral for problem 56.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 1.9922524079504000171D+00

  return
end
subroutine p56_fun ( n, x, fx )

!*****************************************************************************80
!
!! P56_FUN evaluates the integrand for problem 56.
!
!  Interval:
!
!    -1 <= x <= 1
!
!  Integrand:
!
!    1 / ( x^6 + 0.9 )
!
!  Approximate Integral (20 digits):
!
!    1.9922524079504000171...
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software,
!    edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = 1.0D+00 / ( x(1:n)**6 + 0.9D+00 )

  return
end
subroutine p56_lim ( a, b )

!*****************************************************************************80
!
!! P56_LIM returns the integration limits for problem 56.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = - 1.0D+00
  b = 1.0D+00

  return
end
subroutine p57_exact ( exact )

!*****************************************************************************80
!
!! P57_EXACT returns the exact integral for problem 57.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EXACT, the estimated value of the integral.
!
  implicit none

  real ( kind = 8 ) exact

  exact = 0.4D+00

  return
end
subroutine p57_fun ( n, x, fx )

!*****************************************************************************80
!
!! P57_FUN evaluates the integrand for problem 57.
!
!  Interval:
!
!    0 <= x <= 1
!
!  Integrand:
!
!    x^(3/2)
!
!  Antiderivative:
!
!    (2/5) * x^(5/2)
!
!  Exact Integral:
!
!    0.4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner,
!    Comparison of Numerical Quadrature Formulas,
!    in Mathematical Software, edited by John R Rice,
!    Academic Press, 1971.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the evaluation points.
!
!    Output, real ( kind = 8 ) FX(N), the integrand values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fx(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  fx(1:n) = sqrt ( x(1:n)**3 )

  return
end
subroutine p57_lim ( a, b )

!*****************************************************************************80
!
!! P57_LIM returns the integration limits for problem 57.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A, B, the limits of integration.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b

  a = 0.0D+00
  b = 1.0D+00

  return
end
subroutine r89sifg ( x, f, g )

!*****************************************************************************80
!
!! R89SIFG evaluate terms needed in an approximation to SI(X) for 4 <= |X|.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 2000
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which the terms are needed.
!    |X| must be at least 4.
!
!    Output, real ( kind = 8 ) F, G, the values of F and G at X.
!
  implicit none

  real ( kind = 8 ) eta
  real ( kind = 8 ) f
  real ( kind = 8 ), parameter, dimension ( 20 ) :: f1cs = (/ &
  -0.1191081969051363610D+00, -0.0247823144996236248D+00, &
   0.0011910281453357821D+00, -0.0000927027714388562D+00, &
   0.0000093373141568271D+00, -0.0000011058287820557D+00, &
   0.0000001464772071460D+00, -0.0000000210694496288D+00, &
   0.0000000032293492367D+00, -0.0000000005206529618D+00, &
   0.0000000000874878885D+00, -0.0000000000152176187D+00, &
   0.0000000000027257192D+00, -0.0000000000005007053D+00, &
   0.0000000000000940241D+00, -0.0000000000000180014D+00, &
   0.0000000000000035063D+00, -0.0000000000000006935D+00, &
   0.0000000000000001391D+00, -0.0000000000000000282D+00 /)
  real ( kind = 8 ) f2cs(29)
  real ( kind = 8 ) g
  real ( kind = 8 ) g1cs(21)
  real ( kind = 8 ) g2cs(34)
  integer ( kind = 4 ), save :: nf1 = 0
  integer ( kind = 4 ), save :: nf2 = 0
  integer ( kind = 4 ), save :: ng1 = 0
  integer ( kind = 4 ), save :: ng2 = 0
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xbig = 0.0D+00
  real ( kind = 8 ), save :: xbnd = 0.0D+00
  real ( kind = 8 ), save :: xmaxf = 0.0D+00
  real ( kind = 8 ), save :: xmaxg = 0.0D+00
!
! series for f2   on the interval  0.00000D+00 to  2.00000D-02
!                                    with weighted error   4.32D-17
!                                     log weighted error  16.36
!                           significant figures required  14.75
!                                decimal places required  17.10
!
  data f2cs(  1) /    -0.0348409253897013234D+00/
  data f2cs(  2) /    -0.0166842205677959686D+00/
  data f2cs(  3) /     0.0006752901241237738D+00/
  data f2cs(  4) /    -0.0000535066622544701D+00/
  data f2cs(  5) /     0.0000062693421779007D+00/
  data f2cs(  6) /    -0.0000009526638801991D+00/
  data f2cs(  7) /     0.0000001745629224251D+00/
  data f2cs(  8) /    -0.0000000368795403065D+00/
  data f2cs(  9) /     0.0000000087202677705D+00/
  data f2cs( 10) /    -0.0000000022601970392D+00/
  data f2cs( 11) /     0.0000000006324624977D+00/
  data f2cs( 12) /    -0.0000000001888911889D+00/
  data f2cs( 13) /     0.0000000000596774674D+00/
  data f2cs( 14) /    -0.0000000000198044313D+00/
  data f2cs( 15) /     0.0000000000068641396D+00/
  data f2cs( 16) /    -0.0000000000024731020D+00/
  data f2cs( 17) /     0.0000000000009226360D+00/
  data f2cs( 18) /    -0.0000000000003552364D+00/
  data f2cs( 19) /     0.0000000000001407606D+00/
  data f2cs( 20) /    -0.0000000000000572623D+00/
  data f2cs( 21) /     0.0000000000000238654D+00/
  data f2cs( 22) /    -0.0000000000000101714D+00/
  data f2cs( 23) /     0.0000000000000044259D+00/
  data f2cs( 24) /    -0.0000000000000019634D+00/
  data f2cs( 25) /     0.0000000000000008868D+00/
  data f2cs( 26) /    -0.0000000000000004074D+00/
  data f2cs( 27) /     0.0000000000000001901D+00/
  data f2cs( 28) /    -0.0000000000000000900D+00/
  data f2cs( 29) /     0.0000000000000000432D+00/
!
! series for g1   on the interval  2.00000D-02 to  6.25000D-02
!                                    with weighted error   5.48D-17
!                                     log weighted error  16.26
!                           significant figures required  15.47
!                                decimal places required  16.92
!
  data g1cs(  1) /    -0.3040578798253495954D+00/
  data g1cs(  2) /    -0.0566890984597120588D+00/
  data g1cs(  3) /     0.0039046158173275644D+00/
  data g1cs(  4) /    -0.0003746075959202261D+00/
  data g1cs(  5) /     0.0000435431556559844D+00/
  data g1cs(  6) /    -0.0000057417294453025D+00/
  data g1cs(  7) /     0.0000008282552104503D+00/
  data g1cs(  8) /    -0.0000001278245892595D+00/
  data g1cs(  9) /     0.0000000207978352949D+00/
  data g1cs( 10) /    -0.0000000035313205922D+00/
  data g1cs( 11) /     0.0000000006210824236D+00/
  data g1cs( 12) /    -0.0000000001125215474D+00/
  data g1cs( 13) /     0.0000000000209088918D+00/
  data g1cs( 14) /    -0.0000000000039715832D+00/
  data g1cs( 15) /     0.0000000000007690431D+00/
  data g1cs( 16) /    -0.0000000000001514697D+00/
  data g1cs( 17) /     0.0000000000000302892D+00/
  data g1cs( 18) /    -0.0000000000000061400D+00/
  data g1cs( 19) /     0.0000000000000012601D+00/
  data g1cs( 20) /    -0.0000000000000002615D+00/
  data g1cs( 21) /     0.0000000000000000548D+00/
!
! series for g2   on the interval  0.00000D+00 to  2.00000D-02
!                                    with weighted error   5.01D-17
!                                     log weighted error  16.30
!                           significant figures required  15.12
!                                decimal places required  17.07
!
  data g2cs(  1) /    -0.0967329367532432218D+00/
  data g2cs(  2) /    -0.0452077907957459871D+00/
  data g2cs(  3) /     0.0028190005352706523D+00/
  data g2cs(  4) /    -0.0002899167740759160D+00/
  data g2cs(  5) /     0.0000407444664601121D+00/
  data g2cs(  6) /    -0.0000071056382192354D+00/
  data g2cs(  7) /     0.0000014534723163019D+00/
  data g2cs(  8) /    -0.0000003364116512503D+00/
  data g2cs(  9) /     0.0000000859774367886D+00/
  data g2cs( 10) /    -0.0000000238437656302D+00/
  data g2cs( 11) /     0.0000000070831906340D+00/
  data g2cs( 12) /    -0.0000000022318068154D+00/
  data g2cs( 13) /     0.0000000007401087359D+00/
  data g2cs( 14) /    -0.0000000002567171162D+00/
  data g2cs( 15) /     0.0000000000926707021D+00/
  data g2cs( 16) /    -0.0000000000346693311D+00/
  data g2cs( 17) /     0.0000000000133950573D+00/
  data g2cs( 18) /    -0.0000000000053290754D+00/
  data g2cs( 19) /     0.0000000000021775312D+00/
  data g2cs( 20) /    -0.0000000000009118621D+00/
  data g2cs( 21) /     0.0000000000003905864D+00/
  data g2cs( 22) /    -0.0000000000001708459D+00/
  data g2cs( 23) /     0.0000000000000762015D+00/
  data g2cs( 24) /    -0.0000000000000346151D+00/
  data g2cs( 25) /     0.0000000000000159996D+00/
  data g2cs( 26) /    -0.0000000000000075213D+00/
  data g2cs( 27) /     0.0000000000000035970D+00/
  data g2cs( 28) /    -0.0000000000000017530D+00/
  data g2cs( 29) /     0.0000000000000008738D+00/
  data g2cs( 30) /    -0.0000000000000004487D+00/
  data g2cs( 31) /     0.0000000000000002397D+00/
  data g2cs( 32) /    -0.0000000000000001347D+00/
  data g2cs( 33) /     0.0000000000000000801D+00/
  data g2cs( 34) /    -0.0000000000000000501D+00/

  if ( nf1 == 0 ) then

    eta = 0.1D+00 * epsilon ( eta )
    call inits ( f1cs, 20, eta, nf1 )
    call inits ( f2cs, 29, eta, nf2 )
    call inits ( g1cs, 21, eta, ng1 )
    call inits ( g2cs, 34, eta, ng2 )

    xbig = sqrt ( 1.0D+00 / epsilon ( xbig ) )
    xmaxf = exp ( min ( - log ( tiny ( xmaxf ) ), &
      log ( r8_huge ( ) ) ) - 0.01D+00 )
    xmaxg = 1.0D+00 / sqrt( tiny ( xmaxg ) )
    xbnd = sqrt ( 50.0D+00 )

  end if

  if ( x < 4.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R89SIFG - Fatal error!'
    write ( *, '(a)' ) '  Approximation invalid for X < 4.'
    stop
  else if ( x <= xbnd ) then
    call csevl ( ( 1.0D+00 / x**2 - 0.04125D+00 ) / 0.02125D+00, f1cs, nf1, &
      value )
    f = ( 1.0D+00 + value ) / x
    call csevl ( ( 1.0D+00 / x**2 - 0.04125D+00 ) / 0.02125D+00, g1cs, ng1, &
      value )
    g = ( 1.0D+00 + value ) / x**2
  else if ( x <= xbig ) then
    call csevl ( 100.0D+00 / x**2 - 1.0D+00, f2cs, nf2, value )
    f = ( 1.0D+00 + value ) / x
    call csevl ( 100.0D+00 / x**2 - 1.0D+00, g2cs, ng2, value )
    g = ( 1.0D+00 + value ) / x**2
  else

    if ( x < xmaxf ) then
      f = 1.0D+00 / x
    else
      f = 0.0D+00
    end if

    if ( x < xmaxg ) then
      g = 1.0D+00 / x**2
    else
      g = 0.0D+00
    end if

  end if

  return
end
subroutine r89upak ( x, y, n )

!*****************************************************************************80
!
!! R89UPAK rewrites a real value X as Y * 2**N where 0.5 <= |Y| < 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 November 2000
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be decomposed.
!
!    Output, real ( kind = 8 ) Y, the mantissa of the decomposition.
!
!    Output, integer ( kind = 4 ) N, the exponent of the decomposition.
!
  implicit none

  real ( kind = 8 ) absx
  integer ( kind = 4 ) n
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  absx = abs ( x )
  n = 0
  y = 0.0D+00

  if ( x == 0.0D+00 ) then
    return
  end if

  do while ( absx < 0.5D+00 )
    n = n - 1
    absx = absx * 2.0D+00
  end do

  do while ( 1.0D+00 <= absx )
    n = n + 1
    absx = absx * 0.5D+00
  end do

  y = sign ( absx, x )

  return
end
function r8_ci ( x )

!*****************************************************************************80
!
!! R8_CI computes an approximation to the value of the cosine integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2002
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullertony.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the cosine integral.
!
!    Output, real ( kind = 8 ) R8_CI, the value of the cosine integral.
!
  implicit none

  integer ( kind = 4 ), parameter :: ncics = 13

  real ( kind = 8 ), parameter, dimension ( ncics ) :: cics = (/ &
  -0.34004281856055363156D+00, -1.03302166401177456807D+00, &
   0.19388222659917082877D+00, -0.01918260436019865894D+00, &
   0.00110789252584784967D+00, -0.00004157234558247209D+00, &
   0.00000109278524300229D+00, -0.00000002123285954183D+00, &
   0.00000000031733482164D+00, -0.00000000000376141548D+00, &
   0.00000000000003622653D+00, -0.00000000000000028912D+00, &
   0.00000000000000000194D+00 /)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ), save :: nci = 0
  real ( kind = 8 ) r8_ci
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( nci == 0 ) then
    call inits ( cics, ncics, 0.1D+00 * epsilon ( x ), nci )
  end if

  if ( x <= 0.0D+00 ) then

    r8_ci = - r8_huge ( )

  else if ( x <= sqrt ( epsilon ( x ) ) ) then

    y = -1.0D+00
    call csevl ( y, cics, nci, value )
    r8_ci = log ( x ) - 0.5D+00 + value

  else if ( x <= 4.0D+00 ) then

    y = ( x**2 - 8.0D+00 ) * 0.125D+00
    call csevl ( y, cics, nci, value )
    r8_ci = log ( x ) - 0.5D+00 + value

  else

    call r89sifg ( x, f, g )
    r8_ci = f * sin ( x ) - g * cos ( x )

  end if

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA returns the value of the Gamma function at X.
!
!  Definition:
!
!    GAMMA(Z) = Integral ( 0 <= T < +oo ) T^(Z-1) EXP(-T) dT
!
!  Recursion:
!
!    GAMMA(X+1) = X*GAMMA(X)
!
!  Restrictions:
!
!    0 < X ( a software restriction).
!
!  Special values:
!
!    GAMMA(0.5) = sqrt(PI)
!
!    For N a positive integer, GAMMA(N+1) = N!, the standard factorial.
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
!    Input, real ( kind = 8 ) X, the point at which the Gamma function
!    is desired.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the Gamma function of X.
!
  implicit none

  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) x

  r8_gamma = exp ( r8_gamma_log ( x ) )

  return
end
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! R8_GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in references 1 and 2.
!    The program uses rational functions that theoretically approximate
!    log ( GAMMA(X) ) to at least 18 significant decimal digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody and Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Math. Comp.
!    Volume 21, 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
!    John Rice, Henry Thatcher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function.
!    X must be positive.
!
!    Output, real ( kind = 8 ) GAMMA_LOG, the logarithm of the Gamma function
!    of X.  If X <= 0.0, or if overflow would occur, the program returns the
!    largest representable floating point number.
!
!  Mmachine dependent constants:
!
!    BETA   - radix for the floating-point representation.
!
!    MAXEXP - the smallest positive power of BETA that overflows.
!
!    XBIG   - largest argument for which LN(GAMMA(X)) is representable
!             in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!    FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!    Approximate values for some important machines are:
!
!                              BETA      MAXEXP         XBIG
!
!    CRAY-1        (S.P.)        2        8191       9.62D+2461
!    Cyber 180/855
!      under NOS   (S.P.)        2        1070       1.72D+319
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)        2         128       4.08D+36
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)        2        1024       2.55D+305
!    IBM 3033      (D.P.)       16          63       4.29D+73
!    VAX D-Format  (D.P.)        2         127       2.05D+36
!    VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                            FRTBIG
!
!    CRAY-1        (S.P.)   3.13D+615
!    Cyber 180/855
!      under NOS   (S.P.)   6.44D+79
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)   1.42D+9
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)   2.25D+76
!    IBM 3033      (D.P.)   2.56D+18
!    VAX D-Format  (D.P.)   1.20D+9
!    VAX G-Format  (D.P.)   1.89D+76
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) corr
  real ( kind = 8 ), parameter :: d1 = - 5.772156649015328605195174D-01
  real ( kind = 8 ), parameter :: d2 =   4.227843350984671393993777D-01
  real ( kind = 8 ), parameter :: d4 =   1.791759469228055000094023D+00
  real ( kind = 8 ) eps
  real ( kind = 8 ), parameter :: frtbig = 1.42D+09
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = 8 ), parameter :: pnt68 = 0.6796875D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 4.08D+36
  real ( kind = 8 ) xden
  real ( kind = 8 ) xm1
  real ( kind = 8 ) xm2
  real ( kind = 8 ) xm4
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0D+00 .or. xbig < x ) then
    r8_gamma_log = r8_huge ( )
    return
  end if

  eps = epsilon ( eps )

  if ( x <= eps ) then

    res = - log ( x )

  else if ( x <= 1.5D+00 ) then

    if ( x < pnt68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0D+00
      xm1 = ( x - 0.5D+00 ) - 0.5D+00
    end if

    if ( x <= 0.5D+00 .or. pnt68 <= x ) then

      xden = 1.0D+00
      xnum = 0.0D+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5D+00 ) - 0.5D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0D+00 ) then

    xm2 = x - 2.0D+00
    xden = 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0D+00 ) then

    xm4 = x - 4.0D+00
    xden = - 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0D+00

    if ( x <= frtbig ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x + sqrtpi - 0.5D+00 * log ( x ) + x * ( log ( x ) - 1.0D+00 )

  end if

  r8_gamma_log = res

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a large R8.
!
!  Discussion:
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    is more suitable for this purpose.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real ( kind = 8 ) r8_huge

  r8_huge = 1.0E+30

  return
end
function r8_sech ( x )

!*****************************************************************************80
!
!! R8_SECH evaluates the hyperbolic secant, while avoiding COSH overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_SECH, the value of the function.
!
  implicit none

  real ( kind = 8 ) :: log_huge = 80.0D+00
  real ( kind = 8 ) r8_sech
  real ( kind = 8 ) x

  if ( log_huge < abs ( x ) ) then
    r8_sech = 0.0D+00
  else
    r8_sech = 1.0D+00 / cosh ( x )
  end if

  return
end
function r8_si ( x )

!*****************************************************************************80
!
!! R8_SI computes an approximation to the value of the sine integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the sine integral.
!
!    Output, real ( kind = 8 ) R8_SI, the value of the sine integral.
!
  implicit none

  integer ( kind = 4 ), parameter :: nsics = 12

  real ( kind = 8 ) absx
  real ( kind = 8 ) eta
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ), save :: nsi = 0
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_si
  real ( kind = 8 ), save, dimension ( nsics ) :: sics = (/ &
    -0.1315646598184841929D+00, -0.2776578526973601892D+00, &
     0.0354414054866659180D+00, -0.0025631631447933978D+00, &
     0.0001162365390497009D+00, -0.0000035904327241606D+00, &
     0.0000000802342123706D+00, -0.0000000013562997693D+00, &
     0.0000000000179440722D+00, -0.0000000000001908387D+00, &
     0.0000000000000016670D+00, -0.0000000000000000122D+00 /)
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: xsml = 0.0D+00

  if ( nsi == 0 ) then
    eta = 0.1D+00 * epsilon ( eta )
    call inits ( sics, nsics, eta, nsi )
    xsml = sqrt ( epsilon ( xsml ) )
  end if

  absx = abs ( x )

  if ( absx < xsml ) then
    r8_si = x
  else if ( absx <= 4.0D+00 ) then
    call csevl ( ( x**2 - 8.0D+00 ) * 0.125D+00, sics, nsi, value )
    r8_si = x * ( 0.75D+00 + value )
  else
    call r89sifg ( absx, f, g )
    r8_si = 0.5D+00 * pi - f * cos ( absx ) - g * sin ( x )
    if ( x < 0.0D+00 ) then
      r8_si = - r8_si
    end if
  end if

  return
end
subroutine r8vec_uniform ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM returns a scaled pseudorandom R8VEC
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
!    29 January 2005
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
!    Philip Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Examples:
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
