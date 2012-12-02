subroutine bromin ( n, s, tol, xr, xi, wr, wi, eps, ier )

!*****************************************************************************80
!
!! BROMIN calculates a Bromwich quadrature rule.
!
!  this subroutine calculates abscissas and weights of the
!  gaussian quadrature formula of order n for the bromwich
!  integral.  only the abscissas of the first quadrant of
!  the complex plane, the real abscissa (if n is odd) and
!  the corresponding weights are calculated.  the other
!  abscissas and weights are complex conjugates.
!
!  input parameters
!
!    n, order of the quadrature formula.
!      n must be greater than 2.
!    tol, requested relative accuracy of the abscissas.
!    s, parameters of the weight function.
!
!  output parameters
!
!    xr and xi contain the real and imaginary parts of
!      the abscissas.  if n is odd, the real abscissa
!      is xr(1).
!    wr and wi contain the real and imaginary parts of
!      the corresponding weights.
!    eps is a crude estimation of the obtained relative
!      accuracy of the abscissas.
!    ier is an error code.
!      if ier = 0, the computation is carried out to
!        the requested accuracy.
!      if ier.gt.0, the ier-th abscissa is not found.
!      if ier = -1, the computations are carried out,
!        but the requested accuracy is not
!        achieved.
!      if ier = -2, n is less than 3.
!
!  function programs required
!    function gamma(x), which evaluates the gamma
!      function for positive x.
!
  integer ( kind = 4 ) n

  real ( kind = 8 ) ak
  real ( kind = 8 ) an
  real ( kind = 8 ) arg
  real ( kind = 8 ) ci
  real ( kind = 8 ) cr
  real ( kind = 8 ) d
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dgamma
  real ( kind = 8 ) e
  real ( kind = 8 ) eps
  real ( kind = 8 ) fac
  real ( kind = 8 ) facti
  real ( kind = 8 ) factr
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ignal
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) num
  integer ( kind = 4 ) nup
  real ( kind = 8 ) pi
  real ( kind = 8 ) pr
  real ( kind = 8 ) qi
  real ( kind = 8 ) qr
  real ( kind = 8 ) ri
  real ( kind = 8 ) rr
  real ( kind = 8 ) s
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tol
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) wi(n)
  real ( kind = 8 ) wr(n)
  real ( kind = 8 ) xi(n)
  real ( kind = 8 ) xr(n)
  real ( kind = 8 ) yi
  real ( kind = 8 ) yr
  real ( kind = 8 ) z

  if ( n < 3 ) then
    ier = - 2
    return
  end if

  n1 = ( n + 1 ) / 2
  l = n - 1
  an = n
  ier = 0
  eps = tol
  arg = 0.034D+00 * ( 30.0D+00 + an + an ) / ( an - 1.0D+00 )
  factr = cos ( arg )
  facti = sin ( arg )
  fac = 1.0D+00
  ak = 0.0D+00
  do k = 1, l
    ak = ak + 1.0D+00
    fac = - fac * ak
  end do
  fac = fac * ( an + an + s - 2.0D+00 )**2 / ( an * dgamma ( an + s - 1.0D+00 ) )
!
!  Approximate the first abscissa.
!
  yr = 1.333D+00 * an + s - 1.5D+00

  if ( mod ( n, 2 ) .eq. 0 ) then
    yi = 1.6D+00 + 0.07D+00 * s
  else
    yi = 0.0D+00
  end if

  do k = 1, n1

    e = tol
    ignal = 0
    num = 0
    nup = 0
!
!  Newton-Raphson method.
!
    d = yr * yr + yi * yi
    yr = yr / d
    yi = - yi / d

    do

      qr = s * yr - 1.0D+00
      qi = s * yi
      pr = ( s + 1.0D+00 ) * ( ( s + 2.0D+00 ) * ( yr * yr - yi * yi ) &
        - 2.0D+00 * yr ) + 1.0D+00
      pi = 2.0D+00 * ( s + 1.0D+00 ) * yi * ( ( s + 2.0D+00 ) * yr - 1.0D+00 )
      z = 2.0D+00

      do j = 3, n
        rr = qr
        ri = qi
        qr = pr
        qi = pi
        z = z + 1.0D+00
        u = z + s - 2.0D+00
        v = u + z
        d = ( v * yr + ( 2.0D+00 - s ) / ( v - 2.0D+00 ) ) / u
        d1 = ( z - 1.0D+00 ) * v / ( u * ( v - 2.0D+00 ) )
        d2 = v * yi / u
        pr = ( v - 1.0D+00 ) * ( qr * d - qi * d2 ) + d1 * rr
        pi = ( v - 1.0D+00 ) * ( qi * d + qr * d2 ) + d1 * ri
      end do

      if ( ignal == 1 ) then
        exit
      end if

      d = ( yr * yr + yi * yi ) * v
      d1 = ( ( pr + qr ) * yr + ( pi + qi ) * yi ) / d + pr
      d2 = ( ( pi + qi ) * yr - ( pr + qr ) * yi ) / d + pi
      d = ( d1 * d1 + d2 * d2 ) * an
      t1 = pr * yr - pi * yi
      t2 = pi * yr + pr * yi
      cr = ( t1 * d1 + t2 * d2 ) / d
      ci = ( t2 * d1 - t1 * d2 ) / d
      yr = yr - cr
      yi = yi - ci
      num = num + 1
!
!  Test of convergence of iteration process.
!
      if ( cr * cr + ci * ci <= e * e * ( yr * yr + yi * yi ) ) then

        ignal = 1
!
!  Test of number of iteration steps.
!
      else if ( 10 < num ) then

        e = e * 10.0D+00
        ier = - 1
        nup = nup + 1

        if ( 5 < nup ) then
          ier = k
          return
        end if

      end if

    end do
!
!  Calculation of weights.
!
    if ( eps < e ) then
      eps = e
    end if

    d = qr * qr + qi * qi
    d = d * d
    d1 = yr * qr + yi * qi
    d2 = yi * qr - yr * qi
    wr(k) = fac * ( d1 * d1 - d2 * d2 ) / d
    wi(k) = 2.0D+00 * fac * d2 * d1 / d
    d = yr * yr + yi * yi
    xr(k) =   yr / d
    xi(k) = - yi / d

    if ( n1 < k + 1 ) then
      exit
    end if

    if ( n1 == k + 1 ) then
      factr = cos ( 1.5D+00 * arg )
      facti = sin ( 1.5D+00 * arg )
    end if
!
!  Approximate the (K+1)-th abscissa.
!
    yr = ( xr(k) + 0.67D+00 * an ) * factr - xi(k) * facti - 0.67D+00 * an
    yi = ( xr(k) + 0.67D+00 * an ) * facti + xi(k) * factr

  end do

  return
end
function dgamma ( x )

!*****************************************************************************80
!
!! DGAMMA calculates the GAMMA function for a real argument X.
!
!   Computation is based on an algorithm outlined in reference 1.
!   The program uses rational functions that approximate the GAMMA
!   function to at least 20 significant decimal digits.  Coefficients
!   for the approximation over the interval (1,2) are unpublished.
!   Those for the approximation for X .GE. 12 are from reference 2.
!   The accuracy achieved depends on the arithmetic system, the
!   compiler, the intrinsic functions, and proper selection of the
!   machine-dependent constants.
!
!
! Explanation of machine-dependent constants
!
! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows
! XBIG   - the largest argument for which GAMMA(X) is representable
!          in the machine, i.e., the solution to the equation
!                  GAMMA(XBIG) = beta**maxexp
! XINF   - the largest machine representable floating-point number;
!          approximately beta**maxexp
! EPS    - the smallest positive floating-point number such that
!          1.0+EPS .GT. 1.0
! XMININ - the smallest positive floating-point number such that
!          1/XMININ is machine representable
!
!     Approximate values for some important machines are:
!
!                            beta       maxexp        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! Cyber 180/855
!   under NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-Format   (D.P.)        2          127        34.844
! VAX G-Format   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! Cyber 180/855
!   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!
! Error returns
!
!  The program returns the value XINF for singularities or
!     when overflow would occur.  The computation is believed
!     to be free of underflow and overflow.
!
!
! References: "An Overview of Software Development for Special
!              Functions", W. J. Cody, Lecture Notes in Mathematics,
!              506, Numerical Analysis Dundee, 1975, G. A. Watson
!              (ed.), Springer Verlag, Berlin, 1976.
!
!              Computer Approximations, Hart, Et. Al., Wiley and
!              sons, New York, 1968.
!
!  Latest modification: October 12, 1989
!
!  Authors: W. J. Cody and L. Stoltz
!           Applied Mathematics Division
!           Argonne National Laboratory
!           Argonne, IL 60439
!
!
  INTEGER I,N
  LOGICAL PARITY
  real ( kind = 8 ) DGAMMA
  real ( kind = 8 )  C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE
  real ( kind = 8 ) TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
  DIMENSION C(7),P(8),Q(8)
!
!  Mathematical constants
!
  DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D+00,0.5D+00,12.0D+00,2.0D+00,0.0D+00/
  DATA     SQRTPI/0.9189385332046727417803297D+00/
  DATA PI/3.1415926535897932384626434D+00/
!
!  Machine dependent parameters
!
  DATA XBIG,XMININ,EPS/34.844D+00,5.88D-39,1.39D-17/
  DATA     XINF/1.70D+38/
!
!  Numerator and denominator coefficients for rational minimax
!     approximation over (1,2).
!
  DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1, &
        -3.79804256470945635097577D+2,6.29331155312818442661052D+2, &
        8.66966202790413211295064D+2,-3.14512729688483675254357D+4, &
        -3.61444134186911729807069D+4,6.64561438202405440627855D+4/

  DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2, &
       -1.01515636749021914166146D+3,-3.10777167157231109440444D+3, &
         2.25381184209801510330112D+4,4.75584627752788110767815D+3, &
       -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
!
!  Coefficients for minimax approximation over (12, INF).
!
  DATA C/-1.910444077728D-03,8.4171387781295D-04, &
      -5.952379913043012D-04,7.93650793500350248D-04, &
      -2.777777777777681622553D-03,8.333333333333333331554247D-02, &
       5.7083835261D-03/
!
!  Statement functions for conversion between integer and float
!
  CONV(I) = DBLE(I)
  PARITY = .FALSE.
  FACT = ONE
  N = 0
  Y = X
  IF (Y .LE. ZERO) THEN
!
!  Argument is negative
!
        Y = -X
        Y1 = AINT(Y)
        RES = Y - Y1
        IF (RES .NE. ZERO) THEN
              IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
              FACT = -PI / SIN(PI*RES)
              Y = Y + ONE
           ELSE
              RES = XINF
              DGAMMA = RES
              RETURN
        END IF
  END IF
!
!  Argument is positive
!
  IF (Y .LT. EPS) THEN
!
!  Argument .LT. EPS
!
        IF (Y .GE. XMININ) THEN
              RES = ONE / Y
           ELSE
              RES = XINF
              DGAMMA = RES
              RETURN
        END IF
     ELSE IF (Y .LT. TWELVE) THEN
        Y1 = Y
        IF (Y .LT. ONE) THEN
!
!  0.0 .LT. argument .LT. 1.0
!
              Z = Y
              Y = Y + ONE
           ELSE
!
!  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
!
              N = INT(Y) - 1
              Y = Y - CONV(N)
              Z = Y - ONE
        END IF
!
!  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
!
        XNUM = ZERO
        XDEN = ONE
        DO I = 1, 8
           XNUM = (XNUM + P(I)) * Z
           XDEN = XDEN * Z + Q(I)
        end do
        RES = XNUM / XDEN + ONE
        IF (Y1 .LT. Y) THEN
!
!  Adjust result for case  0.0 .LT. argument .LT. 1.0
!
              RES = RES / Y1
           ELSE IF (Y1 .GT. Y) THEN
!
!  Adjust result for case  2.0 .LT. argument .LT. 12.0
!
              DO I = 1, N
                 RES = RES * Y
                 Y = Y + ONE
              end do
        END IF
     ELSE
!
!  Evaluate for argument .GE. 12.0,
!
        IF (Y .LE. XBIG) THEN
              YSQ = Y * Y
              SUM = C(7)
              DO I = 1, 6
                 SUM = SUM / YSQ + C(I)
              end do
              SUM = SUM/Y - Y + SQRTPI
              SUM = SUM + (Y-HALF)*LOG(Y)
              RES = EXP(SUM)
           ELSE
              RES = XINF
              DGAMMA = RES
              RETURN
        END IF
  END IF
!
!  Final adjustments and return
!
  IF (PARITY) RES = -RES
  IF (FACT .NE. ONE) RES = FACT / RES
  DGAMMA = RES

  RETURN
END
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
