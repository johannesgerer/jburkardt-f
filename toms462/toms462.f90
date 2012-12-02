function bivnor ( ah, ak, r )

!*****************************************************************************80
!
!! BIVNOR computes the bivariate normal CDF.
!
!  Discussion:
!
!    BIVNOR computes the probability for two normal variates X and Y
!    whose correlation is R, that AH <= X and AK <= Y.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2012
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Donnelly.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Donnelly,
!    Algorithm 462: Bivariate Normal Distribution,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 10, page 638.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AH, AK, the lower limits of integration.
!
!    Input, real ( kind = 8 ) R, the correlation between X and Y.
!
!    Output, real ( kind = 8 ) BIVNOR, the bivariate normal CDF.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) IDIG, the number of significant digits
!    to the right of the decimal point desired in the answer.
!
  implicit none

  real ( kind = 8 ) a2
  real ( kind = 8 ) ah
  real ( kind = 8 ) ak
  real ( kind = 8 ) ap
  real ( kind = 8 ) b
  real ( kind = 8 ) bivnor
  real ( kind = 8 ) cn
  real ( kind = 8 ) con
  real ( kind = 8 ) conex
  real ( kind = 8 ) ex
  real ( kind = 8 ) g2
  real ( kind = 8 ) gauss
  real ( kind = 8 ) gh
  real ( kind = 8 ) gk
  real ( kind = 8 ) gw
  real ( kind = 8 ) h2
  real ( kind = 8 ) h4
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: idig = 15
  integer ( kind = 4 ) is
  real ( kind = 8 ) r
  real ( kind = 8 ) rr
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) sgn
  real ( kind = 8 ) sn
  real ( kind = 8 ) sp
  real ( kind = 8 ) sqr
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: twopi = 6.283185307179587D+00
  real ( kind = 8 ) w2
  real ( kind = 8 ) wh
  real ( kind = 8 ) wk

  gauss ( t ) = ( 1.0D+00 + erf ( t / sqrt ( 2.0D+00 ) ) ) / 2.0D+00

!  GAUSS is a univariate lower normal tail area calculated here from the
!  central error function ERF.
!
  b = 0.0D+00

  gh = gauss ( - ah ) / 2.0D+00
  gk = gauss ( - ak ) / 2.0D+00

  if ( r == 0.0D+00 ) then
    b = 4.0D+000 * gh * gk
    b = max ( b, 0.0D+00 )
    b = min ( b, 1.0D+00 )
    bivnor = b
    return
  end if

  rr = ( 1.0D+00 + r ) * ( 1.0D+00 - r )

  if ( rr < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BIVNOR - Fatal error!'
    write ( *, '(a)' ) '  1 < |R|.'
    stop
  end if

  if ( rr == 0.0D+00 ) then

    if ( r < 0.0D+00 ) then

      if ( ah + ak < 0.0D+00 ) then
        b = 2.0D+00 * ( gh + gk ) - 1.0D+00
      end if

    else

      if ( ah - ak < 0.0D+00 ) then
        b = 2.0D+00 * gk
      else
        b = 2.0D+00 * gh
      end if

    end if

    b = max ( b, 0.0D+00 )
    b = min ( b, 1.0D+00 )
    bivnor = b
    return

  end if

  sqr = sqrt ( rr )

  if ( idig == 15 ) then
    con = twopi * 1.0D-15 / 2.0D+00
  else
    con = twopi / 2.0D+00
    do i = 1, idig
      con = con / 10.0D+00
    end do
  end if
!
!  (0,0)
!
  if ( ah == 0.0D+00 .and. ak == 0.0D+00 ) then
    b = 0.25D+00 + asin ( r ) / twopi
    b = max ( b, 0.0D+00 )
    b = min ( b, 1.0D+00 )
    bivnor = b
    return
  end if
!
!  (0,nonzero)
!
  if ( ah == 0.0D+00 .and. ak /= 0.0D+00 ) then

    b = gk
    wh = -ak
    wk = ( ah / ak - r ) / sqr
    gw = 2.0D+00 * gk
    is = 1
!
!  (nonzero,0)
!
  else if ( ah /= 0.0D+00 .and. ak == 0.0D+00 ) then

    b = gh
    wh = -ah
    wk = ( ak / ah - r ) / sqr
    gw = 2.0D+00 * gh
    is = -1
!
!  (nonzero,nonzero)
!
  else if ( ah /= 0.0 .and. ak /= 0.0 ) then

    b = gh + gk
    if ( ah * ak < 0.0D+00 ) then
      b = b - 0.5D+00
    end if
    wh = - ah
    wk = ( ak / ah - r ) / sqr
    gw = 2.0D+00 * gh
    is = -1

  end if

  do

    sgn = -1.0D+00
    t = 0.0D+00

    if ( wk /= 0.0D+00 ) then

      if ( abs ( wk ) == 1.0D+00 ) then

        t = wk * gw * ( 1.0D+00 - gw ) / 2.0D+00
        b = b + sgn * t

      else

        if ( 1.0D+00 < abs ( wk ) ) then

          sgn = -sgn
          wh = wh * wk
          g2 = gauss ( wh )
          wk = 1.0D+00 / wk

          if ( wk < 0.0D+00 ) then
            b = b + 0.5D+00
          end if

          b = b - ( gw + g2 ) / 2.0D+00 + gw * g2

        end if

        h2 = wh * wh
        a2 = wk * wk
        h4 = h2 / 2.0D+00
        ex = exp ( - h4 )
        w2 = h4 * ex
        ap = 1.0D+00
        s2 = ap - ex
        sp = ap
        s1 = 0.0D+00
        sn = s1
        conex = abs ( con / wk )

        do

          cn = ap * s2 / ( sn + sp )
          s1 = s1 + cn

          if ( abs ( cn ) <= conex ) then
            exit
          end if

          sn = sp
          sp = sp + 1.0D+00
          s2 = s2 - w2
          w2 = w2 * h4 / sp
          ap = - ap * a2

        end do

        t = ( atan ( wk ) - wk * s1 ) / twopi
        b = b + sgn * t

      end if

    end if

    if ( 0 <= is ) then
      exit
    end if

    if ( ak == 0.0D+00 ) then
      exit
    end if

    wh = -ak
    wk = ( ah / ak - r ) / sqr
    gw = 2.0D+00 * gk
    is = 1

  end do

  b = max ( b, 0.0D+00 )
  b = min ( b, 1.0D+00 )
  bivnor = b

  return
end
subroutine bivariate_normal_cdf_values ( n_data, x, y, r, fxy )

!*****************************************************************************80
!
!! BIVARIATE_NORMAL_CDF_VALUES returns some values of the bivariate normal CDF.
!
!  Discussion:
!
!    FXY is the probability that two variables A and B, which are
!    related by a bivariate normal distribution with correlation R,
!    respectively satisfy A <= X and B <= Y.
!
!    Mathematica can evaluate the bivariate normal CDF via the commands:
!
!      <<MultivariateStatistics`
!      cdf = CDF[MultinormalDistribution[{0,0}{{1,r},{r,1}}],{x,y}]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    National Bureau of Standards,
!    Tables of the Bivariate Normal Distribution and Related Functions,
!    Applied Mathematics Series, Number 50, 1959.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, Y, the parameters of the function.
!
!    Output, real ( kind = 8 ) R, the correlation value.
!
!    Output, real ( kind = 8 ) FXY, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 41

  real ( kind = 8 ) fxy
  real ( kind = 8 ), save, dimension ( n_max ) :: fxy_vec = (/ &
  0.02260327218569867D+00, &
  0.1548729518584100D+00, &
  0.4687428083352184D+00, &
  0.7452035868929476D+00, &
  0.8318608306874188D+00, &
  0.8410314261134202D+00, &
  0.1377019384919464D+00, &
  0.1621749501739030D+00, &
  0.1827411243233119D+00, &
  0.2010067421506235D+00, &
  0.2177751155265290D+00, &
  0.2335088436446962D+00, &
  0.2485057781834286D+00, &
  0.2629747825154868D+00, &
  0.2770729823404738D+00, &
  0.2909261168683812D+00, &
  0.3046406378726738D+00, &
  0.3183113449213638D+00, &
  0.3320262544108028D+00, &
  0.3458686754647614D+00, &
  0.3599150462310668D+00, &
  0.3742210899871168D+00, &
  0.3887706405282320D+00, &
  0.4032765198361344D+00, &
  0.4162100291953678D+00, &
  0.6508271498838664D+00, &
  0.8318608306874188D+00, &
  0.0000000000000000D+00, &
  0.1666666666539970D+00, &
  0.2500000000000000D+00, &
  0.3333333333328906D+00, &
  0.5000000000000000D+00, &
  0.7452035868929476D+00, &
  0.1548729518584100D+00, &
  0.1548729518584100D+00, &
  0.06251409470431653D+00, &
  0.7452035868929476D+00, &
  0.1548729518584100D+00, &
  0.1548729518584100D+00, &
  0.06251409470431653D+00, &
  0.6337020457912916D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r
  real ( kind = 8 ), save, dimension ( n_max ) :: r_vec = (/ &
     0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00, &
     0.500D+00, -0.900D+00, -0.800D+00, -0.700D+00, -0.600D+00, &
    -0.500D+00, -0.400D+00, -0.300D+00, -0.200D+00, -0.100D+00, &
     0.000D+00,  0.100D+00,  0.200D+00,  0.300D+00,  0.400D+00, &
     0.500D+00,  0.600D+00,  0.700D+00,  0.800D+00,  0.900D+00, &
     0.673D+00,  0.500D+00, -1.000D+00, -0.500D+00,  0.000D+00, &
     0.500D+00,  1.000D+00,  0.500D+00,  0.500D+00,  0.500D+00, &
     0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00, &
     0.500D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -2.0D+00, -1.0D+00,  0.0D+00,  1.0D+00,  2.0D+00, &
     3.0D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
    -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
    -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
    -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
     1.0D+00,  2.0D+00,  0.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00,  1.0D+00, -1.0D+00, &
    -1.0D+00,  1.0D+00,  1.0D+00, -1.0D+00, -1.0D+00, &
     0.7071067811865475D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( n_max ) :: y_vec = (/ &
     1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, &
     1.0D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  1.0D+00,  0.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00, -1.0D+00,  1.0D+00, &
    -1.0D+00,  1.0D+00, -1.0D+00,  1.0D+00, -1.0D+00, &
     0.7071067811865475D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    r = 0.0D+00
    x = 0.0D+00
    y = 0.0D+00
    fxy = 0.0D+00
  else
    r = r_vec(n_data)
    x = x_vec(n_data)
    y = y_vec(n_data)
    fxy = fxy_vec(n_data)
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
