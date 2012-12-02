function arc_sine ( s )

!*****************************************************************************80
!
!! ARC_SINE computes the arc sine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ASIN routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant 
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) S, the argument.
!
!    Output, real ( kind = 8 ) ARC_SINE, an angle whose sine is S.
!
  implicit none

  real ( kind = 8 ) arc_sine
  real ( kind = 8 ) s
  real ( kind = 8 ) s2

  s2 = s
  s2 = max ( s2, -1.0D+00 )
  s2 = min ( s2, +1.0D+00 )

  arc_sine = asin ( s2 )

  return
end
subroutine ball_f1_nd ( func, n, center, r, result )

!*****************************************************************************80
!
!! BALL_F1_ND approximates an integral inside a ball in ND.
!
!  Integration region:
!
!    sum ( X(1:N) - CENTER(1:N) )^2 <= R * R.
!
!  Discussion:
!
!    An (N+1)*2^N point 5-th degree formula is used, Stroud number SN:5-6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F at the N-vector X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) CENTER(N), the center of the ball.
!
!    Input, real ( kind = 8 ) R, the radius of the ball.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ball_volume_nd
  real ( kind = 8 ) center(n)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) ktemp
  real ( kind = 8 ) quad
  real ( kind = 8 ) r
  real ( kind = 8 ) result
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) u
  real ( kind = 8 ) u2
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y

  if ( r == 0.0D+00 ) then
    result = 0.0D+00
    return
  end if

  u2 = ( 1.0D+00 - 2.0D+00 * sqrt ( 1.0D+00 / real ( n + 4, kind = 8 ) ) ) &
    / real ( n + 2, kind = 8 )
  u = sqrt ( u2 )
  x(1:n) = center(1:n) - r * u

  w = 1.0D+00 / real ( ( n + 1 ) * 2**n, kind = 8 )

  quad = 0.0D+00
  ihi = 2**n

  do i = 1, ihi

    itemp = i - 1

    do j = 1, n

      u = ( center(j) - x(j) ) / r

      if ( mod ( itemp, 2 ) == 1 ) then
        x(j) = center(j) - abs ( x(j) - center(j) )
      else
        x(j) = center(j) + abs ( x(j) - center(j) )
      end if

      itemp = itemp / 2

    end do

    quad = quad + w * func ( n, x )

  end do

  temp = sqrt ( real ( n + 4, kind = 8 ) )

  t = sqrt ( 2.0D+00 * real ( n + 1, kind = 8 ) / real ( n + 2, kind = 8 ) ) &
    / ( real ( n, kind = 8 ) * temp )

  y = ( 1.0D+00 + 2.0D+00 / ( real ( n, kind = 8 ) * temp ) ) &
    / real ( n + 2, kind = 8 )
  v = sqrt ( y - t )
  u = sqrt ( y + real ( n - 1, kind = 8 ) * t )

  khi = 2**n

  do i = 1, n

    x(1:n) = center(1:n) - r * v

    x(i) = center(i) - r * u

    do k = 1, khi

      ktemp = k - 1

      do j = 1, n

        if ( mod ( ktemp, 2 ) == 1 ) then
          x(j) = center(j) - abs ( x(j) - center(j) )
        else
          x(j) = center(j) + abs ( x(j) - center(j) )
        end if

        ktemp = ktemp / 2

      end do

      quad = quad + w * func ( n, x )

    end do

    x(i) = center(i) - r * v

  end do

  volume = ball_volume_nd ( n, r )
  result = quad * volume

  return
end
subroutine ball_f3_nd ( func, n, center, r, result )

!*****************************************************************************80
!
!! BALL_F3_ND approximates an integral inside a ball in ND.
!
!  Integration region:
!
!    sum ( X(1:N) - CENTER(1:N) )^2 <= R * R.
!
!  Discussion:
!
!    A 2^(N+1)-1 point 5-th degree formula is used, Stroud number SN:5-4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F at the N-vector X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) CENTER(N), the center of the ball.
!
!    Input, real ( kind = 8 ) R, the radius of the ball.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ball_volume_nd
  real ( kind = 8 ) center(n)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtemp
  integer ( kind = 4 ) k
  real ( kind = 8 ) quad
  real ( kind = 8 ) r
  real ( kind = 8 ) result
  real ( kind = 8 ) ri
  real ( kind = 8 ) s
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight
  real ( kind = 8 ) x(n)

  if ( r == 0.0D+00 ) then
    result = 0.0D+00
    return
  end if

  quad = 0.0D+00
!
!  The first point is the center of the ball.
!
  x(1:n) = center(1:n)

  weight = 4.0D+00 / real ( ( n + 2 ) * ( n + 2 ), kind = 8 )
  quad = quad + weight * func ( n, x )

  s = 1.0D+00 / sqrt ( real ( n + 4, kind = 8 ) )

  do i = 1, n

    ri = sqrt ( real ( i + 2, kind = 8 ) / real ( n + 4, kind = 8 ) )
!
!  Set up the first point, with (I-1) zeroes, RI, and then N-I S's.
!
    do j = 1, n

      if ( j < i ) then
        x(j) = center(j)
      else if ( j == i ) then
        x(j) = center(j) + r * ri
      else
        x(j) = center(j) + r * s
      end if

    end do

    weight = 2.0D+00**( i - n ) * real ( n + 4, kind = 8 ) &
      / real ( ( i + 1 ) * ( i + 2 ) * ( n + 2 ), kind = 8 )
!
!  Now go through all sign permutations of the basic point.
!
    do j = 1, 2**(n+1-i)

      jtemp = j - 1

      do k = i, n

        if ( mod ( jtemp, 2 ) == 1 ) then
          x(k) = center(k) - abs ( x(k) - center(k) )
        else
          x(k) = center(k) + abs ( x(k) - center(k) )
        end if

        jtemp = jtemp / 2

      end do

      quad = quad + weight * func ( n, x )

    end do

  end do

  volume = ball_volume_nd ( n, r )
  result = quad * volume

  return
end
function ball_monomial_nd ( n, p, r )

!*****************************************************************************80
!
!! BALL_MONOMIAL_ND integrates a monomial on a ball in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) <= R * R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gerald Folland,
!    How to Integrate a Polynomial Over a Sphere,
!    American Mathematical Monthly,
!    Volume 108, Number 5, May 2001, pages 446-448.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, integer ( kind = 4 ) P(N), the exponents of X(1) through X(N) 
!    in the monomial.  The exponents P(N) must be nonnegative.
!
!    Input, real ( kind = 8 ) R, the radius of the ball.
!
!    Output, real ( kind = 8 ) BALL_MONOMIAL_ND, the integral of
!    X1^P(1) * X2^P(2) * ... * XN^P(N) over the ball.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ball_monomial_nd
  integer ( kind = 4 ) p(n)
  real ( kind = 8 ) power
  real ( kind = 8 ) r
  real ( kind = 8 ) sphere_unit_monomial_nd

  power = real ( sum ( p ) + n, kind = 8 )

  ball_monomial_nd = sphere_unit_monomial_nd ( n, p ) * r**power / power

  return
end
subroutine ball_unit_07_3d ( func, result )

!*****************************************************************************80
!
!! BALL_UNIT_07_3D approximates an integral inside the unit ball in 3D.
!
!  Integration region:
!
!    X*X + Y*Y + Z*Z <= 1.
!
!  Discussion:
!
!    A 64 point 7-th degree formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 4

  real ( kind = 8 ) angle
  real ( kind = 8 ) ball_unit_volume_3d
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) weight1(order)
  real ( kind = 8 ) weight2(order)
  real ( kind = 8 ) weight3(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab1(order)
  real ( kind = 8 ) xtab2(order)
  real ( kind = 8 ) xtab3(order)
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
!  This is the 5 point Gauss-Legendre rule,
!  but with the midpoint deleted, and with different weights.
!
  xtab1(1:4) = (/ &
    -0.906179845938663992797626878299D+00, & 
    -0.538469310105683091036314420700D+00, &
     0.538469310105683091036314420700D+00, &
     0.906179845938663992797626878299D+00 /)

  weight1(1:4) = (/ &
    0.19455533421780251826D+00, &
    0.13877799911553081506D+00, &
    0.13877799911553081506D+00, &
    0.19455533421780251826D+00 /)
!
!  Set XTAB2 and WEIGHT2.
!
  do j = 1, order
    angle = pi * real ( 2 * j - 1, kind = 8 ) &
      / real ( 2 * order, kind = 8 )
    xtab2(j) = cos ( angle )
  end do

  weight2(1:order) = 1.0D+00
!
!  Set XTAB3 and WEIGHT3 for the interval [-1,1].
!
  call legendre_set ( order, xtab3, weight3 )

  w = 3.0D+00 / 16.0D+00

  quad = 0.0D+00

  do i = 1, order
    do j = 1, order
      do k = 1, order

        x = xtab1(i) * sqrt ( 1.0D+00 - xtab2(j) * xtab2(j) ) &
                     * sqrt ( 1.0D+00 - xtab3(k) * xtab3(k) )
        y = xtab1(i) * xtab2(j) * sqrt ( 1.0D+00 - xtab3(k) * xtab3(k) )
        z = xtab1(i) * xtab3(k)

        quad = quad + w * weight1(i) * weight2(j) * weight3(k) &
          * func ( x, y, z )

      end do
    end do
  end do

  volume = ball_unit_volume_3d ( )
  result = quad * volume

  return
end
subroutine ball_unit_14_3d ( func, result )

!*****************************************************************************80
!
!! BALL_UNIT_14_3D approximates an integral inside the unit ball in 3D.
!
!  Integration region:
!
!    X*X + Y*Y + Z*Z <= 1.
!
!  Discussion:
!
!    A 288 point 14-th degree formula is used, Stroud number S3:14-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  real ( kind = 8 ) ball_unit_volume_3d
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) quad
  real ( kind = 8 ), save, dimension ( 4 ) :: r = (/ &
    0.968160240D+00, 0.836031107D+00, 0.613371433D+00, 0.324253423D+00 /)
  real ( kind = 8 ) result
  real ( kind = 8 ) temp
  real ( kind = 8 ) volume
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ), save, dimension ( 4 ) :: weight = (/ &
    0.076181268D+00, 0.126263673D+00, 0.098048133D+00, 0.032840260D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( 5 ) :: xtab = (/ &
    -0.151108275D+00, 0.315838353D+00, 0.346307112D+00, &
    -0.101808787D+00, -0.409228403D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( 5 ) :: ytab = (/ &
    0.155240600D+00, 0.257049387D+00, 0.666277790D+00, &
    0.817386065D+00, 0.501547712D+00 /)
  real ( kind = 8 ) z
  real ( kind = 8 ), save, dimension ( 5 ) :: ztab = (/ &
    0.976251323D+00, 0.913330032D+00, 0.660412970D+00, &
    0.567022920D+00, 0.762221757D+00 /)

  quad = 0.0D+00

  do m = 1, 4

    w1 = 125.0D+00 * weight(m) / 3360.0D+00
    x = 0.525731112D+00 * r(m)
    y = 0.850650808D+00 * r(m)
    z = 0.0D+00

    do j = 1, 2
      x = -x
      do k = 1, 2
        y = -y
        do l = 1, 3
          call r8_swap3 ( x, y, z )
          quad = quad + w1 * func ( x, y, z )
        end do
      end do
    end do

    w2 = 143.0D+00 * weight(m) / 3360.0D+00

    do n = 1, 5

      x = xtab(n) * r(m)
      y = ytab(n) * r(m)
      z = ztab(n) * r(m)

      do i = 1, 3

        temp = x
        x = z
        z = -y
        y = -temp

        do j = 1, 3

          call r8_swap3 ( x, y, z )

          quad = quad + w2 * func ( x, y, z )

        end do

        y = -y
        z = -z
        quad = quad + w2 * func ( x, y, z )

      end do

    end do

  end do

  volume = ball_unit_volume_3d ( )
  result = quad * volume

  return
end
subroutine ball_unit_15_3d ( func, result )

!*****************************************************************************80
!
!! BALL_UNIT_15_3D approximates an integral inside the unit ball in 3D.
!
!  Integration region:
!
!    X*X + Y*Y + Z*Z <= 1.
!
!  Discussion:
!
!    A 512 point 15-th degree formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order1 = 4
  integer ( kind = 4 ), parameter :: order2 = 8

  real ( kind = 8 ) ball_unit_volume_3d
  real ( kind = 8 ) cj
  real ( kind = 8 ) ck
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) sj
  real ( kind = 8 ) sk
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ), save, dimension ( order1 ) :: weight1 = (/ &
    0.0328402599D+00, 0.0980481327D+00, 0.1262636728D+00, 0.0761812678D+00 /)
  real ( kind = 8 ) weight2(order2)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( order1 ) :: xtab1 = (/ &
    0.3242534234D+00, 0.6133714327D+00, 0.8360311073D+00, 0.9681602395D+00 /)
  real ( kind = 8 ) xtab2(order2)
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  call legendre_set ( order2, xtab2, weight2 )

  w = 3.0D+00 / 32.0D+00

  quad = 0.0D+00

  do i = 1, order1

    do j = 1, order2

      sj = xtab2(j)
      cj = sqrt ( 1.0D+00 - sj * sj )

      do k = 1, 16
        sk = sin ( real ( k, kind = 8 ) * pi / 8.0D+00 )
        ck = cos ( real ( k, kind = 8 ) * pi / 8.0D+00 )
        x = xtab1(i) * cj * ck
        y = xtab1(i) * cj * sk
        z = xtab1(i) * sj
        quad = quad + w * weight1(i) * weight2(j) * func ( x, y, z )
      end do

    end do

  end do

  volume = ball_unit_volume_3d ( )
  result = quad * volume

  return
end
subroutine ball_unit_f1_nd ( func, n, result )

!*****************************************************************************80
!
!! BALL_UNIT_F1_ND approximates an integral inside the unit ball in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) <= 1.
!
!  Discussion:
!
!    An (N+1)*2^N point 5-th degree formula is used, Stroud number SN:5-6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F at the N-vector X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ball_unit_volume_nd
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) ktemp
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) u
  real ( kind = 8 ) u2
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y

  u2 = ( 1.0D+00 - 2.0D+00 * sqrt ( 1.0D+00 / real ( n + 4, kind = 8 ) ) ) &
    / real ( n + 2, kind = 8 )
  u = sqrt ( u2 )
  x(1:n) = -u

  w = 1.0D+00 / real ( ( n + 1 ) * 2**n, kind = 8 )

  quad = 0.0D+00
  ihi = 2**n

  do i = 1, ihi

    itemp = i - 1

    do j = 1, n

      if ( mod ( itemp, 2 ) == 1 ) then
        x(j) = -abs ( x(j) )
      else
        x(j) = abs ( x(j) )
      end if

      itemp = itemp / 2

    end do

    quad = quad + w * func ( n, x )

  end do

  temp = sqrt ( real ( n + 4, kind = 8 ) )

  t = sqrt ( 2.0D+00 * real ( n + 1, kind = 8 ) / real ( n + 2, kind = 8 ) ) &
    / ( real ( n, kind = 8 ) * temp )

  y = ( 1.0D+00 + 2.0D+00 / ( real ( n, kind = 8 ) * temp ) ) &
    / real ( n + 2, kind = 8 )
  v = sqrt ( y - t )
  u = sqrt ( y + real ( n - 1, kind = 8 ) * t )

  khi = 2**n

  do i = 1, n

    x(1:n) = -v

    x(i) = -u

    do k = 1, khi

      ktemp = k - 1

      do j = 1, n

        if ( mod ( ktemp, 2 ) == 1 ) then
          x(j) = -abs ( x(j) )
        else
          x(j) = abs ( x(j) )
        end if

        ktemp = ktemp / 2

      end do

      quad = quad + w * func ( n, x )

    end do

    x(i) = -v

  end do

  volume = ball_unit_volume_nd ( n )
  result = quad * volume

  return
end
subroutine ball_unit_f3_nd ( func, n, result )

!*****************************************************************************80
!
!! BALL_UNIT_F3_ND approximates an integral inside the unit ball in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) <= 1.
!
!  Discussion:
!
!    A 2^(N+1)-1 point 5-th degree formula is used, Stroud number SN:5-4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F at the N-vector X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ball_unit_volume_nd
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtemp
  integer ( kind = 4 ) k
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) ri
  real ( kind = 8 ) s
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight
  real ( kind = 8 ) x(n)

  quad = 0.0D+00
!
!  The first point is the center of the ball.
!
  x(1:n) = 0.0D+00

  weight = 4.0D+00 / real ( ( n + 2 ) * ( n + 2 ), kind = 8 )
  quad = quad + weight * func ( n, x )

  s = 1.0D+00 / sqrt ( real ( n + 4, kind = 8 ) )

  do i = 1, n

    ri = sqrt ( real ( i + 2, kind = 8 ) / real ( n + 4, kind = 8 ) )
!
!  Set up the first point, with (I-1) zeroes, RI, and then N-I S's.
!
    do j = 1, n

      if ( j < i ) then
        x(j) = 0.0D+00
      else if ( j == i ) then
        x(j) = ri
      else
        x(j) = s
      end if

    end do

    weight = 2.0D+00**( i - n ) * real ( n + 4, kind = 8 ) &
      / real ( ( i + 1 ) * ( i + 2 ) * ( n + 2 ), kind = 8 )
!
!  Now go through all sign permutations of the basic point.
!
    do j = 1, 2**(n+1-i)

      jtemp = j - 1

      do k = i, n

        if ( mod ( jtemp, 2 ) == 1 ) then
          x(k) = -abs ( x(k) )
        else
          x(k) = abs ( x(k) )
        end if

        jtemp = jtemp / 2

      end do

      quad = quad + weight * func ( n, x )

    end do

  end do

  volume = ball_unit_volume_nd ( n )
  result = quad * volume

  return
end
function ball_unit_volume_3d ( )

!*****************************************************************************80
!
!! BALL_UNIT_VOLUME_3D computes the volume of the unit ball in 3D.
!
!  Integration region:
!
!    X*X + Y*Y + Z*Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) BALL_UNIT_VOLUME_3D, the volume of the ball.
!
  implicit none

  real ( kind = 8 ) ball_unit_volume_3d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  ball_unit_volume_3d = ( 4.0D+00 / 3.0D+00 ) * pi

  return
end
function ball_unit_volume_nd ( n )

!*****************************************************************************80
!
!! BALL_UNIT_VOLUME_ND computes the volume of the unit ball in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) <= 1.
!
!  Discussion:
!
!    N  Volume
!
!    2             PI
!    3  (4/3)    * PI
!    4  (1/2)    * PI^2
!    5  (8/15)   * PI^2
!    6  (1/6)    * PI^3
!    7  (16/105) * PI^3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) BALL_UNIT_VOLUME_ND, the volume of the ball.
!
  implicit none

  real ( kind = 8 ) ball_unit_volume_nd
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) volume

  if ( mod ( n, 2 ) == 0 ) then
    m = n / 2
    volume = ( pi )**m
    do i = 1, m
      volume = volume / real ( i, kind = 8 )
    end do
  else
    m = ( n - 1 ) / 2
    volume = ( pi )**m * 2.0D+00**n
    do i = m+1, 2*m+1
      volume = volume / real ( i, kind = 8 )
    end do
  end if

  ball_unit_volume_nd = volume

  return
end
function ball_volume_3d ( r )

!*****************************************************************************80
!
!! BALL_VOLUME_3D computes the volume of a ball in 3D.
!
!  Integration region:
!
!    X*X + Y*Y + Z*Z <= R * R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the ball.
!
!    Output, real ( kind = 8 ) BALL_VOLUME_3D, the volume of the ball.
!
  implicit none

  real ( kind = 8 ) ball_volume_3d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r

  ball_volume_3d = ( 4.0D+00 / 3.0D+00 ) * pi * r**3

  return
end
function ball_volume_nd ( n, r )

!*****************************************************************************80
!
!! BALL_VOLUME_ND computes the volume of a ball in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) <= R * R
!
!  Discussion:
!
!    N  Volume
!
!    2             PI   * R^2
!    3  (4/3)    * PI   * R^3
!    4  (1/2)    * PI^2 * R^4
!    5  (8/15)   * PI^2 * R^5
!    6  (1/6)    * PI^3 * R^6
!    7  (16/105) * PI^3 * R^7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) R, the radius of the ball.
!
!    Output, real ( kind = 8 ) BALL_VOLUME_ND, the volume of the ball.
!
  implicit none

  real ( kind = 8 ) ball_unit_volume_nd
  real ( kind = 8 ) ball_volume_nd
  integer ( kind = 4 ) n
  real ( kind = 8 ) r

  ball_volume_nd = ball_unit_volume_nd ( n ) * r**n

  return
end
subroutine c1_geg_monomial_integral ( alpha, expon, value )

!*****************************************************************************80
!
!! C1_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on C1.
!
!  Discussion:
!
!    C1_GEG is the interval [-1,+1] with the Gegenbauer weight function
!
!      w(alpha;x) = (1-x^2)^alpha
!
!    with -1.0 < alpha.
!
!    value = integral ( -1 <= x <= +1 ) x^expon (1-x^2)^alpha dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X^2).
!    - 1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) c
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) value
  real ( kind = 8 ) value1

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C1_GEG_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( mod ( expon, 2 ) == 1 ) then
    value = 0.0D+00
    return
  end if

  c = real ( expon, kind = 8 )

  arg1 = - alpha
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + alpha + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )

  value = 2.0D+00 * r8_gamma ( 1.0D+00 + c ) * r8_gamma ( 1.0D+00 + alpha ) &
    * value1 / r8_gamma ( 2.0D+00 + alpha  + c )

  return
end
subroutine c1_jac_monomial_integral ( alpha, beta, expon, value )

!*****************************************************************************80
!
!! C1_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over C1.
!
!  Discussion:
!
!    value = integral ( -1 <= x <= +1 ) x^expon (1-x)^alpha (1+x)^beta dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X) in the weight factor.
!
!    Input, real ( kind = 8 ) BETA, the exponent of (1+X) in the weight factor.
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) beta
  real ( kind = 8 ) c
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) s
  real ( kind = 8 ) value
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2

  c = real ( expon, kind = 8 )

  if ( mod ( expon, 2 ) == 0 ) then
    s = +1.0D+00
  else
    s = -1.0D+00
  end if

  arg1 = - alpha
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + beta + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )

  arg1 = - beta
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + alpha + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value2 )

  value = r8_gamma ( 1.0D+00 + c ) * ( &
      s * r8_gamma ( 1.0D+00 + beta  ) * value1 &
    / r8_gamma ( 2.0D+00 + beta  + c ) &
    +     r8_gamma ( 1.0D+00 + alpha ) * value2 &
    / r8_gamma ( 2.0D+00 + alpha + c ) )

  return
end
subroutine c1_leg_monomial_integral ( expon, value )

!*****************************************************************************80
!
!! C1_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on C1.
!
!  Discussion:
!
!    C1_LEG is the interval [-1,+1] with the Legendre weight function
!
!      w(x) = 1.
!
!    value = integral ( -1 <= x <= +1 ) x^expon dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) expon
  real ( kind = 8 ) value

  if ( mod ( expon, 2 ) == 1 ) then
    value = 0.0D+00
    return
  end if

  if ( expon < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C1_LEG_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  EXPON < 0.'
    stop
  end if

  value = 2.0D+00 / real ( expon + 1, kind = 8 )

  return
end
subroutine circle_annulus ( func, center, radius1, radius2, nr, result )

!*****************************************************************************80
!
!! CIRCLE_ANNULUS approximates an integral in an annulus.
!
!  Integration region:
!
!    RADIUS1^2 <= ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= RADIUS2^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Peirce,
!    Numerical Integration Over the Planar Annulus,
!    Journal of the Society for Industrial and Applied Mathematics,
!    Volume 5, Number 2, June 1957, pages 66-73.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function of two
!    variables which is to be integrated, of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, real ( kind = 8 ) CENTER(2), the center of the circle.
!
!    Input, real ( kind = 8 ) RADIUS1, RADIUS2, the radii of the circles.
!
!    Input, integer ( kind = 4 ) NR, the order of the rule.  This quantity 
!    specifies the number of distinct radii to use.  The number of angles used 
!    will be 4*NR, for a total of 4*NR*NR points.
!
!    Output, real ( kind = 8 ) RESULT, the approximation to the integral.
!
  implicit none

  integer ( kind = 4 ) nr
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) a
  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) circle_annulus_area_2d
  real ( kind = 8 ) d
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nt
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) ra(nr)
  real ( kind = 8 ) radius1
  real ( kind = 8 ) radius2
  real ( kind = 8 ) result
  real ( kind = 8 ) rw(nr)
  real ( kind = 8 ) t
  real ( kind = 8 ) tw
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Choose radial abscissas and weights.
!
  call legendre_set ( nr, ra, rw )
  a = -1.0D+00
  b = +1.0D+00
  c = radius1 * radius1
  d = radius2 * radius2
  call rule_adjust ( a, b, c, d, nr, ra, rw )
  ra(1:nr) = sqrt ( ra(1:nr) )
  rw(1:nr) = rw(1:nr) / ( radius2 - radius1 ) / ( radius2 + radius1 )
!
!  Set angular abscissas and weights.
!
  nt = 4 * nr

  tw = 1.0D+00 / real ( nt, kind = 8 )
!
!  Approximate the integral.
!
  quad = 0.0D+00
  do i = 1, nt
    t = 2.0D+00 * pi * real ( i - 1, kind = 8 ) / real ( nt, kind = 8 )
    do j = 1, nr
      x = center(1) + ra(j) * cos ( t )
      y = center(2) + ra(j) * sin ( t )
      quad = quad + tw * rw(j) * func ( x, y )
    end do
  end do

  area = circle_annulus_area_2d ( radius1, radius2 )
  result = quad * area

  return
end
function circle_annulus_area_2d ( radius1, radius2 )

!*****************************************************************************80
!
!! CIRCLE_ANNULUS_AREA_2D returns the area of a circular annulus in 2D.
!
!  Integration region:
!
!    RADIUS1^2 <= ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= RADIUS2^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RADIUS1, RADIUS2, the radii of the circles.
!
!    Output, real ( kind = 8 ) CIRCLE_ANNULUS_AREA_2D, the area of the annulus.
!
  implicit none

  real ( kind = 8 ) circle_annulus_area_2d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) radius1
  real ( kind = 8 ) radius2

  circle_annulus_area_2d = pi * ( radius1 + radius2 ) &
    * ( radius2 - radius1 )

  return
end
subroutine circle_annulus_sector ( func, center, radius1, radius2, theta1, &
  theta2, nr, result )

!*****************************************************************************80
!
!! CIRCLE_ANNULUS_SECTOR approximates an integral in a circular annulus sector.
!
!  Discussion:
!
!    A circular annulus sector comprises the area between two concentric
!    circles and two concentric rays.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Peirce,
!    Numerical Integration Over the Planar Annulus,
!    Journal of the Society for Industrial and Applied Mathematics,
!    Volume 5, Number 2, June 1957, pages 66-73.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function of two
!    variables which is to be integrated, of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, real ( kind = 8 ) CENTER(2), the center of the circle.
!
!    Input, real ( kind = 8 ) RADIUS1, RADIUS2, the radii of the circles.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angles defining the sector.
!    The sector is measured from THETA1 to THETA2.
!
!    Input, integer ( kind = 4 ) NR, the order of the rule.  This quantity 
!    specifies the number of distinct radii to use.  The number of angles used 
!    will be 4*NR, for a total of 4*NR*NR points.
!
!    Output, real ( kind = 8 ) RESULT, the approximation to the integral.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a
  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) circle_annulus_sector_area_2d
  real ( kind = 8 ) d
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nt
  real ( kind = 8 ) quad
  real ( kind = 8 ) ra(nr)
  real ( kind = 8 ) radius1
  real ( kind = 8 ) radius2
  real ( kind = 8 ) result
  real ( kind = 8 ) rw(nr)
  real ( kind = 8 ) ta(4*nr)
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) tw(4*nr)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Set the radial abscissas and weights.
!
  call legendre_set ( nr, ra, rw )
  a = -1.0D+00
  b = +1.0D+00
  c = radius1 * radius1
  d = radius2 * radius2
  call rule_adjust ( a, b, c, d, nr, ra, rw )
  ra(1:nr) = sqrt ( ra(1:nr) )
  rw(1:nr) = rw(1:nr) / ( radius2 - radius1 ) / ( radius2 + radius1 )
!
!  Pick angles evenly spaced between THETA1 and THETA2, but do not
!  include the endpoints, and use a half interval for the first and last.
!
  nt = 4 * nr

  call tvec_even_bracket3 ( nt, theta1, theta2, ta )
  tw(1:nt) = 1.0D+00 / real ( nt, kind = 8 )
!
!  Approximate the integral.
!
  quad = 0.0D+00
  do i = 1, nt
    do j = 1, nr
      x = center(1) + ra(j) * cos ( ta(i) )
      y = center(2) + ra(j) * sin ( ta(i) )
      quad = quad + tw(i) * rw(j) * func ( x, y )
    end do
  end do

  area = circle_annulus_sector_area_2d ( radius1, radius2, theta1, theta2 )
  result = quad * area

  return
end
function circle_annulus_sector_area_2d ( radius1, radius2, theta1, theta2 )

!*****************************************************************************80
!
!! CIRCLE_ANNULUS_SECTOR_AREA_2D: area of a circular annulus sector in 2D.
!
!  Discussion:
!
!    A circular annulus sector comprises the area between two concentric
!    circles and two concentric rays.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RADIUS1, RADIUS2, the radii of the circles.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angles of the rays.
!    Ordinarily, (THETA2-THETA1) is between 0 and 2*PI.
!
!    Output, real ( kind = 8 ) CIRCLE_ANNULUS_SECTOR_AREA_2D, the area of the
!    circular annulus sector.
!
  implicit none

  real ( kind = 8 ) circle_annulus_sector_area_2d
  real ( kind = 8 ) radius1
  real ( kind = 8 ) radius2
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  circle_annulus_sector_area_2d = 0.5D+00 * ( radius1 + radius2 ) &
    * ( radius2 - radius1 ) * ( theta2 - theta1 )

  return
end
function circle_area_2d ( r )

!*****************************************************************************80
!
!! CIRCLE_AREA_2D returns the area of a circle in 2D.
!
!  Integration region:
!
!    ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= R * R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Output, real ( kind = 8 ) CIRCLE_AREA_2D, the area of the circle.
!
  implicit none

  real ( kind = 8 ) circle_area_2d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r

  circle_area_2d = pi * r * r

  return
end
subroutine circle_cap_area_2d ( r, h, area )

!*****************************************************************************80
!
!! CIRCLE_CAP_AREA_2D computes the area of a circle cap in 2D.
!
!  Discussion:
!
!    Draw any radius R of the circle and denote as P the point where the
!    radius intersects the circle.  Now consider the point Q which lies
!    on the radius and which is H units from P.  The line which is
!    perpendicular to the radius R and passes through Q divides the
!    circle into two pieces.  The piece including the point P is the
!    circular cap of height (or thickness) H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) H, the "height" of the circle cap.  
!
!    Output, real ( kind = 8 ) AREA, the area of the circle cap.
!
  implicit none

  real ( kind = 8 ) arc_sine
  real ( kind = 8 ) area
  real ( kind = 8 ) h
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) theta

  if ( h <= 0.0D+00 ) then

    area = 0.0D+00

  else if ( h <= r ) then

    theta = 2.0D+00 * arc_sine ( sqrt ( h * ( 2.0D+00 * r - h ) ) / r )
    area = r * r * ( theta - sin ( theta ) ) / 2.0D+00

  else if ( h <= 2.0D+00 * r ) then

    theta = 2.0D+00 * arc_sine ( sqrt ( h * ( 2.0D+00 * r - h ) ) / r )
    area = r * r * ( pi - ( theta - sin ( theta ) ) / 2.0D+00 )

  else if ( 2.0D+00 * r <= h ) then

    area = pi * r * r

  end if

  return
end
subroutine circle_cum ( func, center, radius, order, result )

!*****************************************************************************80
!
!! CIRCLE_CUM approximates an integral on the circumference of a circle in 2D.
!
!  Integration region:
!
!    ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= R * R
!
!  Discussion:
!
!    An ORDER point, (ORDER-1)-th degree formula is used, 
!    Stroud number U2:M-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function of two
!    variables which is to be integrated, of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, real ( kind = 8 ) CENTER(2), the coordinates of the center of 
!    the circle.
!
!    Input, real ( kind = 8 ) RADIUS, the radius of the circle.
!
!    Input, integer ( kind = 4 ) ORDER, the number of points to use.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) radius
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  quad = 0.0D+00

  do i = 1, order
    angle = real ( 2 * i, kind = 8 ) * pi / real ( order, kind = 8 )
    x = center(1) + radius * cos ( angle )
    y = center(2) + radius * sin ( angle )
    quad = quad + func ( x, y )
  end do

  quad = quad / real ( order, kind = 8 )

  volume = pi * radius * radius
  result = quad * volume

  return
end
subroutine circle_rt_set ( rule, nr, nt, nc, ra, rw, ta, tw, cw )

!*****************************************************************************80
!
!! CIRCLE_RT_SET sets an R, THETA product quadrature rule in the unit circle.
!
!  Discussion:
!
!    For a given value of RULE, here are the number of points used at the
!    center (NC), the number of points along the radial direction (NR) and
!    the number of points along the theta direction (NT).  The total number
!    of points in the rule will be 
!
!      Total = NC + NR * NT.
!
!    The user, when choosing RULE, must allocate enough space in the arrays
!    RA, RW, TA and TW for the resulting values of NR and NT.
!
!    RULE  NC  NR  NT  Total
!    ----  --  --  --  -----
!       1   1   0   0      1
!       2   0   1   4      4
!       3   1   1   4      5
!       4   1   1   6      7
!       5   1   2   4      9
!       6   0   3   4     12
!       7   1   2  10     21
!       8   0   4  16     64
!       9   0   5  20    120
!
!    The integral of F(X,Y) over the unit circle is approximated by
!
!      Integral ( X*X + Y*Y <= 1 ) F(X,Y) dx dy 
!      = Integral ( 0 <= R <= 1, 0 <= T <= 2PI ) F(R*cos(T),R*sin(T)) r dr dt
!      = approximately
!        CW * F(0,0) 
!        + sum ( 1 <= I <= NR ) Sum ( 1 <= J <= NT )
!        RW(I) * TW(J) * F ( R(I) * cos ( TA(J) ), R(I) * sin ( TA(J) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the rule desired.
!
!    Input, integer ( kind = 4 ) NR, the number of R abscissas.
!
!    Input, integer ( kind = 4 ) NT, the number of Theta abscissas.
!
!    Input, integer ( kind = 4 ) NC, the number of center abscissas (0 or 1 ).
!
!    Output, real ( kind = 8 ) RA(NR), RW(NR), the R abscissas and weights.
!
!    Output, real ( kind = 8 ) TA(NT), TW(NT), the THETA abscissas and weights.
!
!    Output, real ( kind = 8 ) ZW, the weight to use for the center.
!
  implicit none

  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nt

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cw
  real ( kind = 8 ) d
  integer ( kind = 4 ) nc
  real ( kind = 8 ) ra(nr)
  real ( kind = 8 ) rw(nr)
  integer ( kind = 4 ) rule
  real ( kind = 8 ) ta(nt)
  real ( kind = 8 ) tw(nt)
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w


  if ( rule == 1 ) then

    cw = 1.0D+00

  else if ( rule == 2 ) then

    ra(1) = 0.5D+00
    rw(1) = 1.0D+00

    call tvec_even2 ( nt, ta )
    tw(1:nt) = 1.0D+00 / real ( nt, kind = 8 )

    cw = 0.0D+00

  else if ( rule == 3 ) then

    ra(1) = 1.0D+00
    rw(1) = 1.0D+00

    call tvec_even ( nt, ta )
    tw(1:4) = 0.125D+00

    cw = 0.5D+00

  else if ( rule == 4 ) then

    ra(1) = sqrt ( 2.0D+00 / 3.0D+00 )
    rw(1) = 1.0D+00

    call tvec_even ( nt, ta )
    tw(1:nt) = 0.125D+00

    cw = 0.25D+00

  else if ( rule == 5 ) then

    a = 1.0D+00
    b = sqrt ( 2.0D+00 ) / 2.0D+00
    u = 1.0D+00 / 6.0D+00
    v = 4.0D+00 / 6.0D+00

    ra(1:nr) = (/ a, b /)
    rw(1:nr) = (/ u, v /)

    call tvec_even ( nt, ta )
    tw(1:nt) = 1.0D+00 / real ( nt, kind = 8 )

    cw = 4.0D+00 / 24.0D+00

  else if ( rule == 6 ) then

    a = sqrt ( 3.0D+00 ) / 2.0D+00
    b = sqrt ( ( 27.0D+00 - 3.0D+00 * sqrt ( 29.0D+00 ) ) / 52.0D+00 )
    c = sqrt ( ( 27.0D+00 + 3.0D+00 * sqrt ( 29.0D+00 ) ) / 52.0D+00 )

    u = 8.0D+00 / 27.0D+00
    v = ( 551.0D+00 + 41.0D+00 * sqrt ( 29.0D+00 ) ) / 1566.0D+00
    w = ( 551.0D+00 - 41.0D+00 * sqrt ( 29.0D+00 ) ) / 1566.0D+00

    ra(1:nr) = (/ a, b, c /)
    rw(1:nr) = (/ u, v, w /)

    call tvec_even ( nt, ta )
    tw(1:nt) = 1.0D+00 / real ( nt, kind = 8 )

    cw = 0.0D+00

  else if ( rule == 7 ) then

    a = sqrt ( ( 6.0D+00 - sqrt ( 6.0D+00 ) ) / 10.0D+00 )
    b = sqrt ( ( 6.0D+00 + sqrt ( 6.0D+00 ) ) / 10.0D+00 )
    u = ( 16.0D+00 + sqrt ( 6.0D+00 ) ) / 36.0D+00
    v = ( 16.0D+00 - sqrt ( 6.0D+00 ) ) / 36.0D+00

    ra(1:nr) = (/ a, b /)
    rw(1:nr) = (/ u, v /)

    call tvec_even ( nt, ta )
    tw(1:nt) = 1.0D+00 / real ( nt, kind = 8 )

    cw = 1.0D+00 / 9.0D+00

  else if ( rule == 8 ) then

    call legendre_set ( nr, ra, rw )
    a = -1.0D+00
    b = +1.0D+00
    c =  0.0D+00
    d = +1.0D+00
    call rule_adjust ( a, b, c, d, nr, ra, rw )
    ra(1:nr) = sqrt ( ra(1:nr) )

    call tvec_even ( nt, ta )
    tw(1:nt) = 1.0D+00 / real ( nt, kind = 8 )

    cw = 0.0D+00

  else if ( rule == 9 ) then

    call legendre_set ( nr, ra, rw )
    a = -1.0D+00
    b = +1.0D+00
    c =  0.0D+00
    d = +1.0D+00
    call rule_adjust ( a, b, c, d, nr, ra, rw )
    ra(1:nr) = sqrt ( ra(1:nr) )

    call tvec_even ( nt, ta )
    tw(1:nt) = 1.0D+00 / real ( nt, kind = 8 )

    cw = 0.0D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIRCLE_RT_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  There is no rule of index ', rule
    stop

  end if

  return
end
subroutine circle_rt_size ( rule, nr, nt, nc )

!*****************************************************************************80
!
!! CIRCLE_RT_SIZE sizes an R, THETA product quadrature rule in the unit circle.
!
!  Discussion:
!
!    For a given value of RULE, here are the number of points used at the
!    center (NC), the number of points along the radial direction (NR) and
!    the number of points along the theta direction (NT).  The total number
!    of points in the rule will be 
!
!      Total = NC + NR * NT.
!
!    The user, when choosing RULE, must allocate enough space in the arrays
!    RA, RW, TA and TW for the resulting values of NR and NT.
!
!    RULE  NC  NR  NT  Total
!    ----  --  --  --  -----
!       1   1   0   0      1
!       2   0   1   4      4
!       3   1   1   4      5
!       4   1   1   6      7
!       5   1   2   4      9
!       6   0   3   4     12
!       7   1   2  10     21
!       8   0   4  16     64
!       9   0   5  20    120
!
!    The integral of F(X,Y) over the unit circle is approximated by
!
!      Integral ( X*X + Y*Y <= 1 ) F(X,Y) dx dy 
!      = Integral ( 0 <= R <= 1, 0 <= T <= 2PI ) F(R*cos(T),R*sin(T)) r dr dt
!      = approximately
!        ZW * F(0,0) 
!        + sum ( 1 <= I <= NR ) Sum ( 1 <= J <= NT )
!        RW(I) * TW(J) * F ( R(I) * cos ( TA(J) ), R(I) * sin ( TA(J) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the rule desired.
!
!    Output, integer ( kind = 4 ) NR, the number of R abscissas.
!    
!    Output, integer ( kind = 4 ) NT, the number of Theta abscissas.
!
!    Output, integer ( kind = 4 ) NC, the number of center abscissas (0 or 1).
!
  implicit none

  integer ( kind = 4 ) nc
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) rule

  if ( rule == 1 ) then

    nr = 0
    nt = 0
    nc = 1

  else if ( rule == 2 ) then

    nr = 1
    nt = 4
    nc = 0

  else if ( rule == 3 ) then

    nr = 1
    nt = 4
    nc = 1

  else if ( rule == 4 ) then

    nr = 1
    nt = 6
    nc = 1

  else if ( rule == 5 ) then

    nr = 2
    nt = 4
    nc = 1

  else if ( rule == 6 ) then

    nr = 3
    nt = 4
    nc = 0

  else if ( rule == 7 ) then

    nr = 2
    nt = 10
    nc = 1

  else if ( rule == 8 ) then

    nr = 4
    nt = 16
    nc = 0

  else if ( rule == 9 ) then

    nr = 5
    nt = 20
    nc = 0

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIRCLE_RT_SIZE - Fatal error!'
    write ( *, '(a,i8)' ) '  There is no rule of index ', rule
    stop

  end if

  return
end
subroutine circle_rt_sum ( func, center, radius, nr, ra, rw, nt, ta, tw, zw, &
  result )

!*****************************************************************************80
!
!! CIRCLE_RT_SUM applies an R, THETA product quadrature rule inside a circle.
!
!  Integration region:
!
!    (X-CENTER(1))^2 + (Y-CENTER(2))^2 <= RADIUS^2.
!
!  Discussion:
!
!    The product rule is assumed to be have the form:
!
!      Integral_Approx = ZW * F(CENTER(1),CENTER(2)) +
!        sum ( 1 <= IR <= NR ) Sum ( 1 <= IT <= NT )
!        RW(IR) * TW(IT) * F ( CENTER(1) + R(IR) * RADIUS * Cos ( TA(IT) ),
!                              CENTER(2) + R(IR) * RADIUS * Sin ( TA(IT) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of two variables which is to be integrated,
!    of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, real ( kind = 8 ) CENTER(2), the center of the circle.
!
!    Input, real ( kind = 8 ) RADIUS, the radius of the circle.
!
!    Input, integer ( kind = 4 ) NR, the number of R abscissas.
!
!    Input, real ( kind = 8 ) RA(NR), RW(NR), the R abscissas and weights.
!
!    Input, integer ( kind = 4 ) NT, the number of Theta abscissas.
!
!    Input, real ( kind = 8 ) TA(NT), TW(NT), the THETA abscissas and weights.
!
!    Input, real ( kind = 8 ) ZW, the weight to use for the center.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nt

  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) circle_area_2d
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) it
  real ( kind = 8 ) quad
  real ( kind = 8 ) ra(nr)
  real ( kind = 8 ) radius
  real ( kind = 8 ) rct
  real ( kind = 8 ) result
  real ( kind = 8 ) rst
  real ( kind = 8 ) rw(nr)
  real ( kind = 8 ) ta(nt)
  real ( kind = 8 ) tw(nt)
  real ( kind = 8 ) volume
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) zw

  quad = 0.0D+00

  if ( zw /= 0.0D+00 ) then
    x = center(1)
    y = center(2)
    quad = quad + zw * func ( x, y )
  end if

  do it = 1, nt
    rct = radius * cos ( ta(it) )
    rst = radius * sin ( ta(it) )
    do ir = 1, nr
      x = center(1) + ra(ir) * rct
      y = center(2) + ra(ir) * rst
      quad = quad + tw(it) * rw(ir) * func ( x, y )
    end do
  end do

  volume = circle_area_2d ( radius )
  result = quad * volume

  return
end
subroutine circle_sector ( func, center, radius, theta1, theta2, nr, result )

!*****************************************************************************80
!
!! CIRCLE_SECTOR approximates an integral in a circular sector.
!
!  Discussion:
!
!    A sector is contained within a circular arc and the lines joining each
!    endpoint of the arc to the center of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function of two
!    variables which is to be integrated, of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, real ( kind = 8 ) CENTER(2), the center of the circle.
!
!    Input, real ( kind = 8 ) RADIUS, the radius of the circle.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angles defining the sector.
!    The sector is measured from THETA1 to THETA2.
!
!    Input, integer ( kind = 4 ) NR, the number of radial values used in the 
!    approximation of the integral.  NR must be at least 1.  Higher values 
!    improve the accuracy of the integration, at the cost of more function 
!    evaluations.
!
!    Output, real ( kind = 8 ) RESULT, the approximation to the integral.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) nr

  real ( kind = 8 ) a
  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) circle_sector_area_2d
  real ( kind = 8 ) d
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nt
  real ( kind = 8 ) quad
  real ( kind = 8 ) ra(nr)
  real ( kind = 8 ) radius
  real ( kind = 8 ) result
  real ( kind = 8 ) rw(nr)
  real ( kind = 8 ) ta(4*nr)
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) tw(4*nr)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Set the radial abscissas and weights.
!
  call legendre_set ( nr, ra, rw )
  a = -1.0D+00
  b = +1.0D+00
  c =  0.0D+00
  d =  radius * radius
  call rule_adjust ( a, b, c, d, nr, ra, rw )
  ra(1:nr) = sqrt ( ra(1:nr) )
  rw(1:nr) = rw(1:nr) / radius / radius
!
!  Pick angles evenly spaced between THETA1 and THETA2, but do not
!  include the endpoints, and use a half interval for the first and last.
!
  nt = 4 * nr

  call tvec_even_bracket3 ( nt, theta1, theta2, ta )
  tw(1:nt) = 1.0D+00 / real ( nt, kind = 8 )
!
!  Approximate the integral.
!
  quad = 0.0D+00

  do i = 1, nr
    do j = 1, nt
      x = center(1) + ra(i) * cos ( ta(j) )
      y = center(2) + ra(i) * sin ( ta(j) )
      quad = quad + rw(i) * tw(j) * func ( x, y )
    end do
  end do

  area = circle_sector_area_2d ( radius, theta1, theta2 )
  result = quad * area

  return
end
function circle_sector_area_2d ( r, theta1, theta2 )

!*****************************************************************************80
!
!! CIRCLE_SECTOR_AREA_2D returns the area of a circular sector in 2D.
!
!  Discussion:
!
!    A sector is contained within a circular arc and the lines joining each
!    endpoint of the arc to the center of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angles of the rays
!    that delimit the sector.
!
!    Output, real ( kind = 8 ) CIRCLE_SECTOR_AREA_2D, the area of the sector.
!
  implicit none

  real ( kind = 8 ) circle_sector_area_2d
  real ( kind = 8 ) r
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  circle_sector_area_2d = 0.50D+00 * r * r * ( theta2 - theta1 )

  return
end
subroutine circle_sector_rule ( radius, theta1, theta2, nr, nt, ra, rw, ta, tw )

!*****************************************************************************80
!
!! CIRCLE_SECTOR_RULE approximates an integral in a circular sector.
!
!  Discussion:
!
!    A sector is contained within a circular arc and the lines joining each
!    endpoint of the arc to the center of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RADIUS, the radius of the circle.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angles defining the sector.
!    The sector is measured from THETA1 to THETA2.
!
!    Input, integer ( kind = 4 ) NR, the number of radial values.
!
!    Input, integer ( kind = 4 ) NT, the number of angular values.
!
!    Output, real ( kind = 8 ) RA(NR), RW(NR), the radial abscissas and weights.
!
!    Output, real ( kind = 8 ) TA(NT), TW(NT), the angular abscissas 
!    and weights.
!
  implicit none

  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nt

  real ( kind = 8 ) a
  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) ra(nr)
  real ( kind = 8 ) radius
  real ( kind = 8 ) rw(nr)
  real ( kind = 8 ) ta(nt)
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) tw(nt)
!
!  Set the radial abscissas and weights.
!
  call legendre_set ( nr, ra, rw )
  a = -1.0D+00
  b = +1.0D+00
  c =  0.0D+00
  d =  radius * radius
  call rule_adjust ( a, b, c, d, nr, ra, rw )
  ra(1:nr) = sqrt ( ra(1:nr) )
  rw(1:nr) = rw(1:nr) / radius / radius
!
!  Pick angles evenly spaced between THETA1 and THETA2, but do not
!  include the endpoints, and use a half interval for the first and last.
!
  call tvec_even_bracket3 ( nt, theta1, theta2, ta )
  tw(1:nt) = 1.0D+00 / real ( nt, kind = 8 )

  return
end
function circle_triangle_area_2d ( r, theta1, theta2 )

!*****************************************************************************80
!
!! CIRCLE_TRIANGLE_AREA_2D returns the area of a circle triangle in 2D.
!
!  Discussion:
!
!    A circle triangle is formed by drawing a circular arc, and considering
!    the triangle formed by the endpoints of the arc plus the center of
!    the circle.
!
!    The normal situation is that 0 < ( THETA2 - THETA1 ) < PI.  Outside
!    this range, the triangle can actually have NEGATIVE area.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angles of the rays that
!    delimit the arc.
!
!    Output, real ( kind = 8 ) CIRCLE_TRIANGLE_AREA_2D, the (signed) area
!    of the triangle.
!
  implicit none

  real ( kind = 8 ) circle_triangle_area_2d
  real ( kind = 8 ) r
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  circle_triangle_area_2d = 0.5D+00 * r * r * sin ( theta2 - theta1 )

  return
end
subroutine circle_xy_set ( rule, order, xtab, ytab, weight )

!*****************************************************************************80
!
!! CIRCLE_XY_SET sets an XY quadrature rule inside the unit circle in 2D.
!
!  Integration region:
!
!    X*X + Y*Y <= 1.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Frank Lether,
!    A Generalized Product Rule for the Circle,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 2, June 1971, pages 249-253.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the rule desired.
!      1, 1 point 1-st degree;
!      2, 4 point 3-rd degree, Stroud S2:3-1;
!      3, 4 point 3-rd degree, Lether #1;
!      4, 4 point 3-rd degree, Stroud S2:3-2;
!      5, 5 point 3-rd degree;
!      6, 7 point 5-th degree;
!      7, 9 point 5-th degree;
!      8, 9 point 5-th degree, Lether #2;
!      9, 12 point 7-th degree;
!     10, 16 point 7-th degree, Lether #3;
!     11, 21 point 9-th degree, Stroud S2:9-3;
!     12, 25 point 9-th degree, Lether #4 (after correcting error);
!     13, 64 point 15-th degree Gauss product rule.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the desired rule.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas of 
!    the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the ORDER weights of the rule.
!
  implicit none

  integer order

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nr
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) ra(4)
  real ( kind = 8 ) rw(4)
  integer ( kind = 4 ) rule
  real ( kind = 8 ) s
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) w4
  real ( kind = 8 ) w5
  real ( kind = 8 ) w6
  real ( kind = 8 ) w7
  real ( kind = 8 ) w8
  real ( kind = 8 ) w9
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) z

  if ( rule == 1 ) then

    xtab(1) = 0.0D+00
    ytab(1) = 0.0D+00
    weight(1) = 1.0D+00

  else if ( rule == 2 ) then

    a = 0.5D+00
    b = 0.25D+00
    z = 0.0D+00

    xtab(1:4) =   (/  a, -a,  z,  z /)
    ytab(1:4) =   (/  z,  z,  a, -a /)
    weight(1:4) = (/  b,  b,  b,  b /)

  else if ( rule == 3 ) then

    a = 0.5D+00
    b = 0.25D+00

    xtab(1:4) =   (/  a, -a, -a,  a /)
    ytab(1:4) =   (/  a,  a, -a, -a /)
    weight(1:4) = (/  b,  b,  b,  b /)

  else if ( rule == 4 ) then

    a = sqrt ( 2.0D+00 ) / 2.0D+00
    b = 0.25D+00

    xtab(1:4) =   (/  a, -a, -a,  a /)
    ytab(1:4) =   (/  a,  a, -a, -a /)
    weight(1:4) = (/  b,  b,  b,  b /)

  else if ( rule == 5 ) then

    a = 1.0D+00
    b = 0.5D+00
    c = 0.125D+00
    z = 0.0D+00

    xtab(1:5) =   (/ z, a, z, -a,  z /)
    ytab(1:5) =   (/ z, z, a,  z, -a /)
    weight(1:5) = (/ b, c, c,  c,  c /)

  else if ( rule == 6 ) then

    a = sqrt ( 2.0D+00 / 3.0D+00 )
    b = sqrt ( 1.0D+00 / 6.0D+00 )
    c = sqrt ( 2.0D+00 ) / 2.0D+00
    d = 0.125D+00
    e = 0.25D+00
    z = 0.0D+00

    xtab(1:7) =   (/ z, a, -a,  b, -b,  b, -b /)
    ytab(1:7) =   (/ z, z,  z,  c,  c, -c, -c /)
    weight(1:7) = (/ e, d,  d,  d,  d,  d,  d /)

  else if ( rule == 7 ) then

    a = 0.5D+00
    b = 1.0D+00
    c = 4.0D+00 / 24.0D+00
    d = 1.0D+00 / 24.0D+00
    z = 0.0D+00

    xtab(1:9) =   (/ z,  b, -b,  z,  z,  a, -a, -a,  a /)
    ytab(1:9) =   (/ z,  z,  z,  b, -b,  a,  a, -a, -a /)
    weight(1:9) = (/ c,  d,  d,  d,  d,  c,  c,  c,  c /)

  else if ( rule == 8 ) then

    a = sqrt ( 2.0D+00 ) / 2.0D+00
    b = sqrt ( 3.0D+00 / 5.0D+00 )
    c = sqrt ( 3.0D+00 / 10.0D+00 )

    w1 = 16.0D+00 / 72.0D+00
    w2 =  8.0D+00 / 72.0D+00
    w3 = 10.0D+00 / 72.0D+00
    w4 =  5.0D+00 / 72.0D+00

    z = 0.0D+00

    xtab(1:9) =   (/  z,   a,  -a,   z,   z,   a,   a,  -a,  -a /)
    ytab(1:9) =   (/  z,   z,   z,   b,  -b,   c,  -c,   c,  -c /)
    weight(1:9) = (/ w1,  w2,  w2,  w3,  w3,  w4,  w4,  w4,  w4 /)

  else if ( rule == 9 ) then

    a = sqrt ( 3.0D+00 ) / 2.0D+00
    b = sqrt ( ( 27.0D+00 - 3.0D+00 * sqrt ( 29.0D+00 ) ) / 104.0D+00 )
    c = sqrt ( ( 27.0D+00 + 3.0D+00 * sqrt ( 29.0D+00 ) ) / 104.0D+00 )
    u = 2.0D+00 / 27.0D+00
    v = ( 551.0D+00 + 41.0D+00 * sqrt ( 29.0D+00 ) ) / 6264.0D+00
    w = ( 551.0D+00 - 41.0D+00 * sqrt ( 29.0D+00 ) ) / 6264.0D+00
    z = 0.0D+00

    xtab(1:12) =   (/ a, -a,  z,  z,  b, -b,  b, -b,  c,  c, -c, -c /)
    ytab(1:12) =   (/ z,  z,  a, -a,  b,  b, -b, -b,  c, -c,  c, -c /)
    weight(1:12) = (/ u,  u,  u,  u,  v,  v,  v,  v,  w,  w,  w,  w /)

  else if ( rule == 10 ) then

    a = sqrt ( ( 3.0D+00 - sqrt ( 5.0D+00 ) ) / 8.0D+00 )
    b = sqrt ( ( 15.0D+00 + 3.0D+00 * sqrt ( 5.0D+00 ) &
      - 2.0D+00 * sqrt ( 30.0D+00 ) - 2.0D+00 * sqrt ( 6.0D+00 ) ) / 56.0D+00 )
    c = sqrt ( ( 15.0D+00 + 3.0D+00 * sqrt ( 5.0D+00 ) &
      + 2.0D+00 * sqrt ( 30.0D+00 ) + 2.0D+00 * sqrt ( 6.0D+00 ) ) / 56.0D+00 )
    d = sqrt ( ( 3.0D+00 + sqrt ( 5.0D+00 ) ) / 8.0D+00 )
    e = sqrt ( ( 15.0D+00 - 3.0D+00 * sqrt ( 5.0D+00 ) &
      - 2.0D+00 * sqrt ( 30.0D+00 ) + 2.0D+00 * sqrt ( 6.0D+00 ) ) / 56.0D+00 )
    f = sqrt ( ( 15.0D+00 - 3.0D+00 * sqrt ( 5.0D+00 ) &
      + 2.0D+00 * sqrt ( 30.0D+00 ) - 2.0D+00 * sqrt ( 6.0D+00 ) ) / 56.0D+00 )
    w1 = ( 90.0D+00 + 5.0D+00 * sqrt ( 30.0D+00 ) &
       + 18.0D+00 * sqrt ( 5.0D+00 ) &
       + 5.0D+00 * sqrt ( 6.0D+00 ) ) / 1440.0D+00
    w2 = ( 90.0D+00 - 5.0D+00 * sqrt ( 30.0D+00 ) &
       + 18.0D+00 * sqrt ( 5.0D+00 ) &
       - 5.0D+00 * sqrt ( 6.0D+00 ) ) / 1440.0D+00
    w3 = ( 90.0D+00 + 5.0D+00 * sqrt ( 30.0D+00 ) &
       - 18.0D+00 * sqrt ( 5.0D+00 ) &
       - 5.0D+00 * sqrt ( 6.0D+00 ) ) / 1440.0D+00
    w4 = ( 90.0D+00 - 5.0D+00 * sqrt ( 30.0D+00 ) &
       - 18.0D+00 * sqrt ( 5.0D+00 ) &
       + 5.0D+00 * sqrt ( 6.0D+00 ) ) / 1440.0D+00

    xtab(1:order) =   (/  a,  a, -a, -a,  a,  a, -a, -a,  d,  d, -d, -d, &
                           d,  d, -d, -d /)
    ytab(1:order) =   (/  b, -b,  b, -b,  c, -c,  c, -c,  e, -e,  e, -e, &
                           f, -f,  f, -f /)
    weight(1:order) = (/ w1, w1, w1, w1, w2, w2, w2, w2, w3, w3, w3, w3, &
                          w4, w4, w4, w4 /)

  else if ( rule == 11 ) then

    xtab(1) = 0.0D+00
    ytab(1) = 0.0D+00

    weight(1) = 1.0D+00 / 9.0D+00
    weight(2:11) = ( 16.0D+00 + sqrt ( 6.0D+00 ) ) / 360.0D+00
    weight(12:21) = ( 16.0D+00 - sqrt ( 6.0D+00 ) ) / 360.0D+00

    r = sqrt ( ( 6.0D+00 - sqrt ( 6.0D+00 ) ) / 10.0D+00 )

    do i = 1, 10
      a = 2.0D+00 * pi * real ( i, kind = 8 ) / 10.0D+00
      xtab(1+i) = r * cos ( a )
      ytab(1+i) = r * sin ( a )
    end do

    r = sqrt ( ( 6.0D+00 + sqrt ( 6.0D+00 ) ) / 10.0D+00 )

    do i = 1, 10
      a = 2.0D+00 * pi * real ( i, kind = 8 ) / 10.0D+00
      xtab(11+i) = r * cos ( a )
      ytab(11+i) = r * sin ( a )
    end do
!
!  There was apparently a misprint in the Lether paper.  The quantity
!  which here reads "322" was printed there as "332".
!
  else if ( rule == 12 ) then

    a = 0.5D+00
    b = sqrt ( 3.0D+00 ) / 2.0D+00
    c = sqrt ( ( 35.0D+00 + 2.0D+00 * sqrt ( 70.0D+00 ) ) / 252.0D+00 )
    d = sqrt ( ( 35.0D+00 - 2.0D+00 * sqrt ( 70.0D+00 ) ) / 252.0D+00 )
    e = sqrt ( ( 35.0D+00 + 2.0D+00 * sqrt ( 70.0D+00 ) ) / 84.0D+00 )
    f = sqrt ( ( 35.0D+00 - 2.0D+00 * sqrt ( 70.0D+00 ) ) / 84.0D+00 )
    g = sqrt ( ( 35.0D+00 + 2.0D+00 * sqrt ( 70.0D+00 ) ) / 63.0D+00 )
    h = sqrt ( ( 35.0D+00 - 2.0D+00 * sqrt ( 70.0D+00 ) ) / 63.0D+00 )

    w1 = 64.0D+00 / 675.0D+00
    w2 = 16.0D+00 / 225.0D+00
    w3 = 16.0D+00 / 675.0D+00
    w4 = ( 322.0D+00 - 13.0D+00 * sqrt ( 70.0D+00 ) ) / 21600.0D+00
    w5 = ( 322.0D+00 + 13.0D+00 * sqrt ( 70.0D+00 ) ) / 21600.0D+00
    w6 = ( 322.0D+00 - 13.0D+00 * sqrt ( 70.0D+00 ) ) / 7200.0D+00
    w7 = ( 322.0D+00 + 13.0D+00 * sqrt ( 70.0D+00 ) ) / 7200.0D+00
    w8 = ( 322.0D+00 - 13.0D+00 * sqrt ( 70.0D+00 ) ) / 5400.0D+00
    w9 = ( 322.0D+00 + 13.0D+00 * sqrt ( 70.0D+00 ) ) / 5400.0D+00
    z = 0.0D+00

    xtab(1:order) =   (/  z,  a, -a,  b, -b,  b,  b, -b, -b,  b,  b, -b, -b, &
                           a,  a, -a, -a,  a,  a, -a, -a,  z,  z,  z,  z /)
    ytab(1:order) =   (/  z,  z,  z,  z,  z,  c, -c,  c, -c,  d, -d,  d, -d, &
                           e, -e,  e, -e,  f, -f,  f, -f,  g, -g,  h, -h /)
    weight(1:order) = (/ w1, w2, w2, w3, w3, w4, w4, w4, w4, w5, w5, w5, w5, &
                          w6, w6, w6, w6, w7, w7, w7, w7, w8, w8, w9, w9 /)

  else if ( rule == 13 ) then

    nr = 4
    call legendre_set ( nr, ra, rw )
    a = -1.0D+00
    b = +1.0D+00
    c =  0.0D+00
    d = +1.0D+00
    call rule_adjust ( a, b, c, d, nr, ra, rw )
    ra(1:nr) = sqrt ( ra(1:nr) )

    i = 0

    do j = 1, 16

      c = cos ( pi * real ( j, kind = 8 ) / 8.0D+00 )
      s = sin ( pi * real ( j, kind = 8 ) / 8.0D+00 )

      do k = 1, 4

        i = i + 1
        xtab(i) = c * ra(k)
        ytab(i) = s * ra(k)
        weight(i) = rw(k) / 16.0D+00

      end do

    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIRCLE_XY_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  There is no rule of index ', rule
    stop

  end if

  return
end
subroutine circle_xy_size ( rule, order )

!*****************************************************************************80
!
!! CIRCLE_XY_SIZE sizes an XY quadrature rule inside the unit circle in 2D.
!
!  Integration region:
!
!    X*X + Y*Y <= 1.0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Frank Lether,
!    A Generalized Product Rule for the Circle,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 2, June 1971, pages 249-253.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the rule desired.
!      1, 1 point 1-st degree;
!      2, 4 point 3-rd degree, Stroud S2:3-1;
!      3, 4 point 3-rd degree, Lether #1;
!      4, 4 point 3-rd degree, Stroud S2:3-2;
!      5, 5 point 3-rd degree;
!      6, 7 point 5-th degree;
!      7, 9 point 5-th degree;
!      8, 9 point 5-th degree, Lether #2;
!      9, 12 point 7-th degree;
!     10, 16 point 7-th degree, Lether #3;
!     11, 21 point 9-th degree, Stroud S2:9-3;
!     12, 25 point 9-th degree, Lether #4 (after correcting error);
!     13, 64 point 15-th degree Gauss product rule.
!
!    Output, integer ( kind = 4 ) ORDER, the order of the desired rule.
!
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ) rule

  if ( rule == 1 ) then

    order = 1

  else if ( rule == 2 ) then

    order = 4

  else if ( rule == 3 ) then

    order = 4

  else if ( rule == 4 ) then

    order = 4

  else if ( rule == 5 ) then

    order = 5

  else if ( rule == 6 ) then

    order = 7

  else if ( rule == 7 ) then

    order = 9

  else if ( rule == 8 ) then

    order = 9

  else if ( rule == 9 ) then

    order = 12

  else if ( rule == 10 ) then

    order = 16

  else if ( rule == 11 ) then

    order = 21

  else if ( rule == 12 ) then

    order = 25

  else if ( rule == 13 ) then

    order = 64

  else

    order = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIRCLE_XY_SIZE - Fatal error!'
    write ( *, '(a,i8)' ) '  There is no rule of index ', rule
    stop

  end if

  return
end
subroutine circle_xy_sum ( func, center, r, order, xtab, ytab, weight, &
  result )

!*****************************************************************************80
!
!! CIRCLE_XY_SUM applies an XY quadrature rule inside a circle in 2D.
!
!  Integration region:
!
!    (X-CENTER(1))^2 + (Y-CENTER(2))^2 <= R * R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of two variables which is to be integrated,
!    of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, real ( kind = 8 ) CENTER(2), the coordinates of the center of 
!    the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.  The rule is
!    assumed to be defined on the unit circle.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the XY
!    coordinates of the abscissas of the quadrature rule for the unit circle.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) order

  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) circle_area_2d
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) r
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) y
  real ( kind = 8 ) ytab(order)

  quad = 0.0D+00
  do i = 1, order
    x = center(1) + r * xtab(i)
    y = center(2) + r * ytab(i)
    quad = quad + weight(i) * func ( x, y )
  end do

  volume = circle_area_2d ( r )
  result = quad * volume

  return
end
subroutine cn_geg_00_1 ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! CN_GEG_00_1 implements the midpoint rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_00_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  expon = 0
  call c1_geg_monomial_integral ( alpha, expon, volume )
  volume = volume ** n

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = volume

  return
end
subroutine cn_geg_00_1_size ( n, alpha, o )

!*****************************************************************************80
!
!! CN_GEG_00_1_SIZE sizes the midpoint rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_00_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine cn_geg_01_1 ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! CN_GEG_01_1 implements a precision 1 rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_01_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  expon = 0
  call c1_geg_monomial_integral ( alpha, expon, value1 )
  volume = value1 ** n

  expon = 1
  call c1_geg_monomial_integral ( alpha, expon, value2 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = value2 / value1
  w(k) = volume

  return
end
subroutine cn_geg_01_1_size ( n, alpha, o )

!*****************************************************************************80
!
!! CN_GEG_01_1_SIZE sizes a precision 1 rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_01_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine cn_geg_02_xiu ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! CN_GEG_02_XIU implements the Xiu precision 2 rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) coef
  real ( kind = 8 ) delta0
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_02_XIU - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = 1.0D+00
  delta0 = 0.0D+00
  c1 = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  expon = 0
  call c1_geg_monomial_integral ( alpha, expon, volume_1d )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_geg_02_xiu_size ( n, alpha, o )

!*****************************************************************************80
!
!! CN_GEG_02_XIU_SIZE sizes the Xiu rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_02_XIU_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = n + 1

  return
end
subroutine cn_geg_03_xiu ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! CN_GEG_03_XIU implements the Xiu precision 3 rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_03_XIU - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  expon = 0
  call c1_geg_monomial_integral ( alpha, expon, volume )
  volume = volume ** n

  do j = 1, o

    i = 0 
    do r = 1, n / 2
      arg = real ( ( 2 * r - 1 ) * j, kind = 8 ) * pi / real ( n, kind = 8 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) &
        / sqrt ( 2.0D+00 * alpha + 3.0D+00 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg ) &
        / sqrt ( 2.0D+00 * alpha + 3.0D+00 )
    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * r8_mop ( j ) &
        / sqrt ( 2.0D+00 * alpha + 3.0D+00 )
      if ( n == 1 ) then
        x(i,j) = x(i,j) / sqrt ( 2.0D+00 )
      end if
    end if

  end do

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_geg_03_xiu_size ( n, alpha, o )

!*****************************************************************************80
!
!! CN_GEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_GEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_03_XIU_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = 2 * n

  return
end
subroutine cn_geg_monomial_integral ( n, alpha, expon, value )

!*****************************************************************************80
!
!! CN_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on CN.
!
!  Discussion:
!
!    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
!
!      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
!
!    with -1.0 < alpha.
!
!    value = integral ( CN ) 
!      product ( 1 <= i <= n ) x(I)^expon(i) (1-x(i)^2)^alpha dx(i)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X).
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) EXPON(N), the exponents.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) value
  real ( kind = 8 ) value2

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_GEG_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  value = 1.0D+00
  do i = 1, n
    call c1_geg_monomial_integral ( alpha, expon(i), value2 )
    value = value * value2
  end do

  return
end
subroutine cn_jac_00_1 ( n, alpha, beta, o, x, w )

!*****************************************************************************80
!
!! CN_JAC_00_1 implements the midpoint rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
!
!    with -1 < alpha, -1 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) volume
  real ( kind = 8 ) :: w(o)
  real ( kind = 8 ) :: x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_00_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_00_1 - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  expon = 0
  call c1_jac_monomial_integral ( alpha, beta, expon, volume )
  volume = volume ** n

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = volume

  return
end
subroutine cn_jac_00_1_size ( n, alpha, beta, o )

!*****************************************************************************80
!
!! CN_JAC_00_1_SIZE sizes the midpoint rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
!
!    with -1 < alpha, -1 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_00_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_00_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine cn_jac_01_1 ( n, alpha, beta, o, x, w )

!*****************************************************************************80
!
!! CN_JAC_01_1 implements a precision 1 rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha. 
!
!    with -1 < alpha, -1 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2
  real ( kind = 8 ) volume
  real ( kind = 8 ) :: w(o)
  real ( kind = 8 ) :: x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_01_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_01_1 - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  expon = 0
  call c1_jac_monomial_integral ( alpha, beta, expon, value1 )
  volume = value1 ** n

  expon = 1
  call c1_jac_monomial_integral ( alpha, beta, expon, value2 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = value2 / value1
  w(k) = volume

  return
end
subroutine cn_jac_01_1_size ( n, alpha, beta, o )

!*****************************************************************************80
!
!! CN_JAC_01_1_SIZE sizes a precision 1 rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha. 
!
!    with -1 < alpha, -1 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_01_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_01_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine cn_jac_02_xiu ( n, alpha, beta, o, x, w )

!*****************************************************************************80
!
!! CN_JAC_02_XIU implements the Xiu precision 2 rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN_JAC is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
!
!    with -1.0 < alpha, -1.0 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) beta
  real ( kind = 8 ) c1
  real ( kind = 8 ) coef
  real ( kind = 8 ) delta0
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_02_XIU - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_02_XIU - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = ( alpha + beta + 2.0D+00 ) / 2.0D+00
  delta0 = ( alpha - beta ) / 2.0D+00
  c1 = 2.0D+00 * ( alpha + 1.0D+00 ) * ( beta + 1.0D+00 ) &
    / ( alpha + beta + 3.0D+00 ) / ( alpha + beta + 2.0D+00 )

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  expon = 0
  call c1_jac_monomial_integral ( alpha, beta, expon, volume_1d )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_jac_02_xiu_size ( n, alpha, beta, o )

!*****************************************************************************80
!
!! CN_JAC_02_XIU_SIZE sizes the Xiu rule for region CN_JAC.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
!
!      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
!
!    with -1 < alpha, -1 < beta.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters.
!    -1.0 < ALPHA, -1.0 < BETA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_02_XIU_SIZE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_JAC_02_XIU_SIZE - Fatal error!'
    write ( *, '(a)' ) '  BETA <= -1.0'
    stop
  end if

  o = n + 1

  return
end
subroutine cn_jac_monomial_integral ( n, alpha, beta, expon, value )

!*****************************************************************************80
!
!! CN_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over CN.
!
!  Discussion:
!
!    value = integral ( CN ) 
!      product ( 1 <= i <= n ) x(I)^expon(i) (1-x(i))^alpha (1+x(i))^beta dx(i)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X) in the weight factor.
!
!    Input, real ( kind = 8 ) BETA, the exponent of (1+X) in the weight factor.
!
!    Input, integer ( kind = 4 ) EXPON(N), the exponents.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) expon(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) value
  real ( kind = 8 ) value2

  value = 1.0D+00
  do i = 1, n
    call c1_jac_monomial_integral ( alpha, beta, expon(i), value2 )
    value = value * value2
  end do

  return
end
subroutine cn_leg_01_1 ( n, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_01_1 implements the midpoint rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) value1
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  expon = 0
  call c1_leg_monomial_integral ( expon, value1 )
  volume = value1 ** n

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  w(k) = volume

  return
end
subroutine cn_leg_01_1_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_01_1_SIZE sizes the midpoint rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 1

  return
end
subroutine cn_leg_02_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_02_XIU implements the Xiu precision 2 rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) coef
  real ( kind = 8 ) delta0
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = 1.0D+00
  delta0 = 0.0D+00
  c1 = 1.0D+00 / 3.0D+00

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  expon = 0
  call c1_leg_monomial_integral ( expon, volume_1d )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_leg_02_xiu_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_02_XIU_SIZE sizes the Xiu rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n + 1

  return
end
subroutine cn_leg_03_1 ( n, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_03_1 implements the Stroud rule CN:3-1 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!    The necessary treatment of the final coordinate of points when
!    N is odd seems to vary from what Stroud declares! 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit  none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  expon = 0
  call c1_leg_monomial_integral ( expon, volume )
  volume = volume ** n

  do j = 1, o

    i = 0

    do r = 1, ( n / 2 )
      arg = real ( ( 2 * r - 1 ) * j, kind = 8 ) * pi / real ( n, kind = 8 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) / sqrt ( 3.0D+00 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg ) / sqrt ( 3.0D+00 )
    end do
!
!  The following code does not correspond to what Stroud declares.
!
    if ( i < n ) then

      i = i + 1
      if ( n == 1 ) then
        x(i,j) =                r8_mop ( j ) / sqrt ( 3.0D+00 )
      else
        x(i,j) = sqrt ( 2.0 ) * r8_mop ( j ) / sqrt ( 3.0D+00 )
      end if
    end if

  end do

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_leg_03_1_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_03_1_SIZE sizes the Stroud rule CN:3-1 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n

  return
end
subroutine cn_leg_03_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_03_XIU implements the Xiu precision 3 rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  expon = 0
  call c1_leg_monomial_integral ( expon, volume )
  volume = volume ** n

  do j = 1, o

    i = 0 
    do r = 1, n / 2
      arg = real ( ( 2 * r - 1 ) * j, kind = 8 ) * pi / real ( n, kind = 8 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) / sqrt ( 3.0D+00 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg ) / sqrt ( 3.0D+00 )
    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * r8_mop ( j ) / sqrt ( 3.0D+00 )
      if ( n == 1 ) then
        x(i,j) = x(i,j) / sqrt ( 2.0D+00 )
      end if
    end if

  end do

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine cn_leg_03_xiu_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n

  return
end
subroutine cn_leg_05_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_05_1 implements the Stroud rule CN:5-1 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N^2 + N + 2.
!
!    The rule has precision P = 5.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N must be 4, 5, or 6.
!
!    Input, integer ( kind = 4 ) OPTION, is only used in case N = 4 or 5.
!    In that case, OPTION should be 1 or 2 to select the
!    two available variants of the rule.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) eta
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) k
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu
  integer ( kind = 4 ) option
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)
  real ( kind = 8 ) xsi

  if ( n < 4 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_LEG_05_1 - Fatal error!'
    write ( *, '(a)' ) '  The value of N must be 4, 5, or 6.'
    stop
  end if

  if ( n == 4 .or. n == 5 ) then
    if ( option < 1 .or. 2 < option ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CN_LEG_05_1 - Fatal error!'
      write ( *, '(a)' ) '  When N = 4 or 5, OPTION must be 1 or 2.'
      stop
    end if
  end if

  expon = 0
  call c1_leg_monomial_integral ( expon, volume )
  volume = volume ** n

  if ( n == 4 .and. option == 1 ) then

    eta    =   0.778984505799815D+00
    lambda =   1.284565137874656D+00
    xsi =    - 0.713647298819253D+00
    mu =     - 0.715669761974162D+00
    gamma =    0.217089151000943D+00
    a =        0.206186096875899D-01 * volume
    b =        0.975705820221664D-02 * volume
    c =        0.733921929172573D-01 * volume

  else if ( n == 4 .and. option == 2 ) then

    eta    =   0.546190755827425D+00
    lambda =   0.745069130115661D+00
    xsi =    - 0.413927294508700D+00
    mu =     - 0.343989637454535D+00
    gamma =    1.134017894600344D+00
    a =        0.853094758323323D-01 * volume
    b =        0.862099000096395D-01 * volume
    c =        0.116418206881849D-01 * volume

  else if ( n == 5 .and. option == 1 ) then

    eta    =   0.522478547481276D+00
    lambda =   0.936135175985774D+00
    xsi =    - 0.246351362101519D+00
    mu =     - 0.496308106093758D+00
    gamma =    0.827180176822930D+00
    a =        0.631976901960153D-01 * volume
    b =        0.511464127430166D-01 * volume
    c =        0.181070246088902D-01 * volume

  else if ( n == 5 .and. option == 2 ) then

    eta    =   0.798317301388741D+00
    lambda =   0.637344273885728D+00
    xsi =    - 0.455245909918377D+00
    mu =     - 1.063446229997311D+00
    gamma =    0.354482076665770D+00
    a =        0.116952384292206D-01 * volume
    b =        0.701731258612708D-01 * volume
    c =        0.137439132264426D-01 * volume

  else if ( n == 6 ) then

    eta    =   0.660225291773525D+00
    lambda =   1.064581294844754D+00
    xsi =      0.000000000000000D+00
    mu =     - 0.660225291773525D+00
    gamma =    0.660225291773525D+00
    a =        0.182742214532872D-01 * volume
    b =        0.346020761245675D-01 * volume
    c =        0.182742214532872D-01 * volume

  end if

  k = 0

  k = k + 1
  do i = 1, n
    x(i,k) = eta
  end do
  w(k) = a

  k = k + 1
  do i = 1, n
    x(i,k) = - eta
  end do
  w(k) = a

  do i1 = 1, n
    k = k + 1
    do i = 1, n
      x(i,k) = xsi
    end do
    x(i1,k) = lambda
    w(k) = b
  end do

  do i1 = 1, n
    k = k + 1
    do i = 1, n
      x(i,k) = - xsi
    end do
    x(i1,k) = - lambda
    w(k) = b
  end do

  do i1 = 1, n - 1
    do i2 = i1 + 1, n
      k = k + 1
      do i = 1, n
        x(i,k) = gamma
      end do
      x(i1,k) = mu
      x(i2,k) = mu
      w(k) = c
    end do
  end do

  do i1 = 1, n - 1
    do i2 = i1 + 1, n
      k = k + 1
      do i = 1, n
        x(i,k) = - gamma
      end do
      x(i1,k) = - mu
      x(i2,k) = - mu
      w(k) = c
    end do
  end do

  return
end
subroutine cn_leg_05_1_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_05_1_SIZE sizes the Stroud rule CN:5-1 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N^2 + N + 2.
!
!    The rule has precision P = 5.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n * n + n + 2

  return o
end
subroutine cn_leg_05_2 ( n, o, x, w )

!*****************************************************************************80
!
!! CN_LEG_05_2 implements the Stroud rule CN:5-2 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 N^2 + 1.
!
!    The rule has precision P = 5.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N must be at least 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CN_LEG_05_2 - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2.'
    stop
  end if

  expon = 0
  call c1_leg_monomial_integral ( expon, volume )
  volume = volume ** n

  b0 = real ( 25 * n * n - 115 * n + 162, kind = 8 ) * volume / 162.0D+00
  b1 = real ( 70 - 25 * n, kind = 8 ) * volume / 162.0D+00
  b2 = 25.0D+00 * volume / 324.0D+00

  r = sqrt ( 3.0D+00 / 5.0D+00 )

  k = 0

  k = k + 1
  do i = 1, n
    x(i,k) = 0.0D+00
  end do
  w(k) = b0

  do i1 = 1, n

    k = k + 1
    do i = 1, n
      x(i,k) = 0.0D+00
    end do
    x(i1,k) = + r
    w(k) = b1

    k = k + 1
    do i = 1, n
      x(i,k) = 0.0D+00
    end do
    x(i1,k) = - r
    w(k) = b1

  end do

  do i1 = 1, n - 1
    do i2 = i1 + 1, n

      k = k + 1
      do i = 1, n
        x(i,k) = 0.0D+00
      end do
      x(i1,k) = + r
      x(i2,k) = + r
      w(k) = b2

      k = k + 1
      do i = 1, n
        x(i,k) = 0.0D+00
      end do
      x(i1,k) = + r
      x(i2,k) = - r
      w(k) = b2

      k = k + 1
      do i = 1, n
        x(i,k) = 0.0D+00
      end do
      x(i1,k) = - r
      x(i2,k) = + r
      w(k) = b2

      k = k + 1
      do i = 1, n
        x(i,k) = 0.0D+00
      end do
      x(i1,k) = - r
      x(i2,k) = - r
      w(k) = b2

    end do
  end do

  return
end
subroutine cn_leg_05_2_size ( n, o )

!*****************************************************************************80
!
!! CN_LEG_05_2_SIZE sizes the Stroud rule CN:5-2 for region CN_LEG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 N^2 + 1.
!
!    The rule has precision P = 5.
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n * n + 1

  return
end
subroutine cn_leg_monomial_integral ( n, expon, value )

!*****************************************************************************80
!
!! CN_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on CN.
!
!  Discussion:
!
!    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
!
!      w(x) = 1.
!
!    value = integral ( CN ) product ( 1 <= i <= n ) x(I)^expon(i) dx(i)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) EXPON(N), the exponents.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) expon(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) value
  real ( kind = 8 ) value2

  value = 1.0D+00
  do i = 1, n
    call c1_leg_monomial_integral ( expon(i), value2 )
    value = value * value2
  end do

  return
end
subroutine cone_unit_3d ( func, result )

!*****************************************************************************80
!
!! CONE_UNIT_3D approximates an integral inside the unit cone in 3D.
!
!  Integration Region:
!
!      X*X + Y*Y <= 1 - Z  
!
!    and
!
!      0 <= Z <= 1.
!
!  Discussion:
!
!    An 48 point degree 7 formula, Stroud CN:S2:7-1, is used.
!
!    (There is a typographical error in the S2:7-1 formula for B3.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cone_volume_3d
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) r
  real ( kind = 8 ) result
  real ( kind = 8 ), save, dimension ( 4 ) :: u = &
    (/ 0.04850054945D+00, 0.2386007376D+00, &
       0.5170472951D+00,  0.7958514179D+00 /)
  real ( kind = 8 ) volume
  real ( kind = 8 ), save, dimension ( 4 ) :: w1 = &
    (/ 0.1108884156D+00,  0.1434587878D+00, &
       0.06863388717D+00, 0.01035224075D+00 /)
  real ( kind = 8 ) w2(3)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  a = sqrt ( 3.0D+00 ) / 2.0D+00
  b = sqrt ( ( 27.0D+00 - 3.0D+00 * sqrt ( 29.0D+00 ) ) / 104.0D+00 )
  c = sqrt ( ( 27.0D+00 + 3.0D+00 * sqrt ( 29.0D+00 ) ) / 104.0D+00 )
  w2(1:3) = 3.0D+00 * (/ &
    2.0D+00 / 27.0D+00, &
    ( 551.0D+00 + 4.0D+00 * sqrt ( 29.0D+00 ) ) / 6264.0D+00, &
    ( 551.0D+00 - 4.0D+00 * sqrt ( 29.0D+00 ) ) / 6264.0D+00 /)

  quad = 0.0D+00

  do i = 1, 4

    x = a * ( 1.0D+00 - u(i) )
    y = 0.0D+00
    z = u(i)
    quad = quad + w1(i) * w2(1) * func ( x, y, z )

    x = -a * ( 1.0D+00 - u(i) )
    y = 0.0D+00
    z = u(i)
    quad = quad + w1(i) * w2(1) * func ( x, y, z )

    x = 0.0D+00
    y = a * ( 1.0D+00 - u(i) )
    z = u(i)
    quad = quad + w1(i) * w2(1) * func ( x, y, z )

    x = 0.0D+00
    y = -a * ( 1.0D+00 - u(i) )
    z = u(i)
    quad = quad + w1(i) * w2(1) * func ( x, y, z )

  end do

  do i = 1, 4

    x =  b * ( 1.0D+00 - u(i) )
    y =  b * ( 1.0D+00 - u(i) )
    z =  u(i)
    quad = quad + w1(i) * w2(2) * func ( x, y, z )

    x = -b * ( 1.0D+00 - u(i) )
    y =  b * ( 1.0D+00 - u(i) )
    z =  u(i)
    quad = quad + w1(i) * w2(2) * func ( x, y, z )

    x = -b * ( 1.0D+00 - u(i) )
    y = -b * ( 1.0D+00 - u(i) )
    z =  u(i)
    quad = quad + w1(i) * w2(2) * func ( x, y, z )

    x =  b * ( 1.0D+00 - u(i) )
    y = -b * ( 1.0D+00 - u(i) )
    z =  u(i)
    quad = quad + w1(i) * w2(2) * func ( x, y, z )

    x =  c * ( 1.0D+00 - u(i) )
    y =  c * ( 1.0D+00 - u(i) )
    z =  u(i)
    quad = quad + w1(i) * w2(3) * func ( x, y, z )

    x = -c * ( 1.0D+00 - u(i) )
    y =  c * ( 1.0D+00 - u(i) )
    z =  u(i)
    quad = quad + w1(i) * w2(3) * func ( x, y, z )

    x = -c * ( 1.0D+00 - u(i) )
    y = -c * ( 1.0D+00 - u(i) )
    z =  u(i)
    quad = quad + w1(i) * w2(3) * func ( x, y, z )

    x =  c * ( 1.0D+00 - u(i) )
    y = -c * ( 1.0D+00 - u(i) )
    z =  u(i)
    quad = quad + w1(i) * w2(3) * func ( x, y, z )

  end do

  r = 1.0D+00
  h = 1.0D+00

  volume = cone_volume_3d ( r, h )
  result = quad * volume

  return
end
function cone_volume_3d ( r, h )

!*****************************************************************************80
!
!! CONE_VOLUME_3D returns the volume of a cone in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the base of the cone.
!
!    Input, real ( kind = 8 ) H, the height of the cone.
!
!    Output, real ( kind = 8 ) CONE_VOLUME_3D, the volume of the cone.
!
  implicit none

  real ( kind = 8 ) cone_volume_3d
  real ( kind = 8 ) h
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r

  cone_volume_3d = ( pi / 3.0D+00 ) * h * r * r

  return
end
subroutine cube_shell_nd ( func, n, r1, r2, result )

!*****************************************************************************80
!
!! CUBE_SHELL_ND approximates an integral inside a cubic shell in N dimensions.
!
!  Integration region:
!
!    R1 <= abs ( X(1:N) ) <= R2
!
!  Discussion:
!
!    An N*2^N point third degree formula is used, Stroud number CNSHELL:3-4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F at the N-vector X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) R1, R2, the inner and outer radii of the cubical
!    shell.  The outer cube is of side 2*R2, the inner, missing cube of side
!    2*R1.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cube_shell_volume_nd
  logical done
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmin
  real ( kind = 8 ) result
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) x(n)

  if ( r1 == r2 ) then
    result = 0.0D+00
    return
  end if

  rmax = max ( r1, r2 )
  rmin = min ( r1, r2 )    

  u = sqrt ( real ( n, kind = 8 ) * ( rmax**(n+2) - rmin**(n+2) ) &
    / ( real ( n + 2, kind = 8 ) * ( rmax**n - rmin**n ) ) )
  v = u / sqrt ( 3.0D+00 )

  quad = 0.0D+00

  do i = 1, n

    x(1:n) = v
    x(i) = u

    do

      quad = quad + func ( n, x )

      call r8vec_mirror_next ( n, x, done )

      if ( done ) then
        exit
      end if

    end do

  end do

  quad = quad / real ( n * 2**n, kind = 8 )

  volume = cube_shell_volume_nd ( n, r1, r2 )
  result = quad * volume
  
  return
end
function cube_shell_volume_nd ( n, r1, r2 )

!*****************************************************************************80
!
!! CUBE_SHELL_VOLUME_ND computes the volume of a cubic shell in ND.
!
!  Integration region:
!
!    R1 <= abs ( X(1:N) ) <= R2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) R1, R2, the inner and outer radii of the cubic
!    shell.  The outer cube is of side 2*R2, the inner, missing cube of side
!    2*R1.
!
!    Output, real ( kind = 8 ) CUBE_SHELL_VOLUME_ND, the volume of the cubic
!    shell.
!
  implicit none

  real ( kind = 8 ) cube_shell_volume_nd
  integer ( kind = 4 ) n
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2

  cube_shell_volume_nd = ( r2**n - r1**n ) * 2**n

  return
end
subroutine cube_unit_3d ( func, result )

!*****************************************************************************80
!
!! CUBE_UNIT_3D approximates an integral inside the unit cube in 3D.
!
!  Integration region:
!
!      -1 <= X <= 1,
!    and
!      -1 <= Y <= 1,
!    and
!      -1 <= Z <= 1.
!
!  Discussion:
!
!    An 8 point third degree formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  real ( kind = 8 ) cube_unit_volume_nd
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) s
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  s = 1.0D+00 / sqrt ( 3.0D+00 )
  w = 1.0D+00 / 8.0D+00

  x = s
  y = s
  z = s

  quad = w * ( &
      func (  x,  y,  z ) + func (  x,  y, -z ) &
    + func (  x, -y,  z ) + func (  x, -y, -z ) &
    + func ( -x,  y,  z ) + func ( -x,  y, -z ) &
    + func ( -x, -y,  z ) + func ( -x, -y, -z ) )

  volume = cube_unit_volume_nd ( 3 )
  result = quad * volume

  return
end
subroutine cube_unit_nd ( func, qa, qb, n, k )

!*****************************************************************************80
!
!! CUBE_UNIT_ND approximates an integral inside the unit cube in ND.
!
!  Integration region:
!
!    -1 <= X(1:N) <= 1
!
!  Discussion:
!
!    A K^N point product formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Lyness, BJJ McHugh,
!    Integration Over Multidimensional Hypercubes, 
!    A Progressive Procedure,
!    The Computer Journal,
!    Volume 6, 1963, pages 264-270.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates the function, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Output, real ( kind = 8 ) QA(K), QB(K), two sets of estimates for
!    the integral.  The QB entries are obtained from the
!    QA entries by Richardson extrapolation, and QB(K) is
!    the best estimate for the integral.
!
!    Input, integer ( kind = 4 ) N, the dimension of the cube.
!
!    Input, integer ( kind = 4 ) K, the highest order of integration, and the 
!    order of Richardson extrapolation.  K can be no greater than 10.
!
  implicit none

  integer ( kind = 4 ), parameter :: kmax = 10

  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  real ( kind = 8 ), external :: func
  real ( kind = 8 ), save, dimension ( kmax, kmax ) :: g
  integer ( kind = 4 ) i
  real ( kind = 8 ) qa(k)
  real ( kind = 8 ) qb(k)

  g(1:kmax,1:kmax) = 0.0D+00

  g( 1, 1) =  1.0D+00
  g( 2, 1) = -0.3333333333333D+00
  g( 2, 2) =  0.1333333333333D+01
  g( 3, 1) =  0.4166666666667D-01
  g( 3, 2) = -0.1066666666667D+01
  g( 3, 3) =  0.2025000000000D+01
  g( 4, 1) = -0.2777777777778D-02
  g( 4, 2) =  0.3555555555556D+00
  g( 4, 3) = -0.2603571428571D+01
  g( 4, 4) =  0.3250793650794D+01
  g( 5, 1) =  0.1157407407407D-03
  g( 5, 2) = -0.6772486772487D-01
  g( 5, 3) =  0.1464508928571D+01
  g( 5, 4) = -0.5779188712522D+01
  g( 5, 5) =  0.5382288910935D+01
  g( 6, 1) = -0.3306878306878D-05
  g( 6, 2) =  0.8465608465608D-02
  g( 6, 3) = -0.4881696428571D+00
  g( 6, 4) =  0.4623350970018D+01
  g( 6, 5) = -0.1223247479758D+02
  g( 6, 6) =  0.9088831168831D+01
  g( 7, 1) =  0.6889329805996D-07
  g( 7, 2) = -0.7524985302763D-03
  g( 7, 3) =  0.1098381696429D+00
  g( 7, 4) = -0.2241624712736D+01
  g( 7, 5) =  0.1274216124748D+02
  g( 7, 6) = -0.2516907092907D+02
  g( 7, 7) =  0.1555944865432D+02
  g( 8, 1) = -0.1093544413650D-08
  g( 8, 2) =  0.5016656868509D-04
  g( 8, 3) = -0.1797351866883D-01
  g( 8, 4) =  0.7472082375786D+00
  g( 8, 5) = -0.8168052081717D+01
  g( 8, 6) =  0.3236023405166D+02
  g( 8, 7) = -0.5082753227079D+02
  g( 8, 8) =  0.2690606541646D+02
  g( 9, 1) =  0.1366930517063D-10
  g( 9, 2) = -0.2606055516108D-05
  g( 9, 3) =  0.2246689833604D-02
  g( 9, 4) = -0.1839281815578D+00
  g( 9, 5) =  0.3646451822195D+01
  g( 9, 6) = -0.2588818724133D+02
  g( 9, 7) =  0.7782965878964D+02
  g( 9, 8) = -0.1012934227443D+03
  g( 9, 9) =  0.4688718347156D+02
  g(10, 1) = -0.1380737896023D-12
  g(10, 2) =  0.1085856465045D-06
  g(10, 3) = -0.2222000934334D-03
  g(10, 4) =  0.3503393934435D-01
  g(10, 5) = -0.1215483940732D+01
  g(10, 6) =  0.1456210532325D+02
  g(10, 7) = -0.7477751530769D+02
  g(10, 8) =  0.1800771959898D+03
  g(10, 9) = -0.1998874663788D+03
  g(10,10) =  0.8220635246624D+02

  if ( kmax < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBE_UNIT_ND - Fatal error!'
    write ( *, '(a,i8)' ) '  K must be no greater than KMAX = ', kmax
    write ( *, '(a,i8)' ) '  but the input K is ', k
    stop
  end if

  do i = 1, k
    call qmdpt ( func, n, i, qa(i) )
  end do

  qb(1) = qa(1)

  do i = 2, k
    qb(i) = dot_product ( g(i,1:i), qa(1:i) )
  end do

  return
end
function cube_unit_volume_nd ( n )

!*****************************************************************************80
!
!! CUBE_UNIT_VOLUME_ND returns the volume of the unit cube in ND.
!
!  Integration region:
!
!    -1 <= X(1:N) <= 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) CUBE_UNIT_VOLUME_ND, the volume of the unit
!    cube in ND.
!
  implicit none

  real ( kind = 8 ) cube_unit_volume_nd
  integer ( kind = 4 ) n

  cube_unit_volume_nd = 2.0D+00**n

  return
end
function ellipse_area_2d ( r1, r2 )

!*****************************************************************************80
!
!! ELLIPSE_AREA_2D returns the area of an ellipse in 2D.
!
!  Integration region:
!
!    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, the major and minor semi-axes.
!
!    Output, real ( kind = 8 ) ELLIPSE_AREA_2D, the area of the ellipse.
!
  implicit none

  real ( kind = 8 ) ellipse_area_2d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2

  ellipse_area_2d = pi * r1 * r2

  return
end
function ellipse_circumference_2d ( r1, r2 )

!*****************************************************************************80
!
!! ELLIPSE_CIRCUMFERENCE_2D returns the circumference of an ellipse in 2D.
!
!  Discussion:
!
!    There is no closed formula for the circumference of an ellipse.
!
!    Defining the eccentricity by
!
!      E = sqrt ( 1 - ( r2 / r1 )^2 )
!
!    where R1 and R2 are the major and minor axes, then
!
!      circumference
!        = 4 * R1 * E(K,2*PI)
!        = R1 * Integral ( 0 <= T <= 2*PI ) sqrt ( 1 - E * E * sin^2 ( T ) ) dT
!
!    This integral can be approximated by the Gauss-Kummer formula.
!
!  Integration region:
!
!    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Harris, Horst Stocker,
!    Handbook of Mathematics and Computational Science,
!    Springer, 1998,
!    ISBN: 0-387-94746-9,
!    LC: QA40.S76.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, the major and minor semi-axes.
!
!    Output, real ( kind = 8 ) ELLIPSE_CIRCUMFERENCE_2D, the
!    circumference of the ellipse.
!
  implicit none

  real ( kind = 8 ) ellipse_circumference_2d
  real ( kind = 8 ) e
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) term
  real ( kind = 8 ) value

  if ( r1 == r2 ) then
    ellipse_circumference_2d = 2.0D+00 * pi * r1
    return
  end if
!
!  Compute the eccentricity of the ellipse.
!
  e = sqrt ( 1.0D+00 - ( min ( r1, r2 ) / max ( r1, r2 ) )**2 )

  value = 1.0D+00
  term = value
  i = 0

  do

    i = i + 1
    term = term * ( 2 * i - 3 ) * ( 2 * i - 1 ) * e * e &
      / real ( 2 * 2 * i * i, kind = 8 )

    if ( abs ( term ) <= epsilon ( value ) * ( abs ( value ) + 1.0D+00 ) ) then
      exit
    end if

    value = value + term

  end do

  ellipse_circumference_2d = 2.0D+00 * pi * max ( r1, r2 ) * value

  return
end
function ellipse_eccentricity_2d ( r1, r2 )

!*****************************************************************************80
!
!! ELLIPSE_ECCENTRICITY_2D returns the eccentricity of an ellipse in 2D.
!
!  Integration region:
!
!    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, the major and minor semi-axes.
!
!    Output, real ( kind = 8 ) ELLIPSE_ECCENTRICITY_2D, the eccentricity 
!    of the ellipse.
!
  implicit none

  real ( kind = 8 ) ellipse_eccentricity_2d
  real ( kind = 8 ) major
  real ( kind = 8 ) minor
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2

  minor = min ( abs ( r1 ), abs ( r2 ) )
  major = max ( abs ( r1 ), abs ( r2 ) )

  if ( major == 0.0D+00 ) then
    ellipse_eccentricity_2d = - huge ( r1 )
    return
  end if

  ellipse_eccentricity_2d = sqrt ( 1.0D+00 - ( minor / major )**2 )

  return
end
function ellipsoid_volume_3d ( r1, r2, r3 )

!*****************************************************************************80
!
!! ELLIPSOID_VOLUME_3D returns the volume of an ellipsoid in 3d.
!
!  Discussion:
!
!    This is not a general ellipsoid, but one for which each of the 
!    axes lies along a coordinate axis.
!
!  Integration region:
!
!      ( ( X - CENTER(1) ) / R1 )^2 
!    + ( ( Y - CENTER(2) ) / R2 )^2
!    + ( ( Z - CENTER(3) ) / R3 )^2 <= 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, R3, the semi-axes of the ellipsoid.
!
!    Output, real ( kind = 8 ) ELLIPSOID_VOLUME_3D, the volume of the ellipsoid.
!
  implicit none

  real ( kind = 8 ) ellipsoid_volume_3d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3

  ellipsoid_volume_3d = ( 4.0D+00 / 3.0D+00 ) * pi * r1 * r2 * r3

  return
end
subroutine en_r2_01_1 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_01_1 implements the Stroud rule 1.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume = sqrt ( pi**n )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  w(k) = volume

  return
end
subroutine en_r2_01_1_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_01_1_SIZE sizes the Stroud rule 1.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 1

  return
end
subroutine en_r2_02_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_02_XIU implements the Xiu precision 2 rule for region EN_R2.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) delta0
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = 2.0D+00
  delta0 = 0.0D+00
  c1 = 1.0D+00

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine en_r2_02_xiu_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_02_XIU_SIZE sizes the Xiu rule for region EN_R2.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n + 1

  return
end
subroutine en_r2_03_1 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_03_1 implements the Stroud rule 3.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume = sqrt ( pi**n )

  a = volume / real ( o, kind = 8 )
  r = sqrt ( real ( n, kind = 8 ) / 2.0D+00 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = a
    k = k + 1
    x(i,k) = + r
    w(k) = a
  end do

  return
end
subroutine en_r2_03_1_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_03_1_SIZE sizes the Stroud rule 3.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n

  return
end
subroutine en_r2_03_2 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_03_2 implements the Stroud rule 3.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  logical more
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume = sqrt ( pi ** n )

  a = volume / real ( o, kind = 8 )
  r = sqrt ( 0.5D+00 );

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - r
  w(k) = a
  more = .true.

  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < 0.0D+00 ) then
        k = k + 1
        x(1:i-1,k) = x(1:i-1,k-1)
        x(i,k)     = + r
        x(i+1:n,k) = - r
        w(k) = a
        more = .true.
        exit
      end if
    end do
  end do

  return
end
subroutine en_r2_03_2_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_03_2_SIZE sizes the Stroud rule 3.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 ** n;

  return
end
subroutine en_r2_03_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_03_XIU implements the Xiu precision 3 rule for region EN_R2.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  real ( kind = 8 ) arg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) :: w(o)
  real ( kind = 8 ) :: x(n,o)

  volume = sqrt ( pi ** n )

  do j = 1, o

    i = 0 
    do r = 1, n / 2
      arg = real ( ( 2 * r - 1 ) * j, kind = 8 ) * pi / real ( n, kind = 8 )
      i = i + 1
      x(i,j) = cos ( arg )
      i = i + 1
      x(i,j) = sin ( arg )
    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j )
      if ( n == 1 ) then
        x(i,j) = x(i,j) / sqrt ( 2.0D+00 )
      end if
    end if

  end do

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine en_r2_03_xiu_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_03_XIU_SIZE sizes the Xiu precision 3 rule for region EN_R2.
!
!  Discussion:
!
!    The rule has order 
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n

  return
end
subroutine en_r2_05_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_1 implements the Stroud rule 5.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = N^2 + N + 2.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting 
!    the OPTION variable to 1 or 2.
!
!    Versions of this rule are only available for N = 2 through 7.
!
!    There is a typographical error in the reference.
!    For the second version of the rule for N = 2, the line
!      gamma =    0.313300683022281D+00
!    should read
!      gamma =    0.312200683022281D+00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    2 <= N <= 7.
!
!    Input, integer ( kind = 4 ) OPTION, selects option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) eta
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) lambda
  real ( kind = 8 ) mu
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)
  real ( kind = 8 ) xsi

  if ( n < 2 .or. 7 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_1 - Fatal error!'
    write ( *, '(a)' ) '  2 <= N <= 7 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 5 .and. n /= 6 ) then 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_05_1 - Fatal error!'
      write ( *, '(a)' ) '  OPTION = 2 requires N = 3, 5 or 6.'
      stop
    end if
  end if

  volume = sqrt ( pi ** n )

  if ( n == 2 ) then
    eta =      0.446103183094540D+00
    lambda =   0.136602540378444D+01
    xsi =    - 0.366025403784439D+00
    mu =       0.198167882945871D+01
    gamma =    0.000000000000000D+00
    a =        0.328774019778636D+00 * volume
    b =        0.833333333333333D-01 * volume
    c =        0.455931355469736D-02 * volume
  else if ( n == 3 .and. option == 1 ) then
    eta =      0.476731294622796D+00
    lambda =   0.935429018879534D+00
    xsi =    - 0.731237647787132D+00
    mu =       0.433155309477649D+00
    gamma =    0.266922328697744D+01
    a =        0.242000000000000D+00 * volume
    b =        0.810000000000000D-01 * volume
    c =        0.500000000000000D-02 * volume
!
!  The value of gamma that follows corrects an error in the reference.
!
  else if ( n == 3 .and. option == 2 ) then
    eta =      0.476731294622796D+00
    lambda =   0.128679320334269D+01
    xsi =    - 0.379873463323979D+00
    mu =     - 0.192386729447751D+01
    gamma =    0.312200683022281D+00
    a =        0.242000000000000D+00 * volume
    b =        0.810000000000000D-01 * volume
    c =        0.500000000000000D-02 * volume
  else if ( n == 4 ) then
    eta =      0.523945658287507D+00
    lambda =   0.119433782552719D+01
    xsi =    - 0.398112608509063D+00
    mu =     - 0.318569372920112D+00
    gamma =    0.185675837424096D+01
    a =        0.155502116982037D+00 * volume
    b =        0.777510584910183D-01 * volume
    c =        0.558227484231506D-02 * volume
  else if ( n == 5 .and. option == 1 ) then
    eta =      0.214972564378798D+01
    lambda =   0.464252986016289D+01
    xsi =    - 0.623201054093728D+00
    mu =     - 0.447108700673434D+00
    gamma =    0.812171426076311D+00
    a =        0.487749259189752D-03 * volume
    b =        0.487749259189752D-03 * volume
    c =        0.497073504444862D-01 * volume
  else if ( n == 5 .and. option == 2 ) then
    eta =      0.615369528365158D+00
    lambda =   0.132894698387445D+01
    xsi =    - 0.178394363877324D+00
    mu =     - 0.745963266507289D+00
    gamma =    0.135503972310817D+01
    a =        0.726415024414905D-01 * volume
    b =        0.726415024414905D-01 * volume
    c =        0.641509853510569D-02 * volume
  else if ( n == 6 .and. option == 1 ) then
    eta =      0.100000000000000D+01
    lambda =   0.141421356237309D+01
    xsi =      0.000000000000000D+00
    mu =     - 0.100000000000000D+01
    gamma =    0.100000000000000D+01
    a =        0.781250000000000D-02 * volume
    b =        0.625000000000000D-01 * volume
    c =        0.781250000000000D-02 * volume
  else if ( n == 6 .and. option == 2 ) then
    eta =      0.100000000000000D+01
    lambda =   0.942809041582063D+00
    xsi =    - 0.471404520791032D+00
    mu =     - 0.166666666666667D+01
    gamma =    0.333333333333333D+00
    a =        0.781250000000000D-02 * volume
    b =        0.625000000000000D-01 * volume
    c =        0.781250000000000D-02 * volume
  else if ( n == 7 ) then
    eta =      0.000000000000000D+00
    lambda =   0.959724318748357D+00
    xsi =    - 0.772326488820521D+00
    mu =     - 0.141214270131942D+01
    gamma =    0.319908106249452D+00
    a =        0.111111111111111D+00 * volume
    b =        0.138888888888889D-01 * volume
    c =        0.138888888888889D-01 * volume
  end if

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2 points.
!
  k = k + 1
  x(1:n,k) = - eta
  w(k) = a
  k = k + 1
  x(1:n,k) = + eta
  w(k) = a
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(1:n,k) = - xsi
    x(i,k) = - lambda
    w(k) = b
    k = k + 1
    x(1:n,k) = + xsi
    x(i,k) = + lambda
    w(k) = b
  end do
!
!  2 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(1:n,k) = - gamma
      x(i,k) = - mu
      x(j,k) = - mu
      w(k) = c
      k = k + 1
      x(1:n,k) = + gamma
      x(i,k) = + mu
      x(j,k) = + mu
      w(k) = c
    end do
  end do

  return
end
subroutine en_r2_05_1_size ( n, option, o )

!*****************************************************************************80
!
!! EN_R2_05_1_SIZE sizes the Stroud rule 5.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = N^2 + N + 2.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting 
!    the OPTION variable to 1 or 2.
!
!    Versions of this rule are only available for N = 2 through 7.
!
!    There is a typographical error in the reference.
!    For the second version of the rule for N = 2, the line
!      gamma =    0.313300683022281D+00
!    should read
!      gamma =    0.312200683022281D+00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    2 <= N <= 7.
!
!    Input, integer ( kind = 4 ) OPTION, selects option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  if ( n < 2 .or. 7 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  2 <= N <= 7 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 5 .and. n /= 6 ) then 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_05_1_SIZE - Fatal error!'
      write ( *, '(a)' ) '  OPTION = 2 requires N = 3, 5 or 6.'
      stop
    end if
  end if

  o = n * n + n + 2

  return
end
subroutine en_r2_05_2 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_2 implements the Stroud rule 5.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2 * N^2 + 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume = sqrt ( pi ** n )

  a = 2.0D+00 * volume / real ( n + 2, kind = 8 )
  b = real ( 4 - n, kind = 8 ) * volume / 2.0D+00 &
    / real ( ( n + 2 ) * ( n + 2 ), kind = 8 )
  c = volume / real ( ( n + 2 ) * ( n + 2 ), kind = 8 )

  r = sqrt ( real ( n + 2, kind = 8 ) / 2.0D+00 )
  s = sqrt ( real ( n + 2, kind = 8 ) / 4.0D+00 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = a
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = b
    k = k + 1
    x(i,k) = + r
    w(k) = b
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - s
      x(j,k) = - s
      w(k) = c
      k = k + 1
      x(i,k) = - s
      x(j,k) = + s
      w(k) = c
      k = k + 1
      x(i,k) = + s
      x(j,k) = - s
      w(k) = c
      k = k + 1
      x(i,k) = + s
      x(j,k) = + s
      w(k) = c
    end do
  end do

  return
end
subroutine en_r2_05_2_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_05_2_SIZE sizes the Stroud rule 5.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2 * N^2 + 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n * n + 1

  return
end
subroutine en_r2_05_3 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_3 implements the Stroud rule 5.3 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 3 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  logical more
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_3 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  a = 4.0D+00 * volume / real ( ( n + 2 ) * ( n + 2 ), kind = 8 )
  b = real ( ( n - 2 ) * ( n - 2 ), kind = 8 ) * volume &
    / real ( 2**n, kind = 8 ) / real ( ( n + 2 ) * ( n + 2 ), kind = 8 )
  r = sqrt ( real ( n + 2, kind = 8 ) / 4.0D+00 )
  s = sqrt ( real ( n + 2, kind = 8 ) / 2.0D+00 / real ( n - 2, kind = 8 ) )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = a
    k = k + 1
    x(i,k) = + r
    w(k) = a
  end do
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - s
  w(k) = b
  more = .true.
  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < 0.0D+00 ) then
        k = k + 1
        x(1:i-1,k) = x(1:i-1,k-1)
        x(i,k)     = + s
        x(i+1:n,k) = - s
        w(k) = b
        more = .true.
        exit
      end if
    end do
  end do

  return
end
subroutine en_r2_05_3_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_05_3_SIZE sizes the Stroud rule 5.3 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 3 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_3_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    stop
  end if

  o = 2**n + 2 * n

  return
end
subroutine en_r2_05_4 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_4 implements the Stroud rule 5.4 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) - 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume = sqrt ( pi ** n )

  s = sqrt ( 0.5D+00 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2^N + 2^(N-1) + 2^(N-2) + ... + 1 = 2^(N+1)-1 points.
!  but do the last point separately.
!
  do i = 1, n

    r = sqrt ( real ( i + 2, kind = 8 ) / 2.0D+00 )
    b = 2.0D+00 ** ( i - n ) * volume / real ( i + 1, kind = 8 ) &
      / real ( i + 2, kind = 8 )

    k = k + 1
    x(i,k) = - r
    x(i+1:n,k) = - s
    w(k) = b
    more = .true.
    do while ( more )
      more = .false.
      do j = n, i, -1
        if ( x(j,k) < 0.0D+00 ) then
          k = k + 1
          x(1:j-1,k) = x(1:j-1,k-1)
          x(j,k)     =   abs ( x(j,k) )
          x(j+1:n,k) = - abs ( x(j+1:n,k) )
          w(k) = b
          more = .true.
          exit
        end if
      end do
    end do

  end do
!
!  Last point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = 2.0D+00 * volume / real ( n + 2, kind = 8 )

  return
end
subroutine en_r2_05_4_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_05_4_SIZE sizes the Stroud rule 5.4 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) - 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 ** ( n + 1 ) - 1

  return
end
subroutine en_r2_05_5 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_5 implements the Stroud rule 5.5 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = N * 2^N + 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There is a second version of this rule however it results in
!    complex abscissas, and so it has been disabled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  real ( kind = 8 ) n_r8
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  volume = sqrt ( pi ** n )

  n_r8 = real ( n, kind = 8 )

  a = 2.0D+00 * volume / ( n_r8 + 2.0D+00 )
  b =           volume / ( n_r8 + 2.0D+00 ) / ( 2.0D+00 ** n )

  option = 1

  if ( option == 1 ) then
    r = sqrt ( ( n_r8 + 2.0D+00 &
      + ( n_r8 - 1.0D+00 ) * sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
      / 2.0D+00 / n_r8 )
    s = sqrt ( ( n_r8 + 2.0D+00 &
      -                      sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
      / 2.0D+00 / n_r8 )
  else if ( option == 2 ) then
    r = sqrt ( ( n_r8 + 2.0D+00 &
      - ( n_r8 - 1.0D+00 ) * sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
      / 2.0D+00 / n_r8 )
    s = sqrt ( ( n_r8 + 2.0D+00 &
      +                      sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
      / 2.0D+00 / n_r8 )
  end if

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = a
!
!  N * 2^N points:
!  N choices for location of R, 2^N choices of sign pattern.
!
  do i = 1, n

    k = k + 1
    x(1:n,k) = - s
    x(i,k)   = - r
    w(k) = b

    more = .true.

    do while ( more )
      more = .false.
      do j = n, 1, -1
        if ( x(j,k) < 0.0D+00 ) then
          k = k + 1
          x(1:j-1,k) = x(1:j-1,k-1)
          x(j,k)     =   abs ( x(j,k) )
          x(j+1:n,k) = - abs ( x(j+1:n,k) )
          w(k) = b
          more = .true.
          exit
        end if
      end do
    end do

  end do

  return
end
subroutine en_r2_05_5_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_05_5_SIZE sizes the Stroud rule 5.5 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = N * 2^N + 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There is a second version of this rule however it results in
!    complex abscissas, and so it has been disabled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n * 2 ** n + 1

  return
end
subroutine en_r2_05_6 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_6 implements the Stroud rule 5.6 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( N + 1 ) * 2^N.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 5 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    5 <= N.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  real ( kind = 8 ) n_r8
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_6 - Fatal error!'
    write ( *, '(a)' ) '  5 <= N is required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  n_r8 = real ( n, kind = 8 )

  a = volume / ( 2.0D+00 ** n ) / ( n_r8 + 1.0D+00 )

  r = sqrt ( ( n_r8 - sqrt ( 2.0D+00 ) &
    + ( n_r8 - 1.0D+00 ) * sqrt ( 2.0D+00 * ( n_r8 + 1.0D+00 ) ) ) &
    / 2.0D+00 / n_r8 )
  s = sqrt ( ( n_r8 - sqrt ( 2.0D+00 ) &
    -                      sqrt ( 2.0D+00 * ( n_r8 + 1.0D+00 ) ) ) &
    / 2.0D+00 / n_r8 )
  t = sqrt ( ( 1.0D+00 + sqrt ( 2.0D+00 ) ) / 2.0D+00 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  N * 2^N points.
!
  do i = 1, n

    k = k + 1
    x(1:n,k) = - s
    x(i,k)   = - r
    w(k) = a

    more = .true.

    do while ( more )
      more = .false.
      do j = n, 1, -1
        if ( x(j,k) < 0.0D+00 ) then
          k = k + 1
          x(1:j-1,k) = x(1:j-1,k-1)
          x(j,k)     =   abs ( x(j,k) )
          x(j+1:n,k) = - abs ( x(j+1:n,k) )
          w(k) = a
          more = .true.
          exit
        end if
      end do
    end do

  end do
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - t
  w(k) = a
  more = .true.
  do while ( more )
    more = .false.
    do j = n, 1, -1
      if ( x(j,k) < 0.0D+00 ) then
        k = k + 1
        x(1:j-1,k) = x(1:j-1,k-1)
        x(j,k)     =   abs ( x(j,k) )
        x(j+1:n,k) = - abs ( x(j+1:n,k) )
        w(k) = a
        more = .true.
        exit
      end if
    end do
  end do

  return
end
subroutine en_r2_05_6_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_05_6_SIZE sizes the Stroud rule 5.6 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( N + 1 ) * 2^N.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 5 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    5 <= N.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( n < 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_6_SIZE - Fatal error!'
    write ( *, '(a)' ) '  5 <= N is required.'
    stop
  end if

  o = ( 2 ** n ) * ( n + 1 )

  return
end
subroutine en_r2_07_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_R2_07_1 implements the Stroud rule 7.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N^2 + 1.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of the rule, chosen by setting the
!    OPTION variable to 1 or 2.  
!
!    Option 1 is only valid for N = 3, 4, 6 or 7.
!    Option 2 is only valid for N = 3 or 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N = 3, 4, 6 or 7.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  real ( kind = 8 ) n_r8
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  if ( option == 1 ) then
    if ( n /= 3 .and. n /= 4 .and. n /= 6 .and. n /= 7 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_07_1 - Fatal error!'
      write ( *, '(a)' ) '  OPTION 1 requires N =  3, 4, 6 or 7.'
      stop
    end if
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_07_1 - Fatal error!'
      write ( *, '(a)' ) '  OPTION 2 requires N =  3 or 4.'
      stop
    end if
  end if

  volume = sqrt ( pi ** n )

  n_r8 = real ( n, kind = 8 )

  if ( option == 1 ) then
    r = sqrt ( ( 3.0D+00 * ( 8.0D+00 - n_r8 ) - ( n_r8 - 2.0D+00 ) &
      * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 / ( 5.0D+00 - n_r8 ) )
    s = sqrt ( ( 3.0D+00 *             n_r8   -          2.0D+00   &
      * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 &
      / ( 3.0D+00 * n_r8 - 8.0D+00 ) )
    t = sqrt ( ( 6.0D+00 + sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 )
  else if ( option == 2 ) then
    r = sqrt ( ( 3.0D+00 * ( 8.0D+00 - n_r8 ) + ( n_r8 - 2.0D+00 ) &
      * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 / ( 5.0D+00 - n_r8 ) )
    s = sqrt ( ( 3.0D+00 *             n_r8   +          2.0D+00   &
      * sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 &
      / ( 3.0D+00 * n_r8 - 8.0D+00 ) )
    t = sqrt ( ( 6.0D+00 - sqrt ( 3.0D+00 * ( 8.0D+00 - n_r8 ) ) ) / 2.0D+00 )
  end if

  b = ( 8.0D+00 - n_r8 ) * volume / 8.0D+00 / r ** 6
  c = volume / 2.0D+00 ** ( n + 3 ) / s ** 6
  d = volume / 16.0D+00 / t ** 6
  a = volume - 2.0D+00 * n_r8 * b - 2.0D+00 ** n * c - 2.0D+00 * n_r8 &
    * ( n_r8 - 1.0D+00 ) * d

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = a
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = b
    k = k + 1
    x(i,k) = + r
    w(k) = b
  end do
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - s
  w(k) = c
  more = .true.
  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < 0.0D+00 ) then
        k = k + 1
        x(1:i-1,k) = x(1:i-1,k-1)
        x(i,k)     =   abs ( x(i,k) )
        x(i+1:n,k) = - abs ( x(i+1:n,k) )
        w(k) = c
        more = .true.
        exit
      end if
    end do
  end do
!
!  2 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - t
      x(j,k) = - t
      w(k) = d
      k = k + 1
      x(i,k) = - t
      x(j,k) = + t
      w(k) = d
      k = k + 1
      x(i,k) = + t
      x(j,k) = - t
      w(k) = d
      k = k + 1
      x(i,k) = + t
      x(j,k) = + t
      w(k) = d
    end do
  end do

  return
end
subroutine en_r2_07_1_size ( n, option, o )

!*****************************************************************************80
!
!! EN_R2_07_1_SIZE sizes the Stroud rule 7.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N^2 + 1.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of the rule, chosen by setting the
!    OPTION variable to 1 or 2.  
!
!    Option 1 is only valid for N = 3, 4, 6 or 7.
!    Option 2 is only valid for N = 3 or 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N = 3, 4, 6 or 7.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  if ( option == 1 ) then
    if ( n /= 3 .and. n /= 4 .and. n /= 6 .and. n /= 7 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_07_1_SIZE - Fatal error!'
      write ( *, '(a)' ) '  OPTION 1 requires N =  3, 4, 6 or 7.'
      stop
    end if
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_07_1_SIZE - Fatal error!'
      write ( *, '(a)' ) '  OPTION 2 requires N =  3 or 4.'
      stop
    end if
  end if

  o = 2 ** n + 2 * n ** 2 + 1

  return
end
subroutine en_r2_07_2 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_07_2 implements the Stroud rule 7.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) + 4 * N^2.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 3 <= N.
!
!    The reference has a typographical error in the description of this rule.
!    The formula:
!
!      (t,t,t,...,t)FS
!
!    should read
!
!      (t,t,0,...,0)FS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  real ( kind = 8 ) n_r8
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) rho1
  real ( kind = 8 ) rho2
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_2 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  n_r8 = real ( n, kind = 8 )

  rho1 = sqrt ( ( n_r8 + 2.0D+00 - sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
    / 2.0D+00 )
  rho2 = sqrt ( ( n_r8 + 2.0D+00 + sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) &
    / 2.0D+00 )
  a1 = ( n_r8 + 2.0D+00 + sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) / 2.0D+00 &
    / ( n_r8 + 2.0D+00 )
  a2 = ( n_r8 + 2.0D+00 - sqrt ( 2.0D+00 * ( n_r8 + 2.0D+00 ) ) ) / 2.0D+00 &
    / ( n_r8 + 2.0D+00 )

  r = 1.0D+00
  s = sqrt ( 1.0D+00 / n_r8 )
  t = sqrt ( 0.5D+00 )
  b = ( 8.0D+00 - n_r8 ) * volume / n_r8 / ( n_r8 + 2.0D+00 ) &
    / ( n_r8 + 4.0D+00 )
  c = n_r8 ** 3 * volume / 2.0D+00 ** n / n_r8 / ( n_r8 + 2.0D+00 ) &
    / ( n_r8 + 4.0D+00 )
  d = 4.0D+00 * volume / n_r8 / ( n_r8 + 2.0D+00 ) / ( n_r8 + 4.0D+00 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  2 * 2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - rho1 * r
    w(k) = a1 * b
    k = k + 1
    x(i,k) = - rho2 * r
    w(k) = a2 * b
    k = k + 1
    x(i,k) = + rho1 * r
    w(k) = a1 * b
    k = k + 1
    x(i,k) = + rho2 * r
    w(k) = a2 * b
  end do
!
!  2 * 2^N points.
!
  k = k + 1
  x(1:n,k) = - rho1 * s
  w(k) = a1 * c
  k = k + 1
  x(1:n,k) = - rho2 * s
  w(k) = a2 * c
  more = .true.
  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < 0.0D+00 ) then
        k = k + 1
        x(1:i-1,k) =     x(1:i-1,k-2)
        x(i,k)     =   abs ( x(i,k) )
        x(i+1:n,k) = - abs ( x(i+1:n,k) )
        w(k) = a1 * c
        k = k + 1
        x(1:i-1,k) =     x(1:i-1,k-2)
        x(i,k)     =   abs ( x(i,k) )
        x(i+1:n,k) = - abs ( x(i+1:n,k) )
        w(k) = a2 * c
        more = .true.
        exit
      end if
    end do
  end do
!
!  2 * 4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - rho1 * t
      x(j,k) = - rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = - rho1 * t
      x(j,k) = + rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = + rho1 * t
      x(j,k) = - rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = + rho1 * t
      x(j,k) = + rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = - rho2 * t
      x(j,k) = - rho2 * t
      w(k) = a2 * d
      k = k + 1
      x(i,k) = - rho2 * t
      x(j,k) = + rho2 * t
      w(k) = a2 * d
      k = k + 1
      x(i,k) = + rho2 * t
      x(j,k) = - rho2 * t
      w(k) = a2 * d
      k = k + 1
      x(i,k) = + rho2 * t
      x(j,k) = + rho2 * t
      w(k) = a2 * d
    end do
  end do

  return
end
subroutine en_r2_07_2_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_07_2_SIZE sizes the Stroud rule 7.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) + 4 * N^2.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The rule requires 3 <= N.
!
!    The reference has a typographical error in the description of this rule.
!    The formula:
!
!      (t,t,t,...,t)FS
!
!    should read
!
!      (t,t,0,...,0)FS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_2_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    stop
  end if

  o = 2 ** ( n + 1 ) + 4 * n * n

  return
end
subroutine en_r2_07_3 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_R2_07_3 implements the Stroud rule 7.3 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   45
!     4   97
!     5  181
!     6  305
!
!    The reference has a typographical error for N = 5, OPTION 1, B4:
!
!      -(1)0.736330882774831
!
!    should read
!
!      (-1)0.736330882774831
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  real ( kind = 8 ) b4
  real ( kind = 8 ) b5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  logical more
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_3 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_3 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  if ( n == 3 .and. option == 1 ) then
    u =    0.524647623275290D+00
    v =    0.165068012388578D+01
    b0 = - 0.166705761599566D+02
    b1 =   0.100296981655678D+02
    b2 =   0.161699246687754D+00
    b3 = - 0.604719151221535D+01
    b4 =   0.234381399489666D-01
    b5 =   0.417194501880647D+01
  else if ( n == 3 .and. option == 2 ) then
    u =    0.165068012388578D+01
    v =    0.524647623275290D+00
    b0 =   0.166705761599566D+02
    b1 =   0.178903161957074D+00
    b2 = - 0.665808190965810D+01
    b3 =   0.148361823143070D-01
    b4 =   0.229669852539758D+01
    b5 =   0.430097881732984D-02
  else if ( n == 4 .and. option == 1 ) then
    u  =   0.524647623275290D+00
    v  =   0.165068012388578D+01
    b0 = - 0.167539329651562D+03
    b1 =   0.687922329603575D+02
    b2 =   0.203518409659014D+00
    b3 = - 0.255075279116885D+02
    b4 =   0.415430214106084D-01
    b5 =   0.739458001434961D+01
  else if ( n == 4 .and. option == 2 ) then
    u =    0.165068012388578D+01
    v =    0.524647623275290D+00
    b0 =   0.688432856406677D+02
    b1 =   0.294997847268286D+00
    b2 = - 0.199427272118378D+02
    b3 =   0.110498755408511D-01
    b4 =   0.407079214570997D+01
    b5 =   0.762328646743931D-02
  else if ( n == 5 .and. option == 1 ) then
    u  =   0.524647623275290D+00
    v  =   0.165068012388578D+01
    b0 = - 0.826940846964452D+03
    b1 =   0.264779097660331D+03
    b2 =   0.213460812375320D+00
    b3 = - 0.714240197186780D+02
    b4 =   0.736330882774831D-01
    b5 =   0.131065518222629D+02
  else if ( n == 5 .and. option == 2 ) then
    u =    0.165068012388578D+01
    v =    0.524647623275290D+00
    b0 =   0.220502344940121D+03
    b1 =   0.537746975313769D+00
    b2 = - 0.497781460739792D+02
    b3 = - 0.743845245712926D-02
    b4 =   0.721529121489956D+01
    b5 =   0.135119234557687D-01
  else if ( n == 6 .and. option == 1 ) then
    u  =   0.524647623275290D+00
    v  =   0.165068012388578D+01
    b0 = - 0.309679578630802E+04
    b1 =   0.815423321880237D+03
    b2 =   0.117326937169073D+00
    b3 = - 0.173057295296448D+03
    b4 =   0.130511250871491D+00
    b5 =   0.232307582494626D+02
  else if ( n == 6 .and. option == 2 ) then
    u =    0.165068012388578D+01
    v =    0.524647623275290D+00
    b0 =   0.616293651884027D+03
    b1 =   0.107529736766179D+01
    b2 = - 0.113807008098269D+03
    b3 = - 0.610828352270520D-01
    b4 =   0.127887706992535D+02
    b5 =   0.239492607623178D-01
  end if

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = b0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - u
    w(k) = b1
    k = k + 1
    x(i,k) = + u
    w(k) = b1
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - v
    w(k) = b2
    k = k + 1
    x(i,k) = + v
    w(k) = b2
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = - u
      x(j,k) = + u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = + u
      w(k) = b3
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = - v
      x(j,k) = + v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = + v
      w(k) = b4
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b5
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b5
      end do
    end do
  end do

  return
end
subroutine en_r2_07_3_size ( n, option, o )

!*****************************************************************************80
!
!! EN_R2_07_3_SIZE sizes the Stroud rule 7.3 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   45
!     4   97
!     5  181
!     6  305
!
!    The reference has a typographical error for N = 5, OPTION 1, B4:
!
!      -(1)0.736330882774831
!
!    should read
!
!      (-1)0.736330882774831
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_3_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_3_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  o = ( 4 * n ** 3 + 8 * n + 3 ) / 3

  return
end
subroutine en_r2_09_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_R2_09_1 implements the Stroud rule 9.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
!
!    The rule has precision P = 9.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the 
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   77
!     4  193
!     5  421
!     6  825
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  real ( kind = 8 ) b4
  real ( kind = 8 ) b5
  real ( kind = 8 ) b6
  real ( kind = 8 ) b7
  real ( kind = 8 ) b8
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_09_1 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_09_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  if ( n == 3 ) then
    u =    0.202018287045609D+01
    v =    0.958572464613819D+00
    b0 =   0.676448734429924D+00
    b1 =   0.511989106291551D-02
    b2 =   0.448595723493744D+00
    b3 =   0.235223454595606D-03
    b4 =   0.915390713080005D-01
    b5 =   0.139208199920793D-01
    b6 =   0.235223454595606D-03
    b7 =   0.915390713080008D-01
    b8 =   0.000000000000000D+00
  else if ( n == 4 .and. option == 1 ) then
    u =    0.202018287045609D+01
    v =    0.958572464613819D+00
    b0 = - 0.860452945007048D+00
    b1 = - 0.405511998533795D-01
    b2 =   0.107026475449715D+01
    b3 =   0.138974239307092D-03
    b4 = - 0.162248779448181D+00
    b5 =   0.246740110027234D-01
    b6 =   0.138974239307094D-03
    b7 =   0.162248779448181D+00
    b8 =   0.138974239307094D-03
  else if ( n == 4 .and. option == 2 ) then
    u =    0.958572464613819D+00
    v =    0.202018287045609D+01
    b0 =   0.265029088766810D-02
    b1 =   0.637601342635332D+00
    b2 = - 0.394394059389228D-01
    b3 =   0.540829264827264D-01
    b4 = - 0.416922717921281D-03
    b5 =   0.246740110027234D-01
    b6 =   0.540829264827270D-01
    b7 =   0.416922717921281D-03
    b8 =   0.540829264827269D-01
  else if ( n == 5 .and. option == 1 ) then
    u =    0.202018287045609D+01
    v =    0.958572464613819D+00
    b0 = - 0.827347006200826D+01
    b1 = - 0.160820174530905D+00
    b2 =   0.353499863758467D+01
    b3 =   0.738976276909564D-03
    b4 = - 0.862735421812943D+00
    b5 =   0.437335458190621D-01
    b6 = - 0.246325425636523D-03
    b7 =   0.287578473937648D+00
    b8 =   0.246325425636523D-03
  else if ( n == 5 .and. option == 2 ) then
    u =    0.958572464613819D+00
    v =    0.202018287045609D+01
    b0 = - 0.624416791055272D+00
    b1 =   0.467494915583104D+00
    b2 = - 0.152937760910536D+00
    b3 =   0.287578473937646D+00
    b4 = - 0.221692883072871D-02
    b5 =   0.437335458190621D-01
    b6 = - 0.958594913125490D-01
    b7 =   0.738976276909568D-03
    b8 =   0.958594913125492D-01
  else if ( n == 6 .and. option == 1 ) then
    u =    0.202018287045609D+01
    v =    0.958572464613819D+00
    b0 = - 0.361840434143098D+02
    b1 = - 0.447936529138517D+00
    b2 =   0.112077863004144D+02
    b3 =   0.392940404320855D-02
    b4 = - 0.254859786784158D+01
    b5 =   0.775156917007496D-01
    b6 = - 0.130980134773619D-02
    b7 =   0.509719573568315D+00
    b8 =   0.436600449245395D-03
  else if ( n == 6 .and. option == 2 ) then
    u =    0.958572464613819D+00
    v =    0.202018287045609D+01
    b0 =   0.448873836333650D+01
    b1 = - 0.238473566140736D+01
    b2 = - 0.413008493198885D+00
    b3 =   0.152915872070494D+01
    b4 = - 0.654900673868093D-02
    b5 =   0.775156917007496D-01
    b6 = - 0.509719573568314D+00
    b7 =   0.130980134773618D-02
    b8 =   0.169906524522772D+00
  end if

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = b0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - u
    w(k) = b1
    k = k + 1
    x(i,k) = + u
    w(k) = b1
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - v
    w(k) = b2
    k = k + 1
    x(i,k) = + v
    w(k) = b2
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = - u
      x(j,k) = + u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = + u
      w(k) = b3
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = - v
      x(j,k) = + v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = + v
      w(k) = b4
    end do
  end do
!
!  4 * ( N * ( N - 1 ) ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = - u
      x(j,k) = + v
      w(k) = b5
      k = k + 1
      x(i,k) = + u
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = + u
      x(j,k) = + v
      w(k) = b5
      k = k + 1
      x(i,k) = - v
      x(j,k) = - u
      w(k) = b5
      k = k + 1
      x(i,k) = - v
      x(j,k) = + u
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = - u
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = + u
      w(k) = b5
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b6
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b6
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b7
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b7
      end do
    end do
  end do
!
!  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
!
  do i = 1, n - 3
    do j = i + 1, n - 2
      do l = j + 1, n - 1
        do m = l + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
        end do
      end do
    end do
  end do

  return
end
subroutine en_r2_09_1_size ( n, option, o )

!*****************************************************************************80
!
!! EN_R2_09_1_SIZE sizes the Stroud rule 9.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
!
!    The rule has precision P = 9.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the 
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   77
!     4  193
!     5  421
!     6  825
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_09_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_09_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  o = ( 2 * n ** 4 - 4 * n ** 3 + 22 * n ** 2 - 8 * n + 3 ) / 3

  return
end
subroutine en_r2_11_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_R2_11_1 implements the Stroud rule 11.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order 
!
!      O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
!
!    The rule has precision P = 11.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 5.
!
!     N    O
!    __  ___
!     3  151
!     4  417
!     5  983
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 5.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  real ( kind = 8 ) b4
  real ( kind = 8 ) b5
  real ( kind = 8 ) b6
  real ( kind = 8 ) b7
  real ( kind = 8 ) b8
  real ( kind = 8 ) b9
  real ( kind = 8 ) b10
  real ( kind = 8 ) b11
  real ( kind = 8 ) b12
  real ( kind = 8 ) b13
  real ( kind = 8 ) b14
  real ( kind = 8 ) b15
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i5
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) option
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) w2
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( n < 3 .or. 5 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_11_1 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 5 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_11_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  if ( n == 3 .and. option == 1 ) then
    u =     0.235060497367449D+01
    v =     0.436077411927617D+00
    w2 =    0.133584907401370D+01
    b0 =  - 0.881591029957858D+01
    b1 =  - 0.751996143360650D-01
    b2 =    0.621743189471515D+01
    b3 =    0.241426451456494D+00
    b4 =  - 0.120709739276065D-02
    b5 =  - 0.427751221210138D+01
    b6 =    0.550169924840163D-01
    b7 =    0.237084999634707D-01
    b8 =  - 0.169791992887741D-02
    b9 =  - 0.252266276123350D-04
    b10 =   0.326777873717691D+01
    b11 =   0.968469949206802D-02
    b12 =   0.789754514877422D-03
    b13 =   0.000000000000000D+00
    b14 =   0.000000000000000D+00
    b15 =   0.000000000000000D+00
  else if ( n == 3 .and. option == 2 ) then
    u =     0.235060497367449D+01
    v =     0.133584907401370D+01
    w2 =    0.436077411927617D+00
    b0 =  - 0.141214037032900D+02
    b1 =  - 0.803730274707282D-01
    b2 =    0.235546545595906D+00
    b3 =    0.888123191556611D+01
    b4 =    0.142467131155533D-03
    b5 =    0.582993124006494D-01
    b6 =  - 0.561099173155661D+01
    b7 =  - 0.204028691521686D-02
    b8 =    0.252880089932256D-01
    b9 =  - 0.814378678627283D-04
    b10 =   0.804353953375146D-02
    b11 =   0.393451849690453D+01
    b12 =   0.171183493169724D-03
    b13 =   0.000000000000000D+00
    b14 =   0.000000000000000D+00
    b15 =   0.000000000000000D+00
  else if ( n == 4 .and. option == 1 ) then
    u =     0.235060497367449D+01
    v =     0.436077411927617D+00
    w2 =    0.133584907401370D+01
    b0 =    0.241502736147339D+03
    b1 =  - 0.196095938531478D+00
    b2 =  - 0.128675737999280D+03
    b3 =    0.307568784278696D+00
    b4 =  - 0.480908422319460D-02
    b5 =    0.698087019367085D+02
    b6 =    0.631837143743771D-01
    b7 =    0.392226151971179D-01
    b8 =  - 0.300948471646799D-02
    b9 =  - 0.650235306755170D-04
    b10 = - 0.386951974646715D+02
    b11 =   0.171656829095787D-01
    b12 =   0.139980343116450D-02
    b13 =   0.101552487093372D-04
    b14 =   0.222435922356439D+02
    b15 =   0.000000000000000D+00
  else if ( n == 4 .and. option == 2 ) then
    u =     0.235060497367449D+01
    v =     0.133584907401370D+01
    w2 =    0.436077411927617D+00
    b0 =  - 0.151944464736584D+03
    b1 =  - 0.223498438689039D+00
    b2 =    0.243574919068010D+00
    b3 =    0.634373877008693D+02
    b4 =  - 0.782065187814018D-04
    b5 =    0.911833754536616D-01
    b6 =  - 0.238927288245914D+02
    b7 =  - 0.422314408318853D-02
    b8 =    0.448218289217760D-01
    b9 =  - 0.138053374667391D-03
    b10 =   0.607473265800655D-02
    b11 =   0.697375246129742D+01
    b12 =   0.303414841680135D-03
    b13 = - 0.314574391771792D-05
    b14 =   0.409103498175100D-02
    b15 =   0.000000000000000D+00
  else if ( n == 5 .and. option == 1 ) then
    u =     0.235060497367449D+01
    v =     0.436077411927617D+00
    w2 =    0.133584907401370D+01
    b0 =    0.255885269311763E+04
    b1 =  - 0.439598677491526D+00
    b2 =  - 0.106541406144610E+04
    b3 =    0.453540909054264D+00
    b4 =  - 0.132100905623778D-01
    b5 =    0.418606568954203D+03
    b6 =    0.511394563043680D-01
    b7 =    0.645581013845604D-01
    b8 =  - 0.533417277494500D-02
    b9 =  - 0.137981626254496D-03
    b10 = - 0.147436933189884D+03
    b11 =   0.304253807765057D-01
    b12 =   0.248108698207828D-02
    b13 =   0.113652094546015D-04
    b14 =   0.394257407160391D+02
    b15 =   0.331725011358320D-05
  else if ( n == 5 .and. option == 2 ) then
    u =     0.235060497367449D+01
    v =     0.133584907401370D+01
    w2 =    0.436077411927617D+00
    b0 =  - 0.761305347548192D+03
    b1 =  - 0.536360805019297D+00
    b2 =    0.110669832078736D+00
    b3 =    0.246421088923968D+03
    b4 =  - 0.773649327968607D-03
    b5 =    0.169088641205970D+00
    b6 =  - 0.670700680243651D+02
    b7 =  - 0.856090560229205D-02
    b8 =    0.794446232770302D-01
    b9 =  - 0.220272863263544D-03
    b10 = - 0.373515812228225D-02
    b11 =   0.123606544052884D+02
    b12 =   0.537788804557843D-03
    b13 = - 0.122101861480881D-04
    b14 =   0.725117070759373D-02
    b15 =   0.331725011358320D-05
  end if

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = 0.0D+00
  w(k) = b0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - u
    w(k) = b1
    k = k + 1
    x(i,k) = + u
    w(k) = b1
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - v
    w(k) = b2
    k = k + 1
    x(i,k) = + v
    w(k) = b2
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - w2
    w(k) = b3
    k = k + 1
    x(i,k) = + w2
    w(k) = b3
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - u
      w(k) = b4
      k = k + 1
      x(i,k) = - u
      x(j,k) = + u
      w(k) = b4
      k = k + 1
      x(i,k) = + u
      x(j,k) = - u
      w(k) = b4
      k = k + 1
      x(i,k) = + u
      x(j,k) = + u
      w(k) = b4
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - v
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = - v
      x(j,k) = + v
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = + v
      w(k) = b5
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - w2
      x(j,k) = - w2
      w(k) = b6
      k = k + 1
      x(i,k) = - w2
      x(j,k) = + w2
      w(k) = b6
      k = k + 1
      x(i,k) = + w2
      x(j,k) = - w2
      w(k) = b6
      k = k + 1
      x(i,k) = + w2
      x(j,k) = + w2
      w(k) = b6
    end do
  end do
!
!  4 * ( N * ( N - 1 ) ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - v
      w(k) = b7
      k = k + 1
      x(i,k) = - u
      x(j,k) = + v
      w(k) = b7
      k = k + 1
      x(i,k) = + u
      x(j,k) = - v
      w(k) = b7
      k = k + 1
      x(i,k) = + u
      x(j,k) = + v
      w(k) = b7
      k = k + 1
      x(i,k) = - v
      x(j,k) = - u
      w(k) = b7
      k = k + 1
      x(i,k) = - v
      x(j,k) = + u
      w(k) = b7
      k = k + 1
      x(i,k) = + v
      x(j,k) = - u
      w(k) = b7
      k = k + 1
      x(i,k) = + v
      x(j,k) = + u
      w(k) = b7
    end do
  end do
!
!  4 * ( N * ( N - 1 ) ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - w2
      w(k) = b8
      k = k + 1
      x(i,k) = - u
      x(j,k) = + w2
      w(k) = b8
      k = k + 1
      x(i,k) = + u
      x(j,k) = - w2
      w(k) = b8
      k = k + 1
      x(i,k) = + u
      x(j,k) = + w2
      w(k) = b8
      k = k + 1
      x(i,k) = - w2
      x(j,k) = - u
      w(k) = b8
      k = k + 1
      x(i,k) = - w2
      x(j,k) = + u
      w(k) = b8
      k = k + 1
      x(i,k) = + w2
      x(j,k) = - u
      w(k) = b8
      k = k + 1
      x(i,k) = + w2
      x(j,k) = + u
      w(k) = b8
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b9
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b9
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b10
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b10
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - w2
        x(j,k) = - w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = - w2
        x(j,k) = - w2
        x(l,k) = + w2
        w(k) = b11
        k = k + 1
        x(i,k) = - w2
        x(j,k) = + w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = - w2
        x(j,k) = + w2
        x(l,k) = + w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = - w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = - w2
        x(l,k) = + w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = + w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = + w2
        x(l,k) = + w2
        w(k) = b11
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 2 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = - v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = - v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b12
      end do
    end do
  end do
!
!  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
!
  do i = 1, n - 3
    do j = i + 1, n - 2
      do l = j + 1, n - 1
        do m = l + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
        end do
      end do
    end do
  end do
!
!  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
!
  do i = 1, n - 3
    do j = i + 1, n - 2
      do l = j + 1, n - 1
        do m = l + 1, n
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
        end do
      end do
    end do
  end do
!
!  All quintuples UUUUU with 32 sign combinations.
!
  do i1 = 1, n - 4
    do i2 = i1 + 1, n - 3
      do i3 = i2 + 1, n - 2
        do i4 = i3 + 1, n - 1
          do i5 = i4 + 1, n
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
          end do
        end do
      end do
    end do
  end do

  return
end
subroutine en_r2_11_1_size ( n, option, o )

!*****************************************************************************80
!
!! EN_R2_11_1_SIZE sizes the Stroud rule 11.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order 
!
!      O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
!
!    The rule has precision P = 11.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 5.
!
!     N    O
!    __  ___
!     3  151
!     4  417
!     5  983
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 5.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  if ( n < 3 .or. 5 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_11_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 5 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_11_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  o = ( 4 * n ** 5 - 20 * n ** 4 + 140 * n ** 3 - 130 * n ** 2 &
    + 96 * n + 15 ) / 15


  return
end
subroutine en_r2_monomial_integral ( n, alpha, value )

!*****************************************************************************80
!
!! EN_R2_MONOMIAL_INTEGRAL evaluates monomial integrals in EN_R2.
!
!  Discussion:
!
!    ALPHA is the set of polynomial exponents.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
!
!    The integral to be evaluated is
!
!      value = integral ( EN ) x(1)^alpha(1) * x(2)^alpha(2) * ... 
!        * x(n)^alpha(n) * w(x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ALPHA(N), the polynomial exponents.
!    0 <= ALPHA(*).
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) alpha(n)
  real ( kind = 8 ) arg
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) value

  if ( any ( alpha(1:n) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  Some ALPHA(I) < 0.'
    stop
  else if ( any ( mod ( alpha(1:n), 2 ) == 1 ) ) then
    value = 0.0D+00
  else
    value = 1.0D+00
    do i = 1, n
      arg = ( real ( alpha(i) + 1, kind = 8 ) ) / 2.0D+00
      value = value * r8_gamma ( arg )
    end do
  end if

  return
end
subroutine ep1_glg_monomial_integral ( expon, alpha, exact )

!*****************************************************************************80
!
!! EP1_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EP1.
!
!  Discussion:
!
!    EP1_GLG is the interval [0,+oo) with generalized Laguerre weight function:
!
!      w(alpha;x) = x^alpha exp ( - x )
!
!    value = integral ( 0 <= x < +oo ) x^expon x^alpha exp ( - x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) EXACT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_gamma

  arg = alpha + real ( expon + 1, kind = 8 )

  exact = r8_gamma ( arg )

  return
end
subroutine ep1_lag_monomial_integral ( expon, value )

!*****************************************************************************80
!
!! EP1_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EP1.
!
!  Discussion:
!
!    EP1 is the interval [0,+oo) with exponential or Laguerre weight function:
!
!      w(x) = exp ( - x )
!
!    value = integral ( 0 <= x < oo ) x^expon exp ( - x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) expon
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) value

  value = r8_factorial ( expon )

  return
end
subroutine epn_glg_00_1 ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! EPN_GLG_00_1 implements the "midpoint rule" for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  expon = 0
  call ep1_glg_monomial_integral ( expon, alpha, volume )
  volume = volume ** n

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = 1.0D+00
  w(k) = volume

  return
end
subroutine epn_glg_00_1_size ( n, alpha, o )

!*****************************************************************************80
!
!! EPN_GLG_00_1_SIZE sizes the midpoint rule for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine epn_glg_01_1 ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! EPN_GLG_01_1 implements a precision 1 rule for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  expon = 0
  call ep1_glg_monomial_integral ( expon, alpha, value1 )
  volume = value1 ** n

  expon = 1
  call ep1_glg_monomial_integral ( expon, alpha, value2 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = value2 / value1
  w(k) = volume

  return
end
subroutine epn_glg_01_1_size ( n, alpha, o )

!*****************************************************************************80
!
!! EPN_GLG_01_1_SIZE sizes a precision 1 rule for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = 1

  return
end
subroutine epn_glg_02_xiu ( n, alpha, o, x, w )

!*****************************************************************************80
!
!! EPN_GLG_02_XIU implements the Xiu precision 2 rule for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) coef
  real ( kind = 8 ) delta0
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_02_XIU - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = - 1.0D+00
  delta0 = alpha + 1.0D+00
  c1 = - alpha - 1.0D+00

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  expon = 0
  call ep1_glg_monomial_integral ( expon, alpha, volume_1d )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine epn_glg_02_xiu_size ( n, alpha, o )

!*****************************************************************************80
!
!! EPN_GLG_02_XIU_SIZE sizes the Xiu rule for region EPN_GLG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPN_GLG_00_1 - Fatal error!'
    write ( *, '(a)' ) '  ALPHA <= -1.0'
    stop
  end if

  o = n + 1

  return
end
subroutine epn_glg_monomial_integral ( n, expon, alpha, value )

!*****************************************************************************80
!
!! EPN_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EPN.
!
!  Discussion:
!
!    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
!    Laguerre weight function:
!
!      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
!
!    value = integral ( EPN ) 
!      product ( 1 <= i <= n ) x(I)^expon(i) x(i)^alpha exp ( - x(i) ) dx(i)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) EXPON(N), the exponents.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of X in the weight function.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) expon(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) value
  real ( kind = 8 ) value2

  value = 1.0D+00
  do i = 1, n
    call ep1_glg_monomial_integral ( expon(i), alpha, value2 )
    value = value * value2
  end do

  return
end
subroutine epn_lag_00_1 ( n, o, x, w )

!*****************************************************************************80
!
!! EPN_LAG_00_1 implements the "midpoint rule" for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  expon = 0
  call ep1_lag_monomial_integral ( expon, volume )
  volume = volume ** n

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = 1.0D+00
  w(k) = volume

  return
end
subroutine epn_lag_00_1_size ( n, o )

!*****************************************************************************80
!
!! EPN_LAG_00_1_SIZE sizes the midpoint rule for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 0.
!
!    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 1

  return
end
subroutine epn_lag_01_1 ( n, o, x, w )

!*****************************************************************************80
!
!! EPN_LAG_01_1 implements a precision 1 rule for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  integer ( kind = 4 ) expon
  integer ( kind = 4 ) k
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  expon = 0
  call ep1_lag_monomial_integral ( expon, value1 )
  volume = value1 ** n

  expon = 1
  call ep1_lag_monomial_integral ( expon, value2 )

  x(1:n,1:o) = 0.0D+00

  k = 0
!
!  1 point.
!
  k = k + 1
  x(1:n,k) = value2 / value1
  w(k) = volume

  return
end
subroutine epn_lag_01_1_size ( n, o )

!*****************************************************************************80
!
!! EPN_LAG_01_1_SIZE sizes a precision 1 rule for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 1

  return
end
subroutine epn_lag_02_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! EPN_LAG_02_XIU implements the Xiu precision 2 rule for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EPN_LAG is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) coef
  real ( kind = 8 ) delta0
  integer ( kind = 4 ) expon
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  do j = 1, o

    i = 0 
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )

      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg ) 
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do

  gamma0 = - 1.0D+00
  delta0 = 1.0D+00
  c1 = - 1.0D+00

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  expon = 0
  call ep1_lag_monomial_integral ( expon, volume_1d )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind = 8 )

  return
end
subroutine epn_lag_02_xiu_size ( n, o )

!*****************************************************************************80
!
!! EPN_LAG_02_XIU_SIZE sizes the Xiu rule for region EPN_LAG.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n + 1

  return
end
subroutine epn_lag_monomial_integral ( n, expon, value )

!*****************************************************************************80
!
!! EPN_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EPN.
!
!  Discussion:
!
!    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
!    or Laguerre weight function:
!
!      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
!
!    value = integral ( EPN ) 
!      product ( 1 <= i <= n ) x(I)^expon(i) exp ( -x(i) ) dx(i)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) EXPON(N), the exponents.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) expon(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) value
  real ( kind = 8 ) value2

  value = 1.0D+00
  do i = 1, n
    call ep1_lag_monomial_integral ( expon(i), value2 )
    value = value * value2
  end do

  return
end
subroutine gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )

!*****************************************************************************80
!
!! GW_02_XIU implements the Golub-Welsch version of the Xiu rule.
!
!  Discussion:
!
!    The rule has order 
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    It is assumed that the integral is over an N-dimensional region,
!    and has the form
!
!      Integral f(x) w(x) dx
!
!    where w(x) is separable into identical and independent components:
!
!      w(x) = v(x1) * v(x2) * ... * v(xn)
!
!    Associated with the weight function v(x), we assume there is a
!    family of orthogonal polynomials satisfying a three-term recurrence
!    of the form:
!
!      x P(n,x) = An * P(n+1,x) + Bn * P(n,x) + Cn * P(n-1,x)
!
!    with P(0,x) = 1, and P(-1,x) = 0.
!
!    This routine can construct the desired quadrature rule by knowing
!    the values of C1, used in the definition of P2, the values
!    GAMMA0 = 1/A0 and DELTA0 = - B0/A0, for which it is the case that
!    P(1,X) = GAMMA0 * X + DELTA0, and the value of VOLUME_1D, that is,
!    the 1D integral of v(x) over the region.
!
!    Note the values for the following standard polynomial families:
!
!    Chebyshev Type 1
!      V(X) =      1 / sqrt ( 1 - X^2 )
!      Interval =  [-1,+1]
!      GAMMA0 =    1.0
!      DELTA0 =    0.0
!      C1 =        1/2
!      VOLUME_1D = PI
!
!    Chebyshev Type 2
!      V(X) =      sqrt ( 1 - X^2 )
!      Interval =  [-1,+1]
!      GAMMA0 =    2.0
!      DELTA0 =    0.0
!      C1 =        1/2
!      VOLUME_1D = PI / 2
!
!    Gegenbauer
!      V(X) =      ( 1 - X^2 )^A
!      Interval =  [-1,+1]
!      GAMMA0 =    2 * A + 1
!      DELTA0 =    0.0
!      C1 =        ( 2 * A + 1 ) / ( 2 A + 3 )
!      VOLUME_1D = sqrt ( PI ) * Gamma(A+1) / Gamma(A+3/2)
!
!    Gegenbauer* (Removes singularity at ALPHA = -0.5):
!      V(X) =      ( 1 - X^2 )^A
!      Interval =  [-1,+1]
!      GAMMA0 =    1
!      DELTA0 =    0.0
!      C1 =        1 / ( 2 A + 3 )
!      VOLUME_1D = sqrt ( PI ) * Gamma(A+1) / Gamma(A+3/2)
!
!    Generalized Hermite
!      V(X) = |x|^A exp ( - x^2 )
!      Interval = (-oo,+oo)
!      GAMMA0 =    2
!      DELTA0 =    0
!      C1 =        2+2A
!      VOLUME_1D = Gamma((A+1)/2)
!
!    Generalized Laguerre
!      V(X) =       x^A exp ( - x )
!      Interval =  [0,+oo)
!      GAMMA0 =    -1.0
!      DELTA0 =     A+1.0
!      C1 =        -A-1.0
!      VOLUME_1D =  Gamma(A+1)
!
!    Hermite (physicist)
!      V(X) =      exp ( - x^2 )
!      Interval =  (-oo,+oo)
!      GAMMA0 =    2.0
!      DELTA0 =    0.0
!      C1 =        1.0
!      VOLUME_1D = sqrt ( PI )
!
!    Hermite (probabilist)
!      V(X) =      exp ( - x^2 / 2 )
!      Interval =  (-oo,+oo)
!      GAMMA0 =    1.0
!      DELTA0 =    0.0
!      C1 =        1.0
!      VOLUME_1D = sqrt ( 2 PI )
!
!    Jacobi
!      V(X) =      (1-x)^A (1+x)^B
!      Interval =  [-1,+1]
!      GAMMA0 =    (A+B+2)/2  
!      DELTA0 =    (A-B)/2
!      C1 =        2(A+1)(B+1)/(A+B+3)/(A+B+2)
!      VOLUME_1D = 2^(A+B+1) * Gamma(A+1) * Gamma(B+1) / ( A+B+1) / Gamma(A+B+1)
!
!    Laguerre
!      V(X) =       exp ( - x )
!      Interval =  [0,+oo)
!      GAMMA0 =    -1.0
!      DELTA0 =     1.0
!      C1 =        -1.0
!      VOLUME_1D =  1.0
!
!    Legendre
!      V(X) =      1.0
!      Interval =  [-1,+1]
!      GAMMA0 =    1.0
!      DELTA0 =    0.0
!      C1 =        1/3
!      VOLUME_1D = 2.0
!                                  
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Input, real ( kind = 8 ) GAMMA0, the ratio 1 / A0.
!
!    Input, real ( kind = 8 ) DELTA0, the ratio B0 / A0.
!
!    Input, real ( kind = 8 ) C1, the coefficient of P(0,X) in 
!    the definition of P(2,X).
!
!    Input, real ( kind = 8 ) VOLUME_1D, the 1D integral of V(X).
!
!    Output, real ( kind = 8 ) X(N,O), the abscissas.
!
!    Output, real ( kind = 8 ) W(O), the weights.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind = 8 ) arg
  real ( kind = 8 ) c1
  real ( kind = 8 ) delta0
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) r
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ) w(o)
  real ( kind = 8 ) x(n,o)

  do j = 1, o

    i = 0
    do r = 1, ( n / 2 )
      arg = real ( 2 * r * ( j - 1 ), kind = 8 ) * pi / real ( n + 1, kind = 8 )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * cos ( arg )
      i = i + 1
      x(i,j) = sqrt ( 2.0D+00 ) * sin ( arg )
    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = r8_mop ( j - 1 )
    end if

  end do
!
!  Adjust for the GW rule.
!
  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0
!
!  The weights are equal.
!
  w(1:o) = volume_1d ** n  / real ( o, kind = 8 )

  return
end
subroutine gw_02_xiu_size ( n, o )

!*****************************************************************************80
!
!! GW_02_XIU_SIZE sizes the Golub Welsch version of the Xiu rule.
!
!  Discussion:
!
!    The rule has order O = N + 1;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n + 1

  return
end
function hexagon_area_2d ( r )

!*****************************************************************************80
!
!! HEXAGON_AREA_2D returns the area of a regular hexagon in 2D.
!
!  Discussion:
!
!    The formula for the area only requires the radius, and does
!    not depend on the location of the center, or the orientation
!    of the hexagon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the hexagon.
!
!    Output, real ( kind = 8 ) HEXAGON_AREA_2D, the area of the hexagon.
!
  implicit none

  real ( kind = 8 ) hexagon_area_2d
  real ( kind = 8 ) hexagon_unit_area_2d
  real ( kind = 8 ) r

  hexagon_area_2d = r * r * hexagon_unit_area_2d ( )

  return
end
subroutine hexagon_sum ( func, center, r, order, xtab, ytab, weight, &
  result )

!*****************************************************************************80
!
!! HEXAGON_SUM applies a quadrature rule inside a hexagon in 2D.
!
!  Discussion:
!
!    The input quadrature rule is assumed to be defined for a unit hexagon.
!
!    The input quadrature rule may be defined by calling HEXAGON_UNIT_SET.
!
!  Integration region:
!
!    The definition is given in terms of THETA, the angle in degrees of the
!    vector (X-CENTER(1),Y-CENTER(2)).  The following six conditions apply,
!    respectively, between the bracketing values of THETA of 0, 60, 120, 
!    180, 240, 300, and 360.
!
!      0 <= Y-CENTER(2) <= -SQRT(3) * (X-CENTER(1)) + R * SQRT(3)
!      0 <= Y-CENTER(2) <=                     R * SQRT(3)/2
!      0 <= Y-CENTER(2) <=  SQRT(3) * (X-CENTER(1)) + R * SQRT(3) 
!      -SQRT(3) * (X-CENTER(1)) - R * SQRT(3)	<= Y-CENTER(2) <= 0
!                        - R * SQRT(3)/2 <= Y-CENTER(2) <= 0
!       SQRT(3) * (X-CENTER(1)) - R * SQRT(3)   <= Y-CENTER(2) <= 0
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
!    Input, external FUNC, the name of the user supplied
!    function of two variables which is to be integrated,
!    of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, real ( kind = 8 ) CENTER(2), the center of the hexagon.
!
!    Input, real ( kind = 8 ) R, the radius of the hexagon.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) order

  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) hexagon_area_2d
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) r
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) y
  real ( kind = 8 ) ytab(order)

  quad = 0.0D+00
  do i = 1, order
    x = center(1) + r * xtab(i)
    y = center(2) + r * ytab(i)
    quad = quad + weight(i) * func ( x, y )
  end do

  volume = hexagon_area_2d ( r )
  result = quad * volume

  return
end
function hexagon_unit_area_2d ( )

!*****************************************************************************80
!
!! HEXAGON_UNIT_AREA_2D returns the area of the unit regular hexagon in 2D.
!
!  Integration region:
!
!    The definition is given in terms of THETA, the angle in degrees of the
!    vector (X,Y).  The following six conditions apply, respectively,
!    between the bracketing values of THETA of 0, 60, 120, 180, 240,
!    300, and 360.
!
!                              0 <= Y <= -SQRT(3) * X + SQRT(3)
!                              0 <= Y <=                 SQRT(3)/2
!                              0 <= Y <=  SQRT(3) * X + SQRT(3)
!      - SQRT(3) * X - SQRT(3)   <= Y <= 0
!                    - SQRT(3)/2 <= Y <= 0
!        SQRT(3) * X - SQRT(3)   <= Y <= 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) HEXAGON_UNIT_AREA_2D, the area of the hexagon.
!
  implicit none

  real ( kind = 8 ) hexagon_unit_area_2d

  hexagon_unit_area_2d = 3.0D+00 * sqrt ( 3.0D+00 ) / 2.0D+00

  return
end
subroutine hexagon_unit_set ( rule, order, xtab, ytab, weight )

!*****************************************************************************80
!
!! HEXAGON_UNIT_SET sets a quadrature rule inside the unit hexagon in 2D.
!
!  Integration region:
!
!    The definition is given in terms of THETA, the angle in degrees of the
!    vector (X,Y).  The following six conditions apply, respectively,
!    between the bracketing values of THETA of 0, 60, 120, 180, 240,
!    300, and 360.
!
!                              0 <= Y <= -SQRT(3) * X + SQRT(3)
!                              0 <= Y <=                 SQRT(3)/2
!                              0 <= Y <=  SQRT(3) * X + SQRT(3)
!       -SQRT(3) * X - SQRT(3)   <= Y <= 0
!                    - SQRT(3)/2 <= Y <= 0
!        SQRT(3) * X - SQRT(3)   <= Y <= 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the rule desired.
!      1, 1 point,  degree 1;
!      2, 4 points, degree 3;
!      3, 7 points, degree 3;
!      4, 7 points, degree 5;
!
!    Input, integer ( kind = 4 ) ORDER, the order of the desired rule.
!
!    Output, real ( kind = 8 ) XTAB(*), YTAB(*), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(*), the ORDER weights of the rule.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  integer ( kind = 4 ) rule
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) z

  if ( rule == 1 ) then

    xtab(1) = 0.0D+00
    ytab(1) = 0.0D+00
    weight(1) = 1.0D+00
!
!  Stroud rule H2:3-1.
!
  else if ( rule == 2 ) then

    a = sqrt ( 5.0D+00 / 12.0D+00 )
    b = 1.0D+00 / 4.0D+00
    z = 0.0D+00

    xtab(1:order) =   (/  a, -a,  z,  z /)
    ytab(1:order) =   (/  z,  z,  a, -a /)
    weight(1:order) = (/  b,  b,  b,  b /)
!
!  Stroud rule H2:3-2.
!
  else if ( rule == 3 ) then

    a = sqrt ( 3.0D+00 ) / 2.0D+00
    b =  0.5D+00
    c =  1.0D+00
    d =  5.0D+00 / 72.0D+00
    e = 42.0D+00 / 72.0D+00
    z =  0.0D+00

    xtab(1:order) =   (/  z,  c, -c,  b, -b,  b, -b /)
    ytab(1:order) =   (/  z,  z,  z,  a,  a, -a, -a /)
    weight(1:order) = (/  e,  d,  d,  d,  d,  d,  d /)
!
!  Stroud rule H2:5-1.
!
  else if ( rule == 4 ) then

    a = sqrt ( 14.0D+00 ) / 5.0D+00
    b = sqrt ( 14.0D+00 ) / 10.0D+00
    c = sqrt ( 42.0D+00 ) / 10.0D+00
    d = 125.0D+00 / 1008.0D+00
    e = 258.0D+00 / 1008.0D+00
    z = 0.0D+00

    xtab(1:order) =   (/ z,  a, -a,  b, -b,  b, -b /)
    ytab(1:order) =   (/ z,  z,  z,  c,  c, -c, -c /)
    weight(1:order) = (/ e,  d,  d,  d,  d,  d,  d /)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEXAGON_UNIT_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of RULE = ', rule
    stop

  end if

  return
end
subroutine hexagon_unit_size ( rule, order )

!*****************************************************************************80
!
!! HEXAGON_UNIT_SIZE sizes a quadrature rule inside the unit hexagon in 2D.
!
!  Integration region:
!
!    The definition is given in terms of THETA, the angle in degrees of the
!    vector (X,Y).  The following six conditions apply, respectively,
!    between the bracketing values of THETA of 0, 60, 120, 180, 240,
!    300, and 360.
!
!                              0 <= Y <= -SQRT(3) * X + SQRT(3)
!                              0 <= Y <=                 SQRT(3)/2
!                              0 <= Y <=  SQRT(3) * X + SQRT(3)
!       -SQRT(3) * X - SQRT(3)   <= Y <= 0
!                    - SQRT(3)/2 <= Y <= 0
!        SQRT(3) * X - SQRT(3)   <= Y <= 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the rule desired.
!      1, 1 point,  degree 1;
!      2, 4 points, degree 3;
!      3, 7 points, degree 3;
!      4, 7 points, degree 5;
!
!    Output, integer ( kind = 4 ) ORDER, the order of the desired rule.
!    If RULE is not legal, then ORDER is returned as -1.
!
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ) rule

  if ( rule == 1 ) then

    order = 1
!
!  Stroud rule H2:3-1.
!
  else if ( rule == 2 ) then

    order = 4
!
!  Stroud rule H2:3-2.
!
  else if ( rule == 3 ) then

    order = 7
!
!  Stroud rule H2:5-1.
!
  else if ( rule == 4 ) then

    order = 7

  else

    order = -1

  end if

  return
end
function i4_factorial ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL computes N! (for small values of N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument.
!
!    Output, integer ( kind = 4 ) I4_FACTORIAL, the value of N!.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) i4_factorial

  i4_factorial = 1
  do i = 1, n
    i4_factorial = i4_factorial * i
  end do

  return
end
function i4_factorial2 ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the double factorial 
!    function.  If N is less than 1, I4_FACTORIAL2 is returned as 1.
!
!    Output, integer ( kind = 4 ) I4_FACTORIAL2, the value of N!!.
!
  implicit none

  integer ( kind = 4 ) i4_factorial2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_copy

  if ( n < 1 ) then
    i4_factorial2 = 1
    return
  end if

  n_copy = n
  i4_factorial2 = 1

  do while ( 1 < n_copy )
    i4_factorial2 = i4_factorial2 * n_copy
    n_copy = n_copy - 2
  end do

  return
end
subroutine ksub_next2 ( n, k, iarray, in, iout )

!*****************************************************************************80
!
!! KSUB_NEXT2 computes the next K subset of an N set.
!
!  Discussion:
!
!    This routine uses the revolving door method.  It has no "memory".
!    It simply calculates the successor of the input set,
!    and will start from the beginning after the last set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 March 2001
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
!    Second edition,
!    Academic Press, 1978,
!    ISBN 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set from which subsets 
!    are drawn.  N must be positive.
!
!    Input, integer ( kind = 4 ) K, the size of the desired subset.  K must be
!    between 0 and N.
!
!    Input/output, integer ( kind = 4 ) IARRAY(K).  On input, the user must
!    supply a subset of size K in IARRAY.  That is, IARRAY must
!    contain K unique numbers, in order, between 1 and N.  On
!    output, IARRAY(I) is the I-th element of the output subset.
!    The output array is also in sorted order.
!
!    Output, integer ( kind = 4 ) IN, the element of the output subset which
!    was not in the input set.  Each new subset differs from the
!    last one by adding one element and deleting another.
!
!    Output, integer ( kind = 4 ) IOUT, the element of the input subset which
!    is not in the output subset.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) iarray(k)
  integer ( kind = 4 ) in
  integer ( kind = 4 ) iout
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT2 - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a)' ) '  but 0 < N is required!'
    stop
  end if

  if ( k < 0 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT2 - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if

  j = 0

  do

    if ( 0 < j .or. mod ( k, 2 ) == 0 ) then

      j = j + 1

      if ( k < j ) then
        iarray(k) = k
        in = k
        iout = n
        return
      end if

      if ( iarray(j) /= j ) then

        iout = iarray(j)
        in = iout - 1
        iarray(j) = in

        if ( j /= 1 ) then
          in = j - 1
          iarray(j-1) = in
        end if

        return

      end if

    end if

    j = j + 1
    m = n

    if ( j < k ) then
      m = iarray(j+1) - 1
    end if

    if ( m /= iarray(j) ) then
      exit
    end if

  end do

  in = iarray(j) + 1
  iarray(j) = in
  iout = in - 1

  if ( j /= 1 ) then
    iarray(j-1) = iout
    iout = j - 1
  end if

  return
end
subroutine legendre_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = 1.0.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X^(2*N-1).
!
!    The abscissas of the rule are the zeroes of the Legendre polynomial
!    P(N)(X).
!
!    The integral produced by a Gauss-Legendre rule is equal to the
!    integral of the unique polynomial of degree N-1 which
!    agrees with the function at the ORDER abscissas of the rule.
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
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    Dover, 2006,
!    ISBN: 0486445798,
!    LC: QA311.K713.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!    N must be between 1 and 33 or 63, 64, 65, 127 or 255.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) WN), the weights.
!    The weights are positive, symmetric and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) =   0.0D+00

    w(1) = 2.0D+00

  else if ( n == 2 ) then

    x(1) = -0.577350269189625764509148780502D+00
    x(2) =  0.577350269189625764509148780502D+00

    w(1) = 1.0D+00
    w(2) = 1.0D+00

  else if ( n == 3 ) then

    x(1) = -0.774596669241483377035853079956D+00
    x(2) =  0.000000000000000000000000000000D+00
    x(3) =  0.774596669241483377035853079956D+00

    w(1) = 5.0D+00 / 9.0D+00
    w(2) = 8.0D+00 / 9.0D+00
    w(3) = 5.0D+00 / 9.0D+00

  else if ( n == 4 ) then

    x(1) = -0.861136311594052575223946488893D+00
    x(2) = -0.339981043584856264802665759103D+00
    x(3) =  0.339981043584856264802665759103D+00
    x(4) =  0.861136311594052575223946488893D+00

    w(1) = 0.347854845137453857373063949222D+00
    w(2) = 0.652145154862546142626936050778D+00
    w(3) = 0.652145154862546142626936050778D+00
    w(4) = 0.347854845137453857373063949222D+00

  else if ( n == 5 ) then

    x(1) = -0.906179845938663992797626878299D+00
    x(2) = -0.538469310105683091036314420700D+00
    x(3) =  0.000000000000000000000000000000D+00
    x(4) =  0.538469310105683091036314420700D+00
    x(5) =  0.906179845938663992797626878299D+00

    w(1) = 0.236926885056189087514264040720D+00
    w(2) = 0.478628670499366468041291514836D+00
    w(3) = 0.568888888888888888888888888889D+00
    w(4) = 0.478628670499366468041291514836D+00
    w(5) = 0.236926885056189087514264040720D+00

  else if ( n == 6 ) then

    x(1) = - 0.932469514203152027812301554494D+00
    x(2) = - 0.661209386466264513661399595020D+00
    x(3) = - 0.238619186083196908630501721681D+00
    x(4) =   0.238619186083196908630501721681D+00
    x(5) =   0.661209386466264513661399595020D+00
    x(6) =   0.932469514203152027812301554494D+00

    w(1) = 0.171324492379170345040296142173D+00
    w(2) = 0.360761573048138607569833513838D+00
    w(3) = 0.467913934572691047389870343990D+00
    w(4) = 0.467913934572691047389870343990D+00
    w(5) = 0.360761573048138607569833513838D+00
    w(6) = 0.171324492379170345040296142173D+00

  else if ( n == 7 ) then

    x(1) = - 0.949107912342758524526189684048D+00
    x(2) = - 0.741531185599394439863864773281D+00
    x(3) = - 0.405845151377397166906606412077D+00
    x(4) =   0.0D+00
    x(5) =   0.405845151377397166906606412077D+00
    x(6) =   0.741531185599394439863864773281D+00
    x(7) =   0.949107912342758524526189684048D+00

    w(1) = 0.129484966168869693270611432679D+00
    w(2) = 0.279705391489276667901467771424D+00
    w(3) = 0.381830050505118944950369775489D+00
    w(4) = 0.417959183673469387755102040816D+00
    w(5) = 0.381830050505118944950369775489D+00
    w(6) = 0.279705391489276667901467771424D+00
    w(7) = 0.129484966168869693270611432679D+00

  else if ( n == 8 ) then

    x(1) = - 0.960289856497536231683560868569D+00
    x(2) = - 0.796666477413626739591553936476D+00
    x(3) = - 0.525532409916328985817739049189D+00
    x(4) = - 0.183434642495649804939476142360D+00
    x(5) =   0.183434642495649804939476142360D+00
    x(6) =   0.525532409916328985817739049189D+00
    x(7) =   0.796666477413626739591553936476D+00
    x(8) =   0.960289856497536231683560868569D+00

    w(1) = 0.101228536290376259152531354310D+00
    w(2) = 0.222381034453374470544355994426D+00
    w(3) = 0.313706645877887287337962201987D+00
    w(4) = 0.362683783378361982965150449277D+00
    w(5) = 0.362683783378361982965150449277D+00
    w(6) = 0.313706645877887287337962201987D+00
    w(7) = 0.222381034453374470544355994426D+00
    w(8) = 0.101228536290376259152531354310D+00

  else if ( n == 9 ) then

    x(1) = - 0.968160239507626089835576202904D+00
    x(2) = - 0.836031107326635794299429788070D+00
    x(3) = - 0.613371432700590397308702039341D+00
    x(4) = - 0.324253423403808929038538014643D+00
    x(5) =   0.0D+00
    x(6) =   0.324253423403808929038538014643D+00
    x(7) =   0.613371432700590397308702039341D+00
    x(8) =   0.836031107326635794299429788070D+00
    x(9) =   0.968160239507626089835576202904D+00

    w(1) = 0.812743883615744119718921581105D-01
    w(2) = 0.180648160694857404058472031243D+00
    w(3) = 0.260610696402935462318742869419D+00
    w(4) = 0.312347077040002840068630406584D+00
    w(5) = 0.330239355001259763164525069287D+00
    w(6) = 0.312347077040002840068630406584D+00
    w(7) = 0.260610696402935462318742869419D+00
    w(8) = 0.180648160694857404058472031243D+00
    w(9) = 0.812743883615744119718921581105D-01

  else if ( n == 10 ) then

    x(1) =  - 0.973906528517171720077964012084D+00
    x(2) =  - 0.865063366688984510732096688423D+00
    x(3) =  - 0.679409568299024406234327365115D+00
    x(4) =  - 0.433395394129247190799265943166D+00
    x(5) =  - 0.148874338981631210884826001130D+00
    x(6) =    0.148874338981631210884826001130D+00
    x(7) =    0.433395394129247190799265943166D+00
    x(8) =    0.679409568299024406234327365115D+00
    x(9) =    0.865063366688984510732096688423D+00
    x(10) =   0.973906528517171720077964012084D+00

    w(1) =  0.666713443086881375935688098933D-01
    w(2) =  0.149451349150580593145776339658D+00
    w(3) =  0.219086362515982043995534934228D+00
    w(4) =  0.269266719309996355091226921569D+00
    w(5) =  0.295524224714752870173892994651D+00
    w(6) =  0.295524224714752870173892994651D+00
    w(7) =  0.269266719309996355091226921569D+00
    w(8) =  0.219086362515982043995534934228D+00
    w(9) =  0.149451349150580593145776339658D+00
    w(10) = 0.666713443086881375935688098933D-01

  else if ( n == 11 ) then

    x(1) =  - 0.978228658146056992803938001123D+00
    x(2) =  - 0.887062599768095299075157769304D+00
    x(3) =  - 0.730152005574049324093416252031D+00
    x(4) =  - 0.519096129206811815925725669459D+00
    x(5) =  - 0.269543155952344972331531985401D+00
    x(6) =    0.0D+00
    x(7) =    0.269543155952344972331531985401D+00
    x(8) =    0.519096129206811815925725669459D+00
    x(9) =    0.730152005574049324093416252031D+00
    x(10) =   0.887062599768095299075157769304D+00
    x(11) =   0.978228658146056992803938001123D+00

    w(1) =  0.556685671161736664827537204425D-01
    w(2) =  0.125580369464904624634694299224D+00
    w(3) =  0.186290210927734251426097641432D+00
    w(4) =  0.233193764591990479918523704843D+00
    w(5) =  0.262804544510246662180688869891D+00
    w(6) =  0.272925086777900630714483528336D+00
    w(7) =  0.262804544510246662180688869891D+00
    w(8) =  0.233193764591990479918523704843D+00
    w(9) =  0.186290210927734251426097641432D+00
    w(10) = 0.125580369464904624634694299224D+00
    w(11) = 0.556685671161736664827537204425D-01

  else if ( n == 12 ) then

    x(1) =  - 0.981560634246719250690549090149D+00
    x(2) =  - 0.904117256370474856678465866119D+00
    x(3) =  - 0.769902674194304687036893833213D+00
    x(4) =  - 0.587317954286617447296702418941D+00
    x(5) =  - 0.367831498998180193752691536644D+00
    x(6) =  - 0.125233408511468915472441369464D+00
    x(7) =    0.125233408511468915472441369464D+00
    x(8) =    0.367831498998180193752691536644D+00
    x(9) =    0.587317954286617447296702418941D+00
    x(10) =   0.769902674194304687036893833213D+00
    x(11) =   0.904117256370474856678465866119D+00
    x(12) =   0.981560634246719250690549090149D+00

    w(1) =  0.471753363865118271946159614850D-01
    w(2) =  0.106939325995318430960254718194D+00
    w(3) =  0.160078328543346226334652529543D+00
    w(4) =  0.203167426723065921749064455810D+00
    w(5) =  0.233492536538354808760849898925D+00
    w(6) =  0.249147045813402785000562436043D+00
    w(7) =  0.249147045813402785000562436043D+00
    w(8) =  0.233492536538354808760849898925D+00
    w(9) =  0.203167426723065921749064455810D+00
    w(10) = 0.160078328543346226334652529543D+00
    w(11) = 0.106939325995318430960254718194D+00
    w(12) = 0.471753363865118271946159614850D-01

  else if ( n == 13 ) then

    x(1) =  - 0.984183054718588149472829448807D+00
    x(2) =  - 0.917598399222977965206547836501D+00
    x(3) =  - 0.801578090733309912794206489583D+00
    x(4) =  - 0.642349339440340220643984606996D+00
    x(5) =  - 0.448492751036446852877912852128D+00
    x(6) =  - 0.230458315955134794065528121098D+00
    x(7) =    0.0D+00
    x(8) =    0.230458315955134794065528121098D+00
    x(9) =    0.448492751036446852877912852128D+00
    x(10) =   0.642349339440340220643984606996D+00
    x(11) =   0.801578090733309912794206489583D+00
    x(12) =   0.917598399222977965206547836501D+00
    x(13) =   0.984183054718588149472829448807D+00

    w(1) =  0.404840047653158795200215922010D-01
    w(2) =  0.921214998377284479144217759538D-01
    w(3) =  0.138873510219787238463601776869D+00
    w(4) =  0.178145980761945738280046691996D+00
    w(5) =  0.207816047536888502312523219306D+00
    w(6) =  0.226283180262897238412090186040D+00
    w(7) =  0.232551553230873910194589515269D+00
    w(8) =  0.226283180262897238412090186040D+00
    w(9) =  0.207816047536888502312523219306D+00
    w(10) = 0.178145980761945738280046691996D+00
    w(11) = 0.138873510219787238463601776869D+00
    w(12) = 0.921214998377284479144217759538D-01
    w(13) = 0.404840047653158795200215922010D-01

  else if ( n == 14 ) then

    x(1) =  - 0.986283808696812338841597266704D+00
    x(2) =  - 0.928434883663573517336391139378D+00
    x(3) =  - 0.827201315069764993189794742650D+00
    x(4) =  - 0.687292904811685470148019803019D+00
    x(5) =  - 0.515248636358154091965290718551D+00
    x(6) =  - 0.319112368927889760435671824168D+00
    x(7) =  - 0.108054948707343662066244650220D+00
    x(8) =    0.108054948707343662066244650220D+00
    x(9) =    0.319112368927889760435671824168D+00
    x(10) =   0.515248636358154091965290718551D+00
    x(11) =   0.687292904811685470148019803019D+00
    x(12) =   0.827201315069764993189794742650D+00
    x(13) =   0.928434883663573517336391139378D+00
    x(14) =   0.986283808696812338841597266704D+00

    w(1) =  0.351194603317518630318328761382D-01
    w(2) =  0.801580871597602098056332770629D-01
    w(3) =  0.121518570687903184689414809072D+00
    w(4) =  0.157203167158193534569601938624D+00
    w(5) =  0.185538397477937813741716590125D+00
    w(6) =  0.205198463721295603965924065661D+00
    w(7) =  0.215263853463157790195876443316D+00
    w(8) =  0.215263853463157790195876443316D+00
    w(9) =  0.205198463721295603965924065661D+00
    w(10) = 0.185538397477937813741716590125D+00
    w(11) = 0.157203167158193534569601938624D+00
    w(12) = 0.121518570687903184689414809072D+00
    w(13) = 0.801580871597602098056332770629D-01
    w(14) = 0.351194603317518630318328761382D-01

  else if ( n == 15 ) then

    x(1) =  - 0.987992518020485428489565718587D+00
    x(2) =  - 0.937273392400705904307758947710D+00
    x(3) =  - 0.848206583410427216200648320774D+00
    x(4) =  - 0.724417731360170047416186054614D+00
    x(5) =  - 0.570972172608538847537226737254D+00
    x(6) =  - 0.394151347077563369897207370981D+00
    x(7) =  - 0.201194093997434522300628303395D+00
    x(8) =    0.0D+00
    x(9) =    0.201194093997434522300628303395D+00
    x(10) =   0.394151347077563369897207370981D+00
    x(11) =   0.570972172608538847537226737254D+00
    x(12) =   0.724417731360170047416186054614D+00
    x(13) =   0.848206583410427216200648320774D+00
    x(14) =   0.937273392400705904307758947710D+00
    x(15) =   0.987992518020485428489565718587D+00

    w(1) =  0.307532419961172683546283935772D-01
    w(2) =  0.703660474881081247092674164507D-01
    w(3) =  0.107159220467171935011869546686D+00
    w(4) =  0.139570677926154314447804794511D+00
    w(5) =  0.166269205816993933553200860481D+00
    w(6) =  0.186161000015562211026800561866D+00
    w(7) =  0.198431485327111576456118326444D+00
    w(8) =  0.202578241925561272880620199968D+00
    w(9) =  0.198431485327111576456118326444D+00
    w(10) = 0.186161000015562211026800561866D+00
    w(11) = 0.166269205816993933553200860481D+00
    w(12) = 0.139570677926154314447804794511D+00
    w(13) = 0.107159220467171935011869546686D+00
    w(14) = 0.703660474881081247092674164507D-01
    w(15) = 0.307532419961172683546283935772D-01

  else if ( n == 16 ) then

    x(1) =  - 0.989400934991649932596154173450D+00
    x(2) =  - 0.944575023073232576077988415535D+00
    x(3) =  - 0.865631202387831743880467897712D+00
    x(4) =  - 0.755404408355003033895101194847D+00
    x(5) =  - 0.617876244402643748446671764049D+00
    x(6) =  - 0.458016777657227386342419442984D+00
    x(7) =  - 0.281603550779258913230460501460D+00
    x(8) =  - 0.950125098376374401853193354250D-01
    x(9) =    0.950125098376374401853193354250D-01
    x(10) =   0.281603550779258913230460501460D+00
    x(11) =   0.458016777657227386342419442984D+00
    x(12) =   0.617876244402643748446671764049D+00
    x(13) =   0.755404408355003033895101194847D+00
    x(14) =   0.865631202387831743880467897712D+00
    x(15) =   0.944575023073232576077988415535D+00
    x(16) =   0.989400934991649932596154173450D+00

    w(1) =  0.271524594117540948517805724560D-01
    w(2) =  0.622535239386478928628438369944D-01
    w(3) =  0.951585116824927848099251076022D-01
    w(4) =  0.124628971255533872052476282192D+00
    w(5) =  0.149595988816576732081501730547D+00
    w(6) =  0.169156519395002538189312079030D+00
    w(7) =  0.182603415044923588866763667969D+00
    w(8) =  0.189450610455068496285396723208D+00
    w(9) =  0.189450610455068496285396723208D+00
    w(10) = 0.182603415044923588866763667969D+00
    w(11) = 0.169156519395002538189312079030D+00
    w(12) = 0.149595988816576732081501730547D+00
    w(13) = 0.124628971255533872052476282192D+00
    w(14) = 0.951585116824927848099251076022D-01
    w(15) = 0.622535239386478928628438369944D-01
    w(16) = 0.271524594117540948517805724560D-01

  else if ( n == 17 ) then

    x(1) =  - 0.990575475314417335675434019941D+00
    x(2) =  - 0.950675521768767761222716957896D+00
    x(3) =  - 0.880239153726985902122955694488D+00
    x(4) =  - 0.781514003896801406925230055520D+00
    x(5) =  - 0.657671159216690765850302216643D+00
    x(6) =  - 0.512690537086476967886246568630D+00
    x(7) =  - 0.351231763453876315297185517095D+00
    x(8) =  - 0.178484181495847855850677493654D+00
    x(9) =    0.0D+00
    x(10) =   0.178484181495847855850677493654D+00
    x(11) =   0.351231763453876315297185517095D+00
    x(12) =   0.512690537086476967886246568630D+00
    x(13) =   0.657671159216690765850302216643D+00
    x(14) =   0.781514003896801406925230055520D+00
    x(15) =   0.880239153726985902122955694488D+00
    x(16) =   0.950675521768767761222716957896D+00
    x(17) =   0.990575475314417335675434019941D+00

    w(1) =  0.241483028685479319601100262876D-01
    w(2) =  0.554595293739872011294401653582D-01
    w(3) =  0.850361483171791808835353701911D-01
    w(4) =  0.111883847193403971094788385626D+00
    w(5) =  0.135136368468525473286319981702D+00
    w(6) =  0.154045761076810288081431594802D+00
    w(7) =  0.168004102156450044509970663788D+00
    w(8) =  0.176562705366992646325270990113D+00
    w(9) =  0.179446470356206525458265644262D+00
    w(10) = 0.176562705366992646325270990113D+00
    w(11) = 0.168004102156450044509970663788D+00
    w(12) = 0.154045761076810288081431594802D+00
    w(13) = 0.135136368468525473286319981702D+00
    w(14) = 0.111883847193403971094788385626D+00
    w(15) = 0.850361483171791808835353701911D-01
    w(16) = 0.554595293739872011294401653582D-01
    w(17) = 0.241483028685479319601100262876D-01

  else if ( n == 18 ) then

    x(1) =  - 0.991565168420930946730016004706D+00
    x(2) =  - 0.955823949571397755181195892930D+00
    x(3) =  - 0.892602466497555739206060591127D+00
    x(4) =  - 0.803704958972523115682417455015D+00
    x(5) =  - 0.691687043060353207874891081289D+00
    x(6) =  - 0.559770831073947534607871548525D+00
    x(7) =  - 0.411751161462842646035931793833D+00
    x(8) =  - 0.251886225691505509588972854878D+00
    x(9) =  - 0.847750130417353012422618529358D-01
    x(10) =   0.847750130417353012422618529358D-01
    x(11) =   0.251886225691505509588972854878D+00
    x(12) =   0.411751161462842646035931793833D+00
    x(13) =   0.559770831073947534607871548525D+00
    x(14) =   0.691687043060353207874891081289D+00
    x(15) =   0.803704958972523115682417455015D+00
    x(16) =   0.892602466497555739206060591127D+00
    x(17) =   0.955823949571397755181195892930D+00
    x(18) =   0.991565168420930946730016004706D+00

    w(1) =  0.216160135264833103133427102665D-01
    w(2) =  0.497145488949697964533349462026D-01
    w(3) =  0.764257302548890565291296776166D-01
    w(4) =  0.100942044106287165562813984925D+00
    w(5) =  0.122555206711478460184519126800D+00
    w(6) =  0.140642914670650651204731303752D+00
    w(7) =  0.154684675126265244925418003836D+00
    w(8) =  0.164276483745832722986053776466D+00
    w(9) =  0.169142382963143591840656470135D+00
    w(10) = 0.169142382963143591840656470135D+00
    w(11) = 0.164276483745832722986053776466D+00
    w(12) = 0.154684675126265244925418003836D+00
    w(13) = 0.140642914670650651204731303752D+00
    w(14) = 0.122555206711478460184519126800D+00
    w(15) = 0.100942044106287165562813984925D+00
    w(16) = 0.764257302548890565291296776166D-01
    w(17) = 0.497145488949697964533349462026D-01
    w(18) = 0.216160135264833103133427102665D-01

  else if ( n == 19 ) then

    x(1) =  - 0.992406843843584403189017670253D+00
    x(2) =  - 0.960208152134830030852778840688D+00
    x(3) =  - 0.903155903614817901642660928532D+00
    x(4) =  - 0.822714656537142824978922486713D+00
    x(5) =  - 0.720966177335229378617095860824D+00
    x(6) =  - 0.600545304661681023469638164946D+00
    x(7) =  - 0.464570741375960945717267148104D+00
    x(8) =  - 0.316564099963629831990117328850D+00
    x(9) =  - 0.160358645640225375868096115741D+00
    x(10) =   0.0D+00
    x(11) =   0.160358645640225375868096115741D+00
    x(12) =   0.316564099963629831990117328850D+00
    x(13) =   0.464570741375960945717267148104D+00
    x(14) =   0.600545304661681023469638164946D+00
    x(15) =   0.720966177335229378617095860824D+00
    x(16) =   0.822714656537142824978922486713D+00
    x(17) =   0.903155903614817901642660928532D+00
    x(18) =   0.960208152134830030852778840688D+00
    x(19) =   0.992406843843584403189017670253D+00

    w(1) =  0.194617882297264770363120414644D-01
    w(2) =  0.448142267656996003328381574020D-01
    w(3) =  0.690445427376412265807082580060D-01
    w(4) =  0.914900216224499994644620941238D-01
    w(5) =  0.111566645547333994716023901682D+00
    w(6) =  0.128753962539336227675515784857D+00
    w(7) =  0.142606702173606611775746109442D+00
    w(8) =  0.152766042065859666778855400898D+00
    w(9) =  0.158968843393954347649956439465D+00
    w(10) = 0.161054449848783695979163625321D+00
    w(11) = 0.158968843393954347649956439465D+00
    w(12) = 0.152766042065859666778855400898D+00
    w(13) = 0.142606702173606611775746109442D+00
    w(14) = 0.128753962539336227675515784857D+00
    w(15) = 0.111566645547333994716023901682D+00
    w(16) = 0.914900216224499994644620941238D-01
    w(17) = 0.690445427376412265807082580060D-01
    w(18) = 0.448142267656996003328381574020D-01
    w(19) = 0.194617882297264770363120414644D-01

  else if ( n == 20 ) then

    x(1) =  - 0.993128599185094924786122388471D+00
    x(2) =  - 0.963971927277913791267666131197D+00
    x(3) =  - 0.912234428251325905867752441203D+00
    x(4) =  - 0.839116971822218823394529061702D+00
    x(5) =  - 0.746331906460150792614305070356D+00
    x(6) =  - 0.636053680726515025452836696226D+00
    x(7) =  - 0.510867001950827098004364050955D+00
    x(8) =  - 0.373706088715419560672548177025D+00
    x(9) =  - 0.227785851141645078080496195369D+00
    x(10) = - 0.765265211334973337546404093988D-01
    x(11) =   0.765265211334973337546404093988D-01
    x(12) =   0.227785851141645078080496195369D+00
    x(13) =   0.373706088715419560672548177025D+00
    x(14) =   0.510867001950827098004364050955D+00
    x(15) =   0.636053680726515025452836696226D+00
    x(16) =   0.746331906460150792614305070356D+00
    x(17) =   0.839116971822218823394529061702D+00
    x(18) =   0.912234428251325905867752441203D+00
    x(19) =   0.963971927277913791267666131197D+00
    x(20) =   0.993128599185094924786122388471D+00

    w(1) =  0.176140071391521183118619623519D-01
    w(2) =  0.406014298003869413310399522749D-01
    w(3) =  0.626720483341090635695065351870D-01
    w(4) =  0.832767415767047487247581432220D-01
    w(5) =  0.101930119817240435036750135480D+00
    w(6) =  0.118194531961518417312377377711D+00
    w(7) =  0.131688638449176626898494499748D+00
    w(8) =  0.142096109318382051329298325067D+00
    w(9) =  0.149172986472603746787828737002D+00
    w(10) = 0.152753387130725850698084331955D+00
    w(11) = 0.152753387130725850698084331955D+00
    w(12) = 0.149172986472603746787828737002D+00
    w(13) = 0.142096109318382051329298325067D+00
    w(14) = 0.131688638449176626898494499748D+00
    w(15) = 0.118194531961518417312377377711D+00
    w(16) = 0.101930119817240435036750135480D+00
    w(17) = 0.832767415767047487247581432220D-01
    w(18) = 0.626720483341090635695065351870D-01
    w(19) = 0.406014298003869413310399522749D-01
    w(20) = 0.176140071391521183118619623519D-01

  else if ( n == 21 ) then

    x( 1) =  -0.9937521706203896D+00
    x( 2) =  -0.9672268385663063D+00
    x( 3) =  -0.9200993341504008D+00
    x( 4) =  -0.8533633645833173D+00
    x( 5) =  -0.7684399634756779D+00
    x( 6) =  -0.6671388041974123D+00
    x( 7) =  -0.5516188358872198D+00
    x( 8) =  -0.4243421202074388D+00
    x( 9) =  -0.2880213168024011D+00
    x(10) =  -0.1455618541608951D+00
    x(11) =   0.0000000000000000D+00
    x(12) =   0.1455618541608951D+00
    x(13) =   0.2880213168024011D+00
    x(14) =   0.4243421202074388D+00
    x(15) =   0.5516188358872198D+00
    x(16) =   0.6671388041974123D+00
    x(17) =   0.7684399634756779D+00
    x(18) =   0.8533633645833173D+00
    x(19) =   0.9200993341504008D+00
    x(20) =   0.9672268385663063D+00
    x(21) =   0.9937521706203896D+00 
   
    w( 1) =   0.1601722825777420D-01
    w( 2) =   0.3695378977085242D-01
    w( 3) =   0.5713442542685715D-01
    w( 4) =   0.7610011362837928D-01
    w( 5) =   0.9344442345603393D-01
    w( 6) =   0.1087972991671484D+00
    w( 7) =   0.1218314160537285D+00
    w( 8) =   0.1322689386333373D+00
    w( 9) =   0.1398873947910731D+00
    w(10) =   0.1445244039899700D+00
    w(11) =   0.1460811336496904D+00
    w(12) =   0.1445244039899700D+00
    w(13) =   0.1398873947910731D+00
    w(14) =   0.1322689386333373D+00
    w(15) =   0.1218314160537285D+00
    w(16) =   0.1087972991671484D+00
    w(17) =   0.9344442345603393D-01
    w(18) =   0.7610011362837928D-01
    w(19) =   0.5713442542685715D-01
    w(20) =   0.3695378977085242D-01
    w(21) =   0.1601722825777420D-01

  else if ( n == 22 ) then

    x( 1) =  -0.9942945854823994D+00
    x( 2) =  -0.9700604978354287D+00
    x( 3) =  -0.9269567721871740D+00
    x( 4) =  -0.8658125777203002D+00
    x( 5) =  -0.7878168059792081D+00
    x( 6) =  -0.6944872631866827D+00
    x( 7) =  -0.5876404035069116D+00
    x( 8) =  -0.4693558379867570D+00
    x( 9) =  -0.3419358208920842D+00
    x(10) =  -0.2078604266882213D+00
    x(11) =  -0.6973927331972223D-01
    x(12) =   0.6973927331972223D-01
    x(13) =   0.2078604266882213D+00
    x(14) =   0.3419358208920842D+00
    x(15) =   0.4693558379867570D+00
    x(16) =   0.5876404035069116D+00
    x(17) =   0.6944872631866827D+00
    x(18) =   0.7878168059792081D+00
    x(19) =   0.8658125777203002D+00
    x(20) =   0.9269567721871740D+00
    x(21) =   0.9700604978354287D+00
    x(22) =   0.9942945854823994D+00
 
    w( 1) =   0.1462799529827203D-01
    w( 2) =   0.3377490158481413D-01
    w( 3) =   0.5229333515268327D-01
    w( 4) =   0.6979646842452038D-01
    w( 5) =   0.8594160621706777D-01
    w( 6) =   0.1004141444428809D+00
    w( 7) =   0.1129322960805392D+00
    w( 8) =   0.1232523768105124D+00
    w( 9) =   0.1311735047870623D+00
    w(10) =   0.1365414983460152D+00
    w(11) =   0.1392518728556321D+00
    w(12) =   0.1392518728556321D+00
    w(13) =   0.1365414983460152D+00
    w(14) =   0.1311735047870623D+00
    w(15) =   0.1232523768105124D+00
    w(16) =   0.1129322960805392D+00
    w(17) =   0.1004141444428809D+00
    w(18) =   0.8594160621706777D-01
    w(19) =   0.6979646842452038D-01
    w(20) =   0.5229333515268327D-01
    w(21) =   0.3377490158481413D-01
    w(22) =   0.1462799529827203D-01

  else if ( n == 23 ) then

    x( 1) =  -0.9947693349975522D+00
    x( 2) =  -0.9725424712181152D+00
    x( 3) =  -0.9329710868260161D+00
    x( 4) =  -0.8767523582704416D+00
    x( 5) =  -0.8048884016188399D+00
    x( 6) =  -0.7186613631319502D+00
    x( 7) =  -0.6196098757636461D+00
    x( 8) =  -0.5095014778460075D+00
    x( 9) =  -0.3903010380302908D+00
    x(10) =  -0.2641356809703449D+00
    x(11) =  -0.1332568242984661D+00
    x(12) =   0.0000000000000000D+00
    x(13) =   0.1332568242984661D+00
    x(14) =   0.2641356809703449D+00
    x(15) =   0.3903010380302908D+00
    x(16) =   0.5095014778460075D+00
    x(17) =   0.6196098757636461D+00
    x(18) =   0.7186613631319502D+00
    x(19) =   0.8048884016188399D+00
    x(20) =   0.8767523582704416D+00
    x(21) =   0.9329710868260161D+00
    x(22) =   0.9725424712181152D+00
    x(23) =   0.9947693349975522D+00
 
    w( 1) =   0.1341185948714167D-01
    w( 2) =   0.3098800585697944D-01
    w( 3) =   0.4803767173108464D-01
    w( 4) =   0.6423242140852586D-01
    w( 5) =   0.7928141177671895D-01
    w( 6) =   0.9291576606003514D-01
    w( 7) =   0.1048920914645414D+00
    w( 8) =   0.1149966402224114D+00
    w( 9) =   0.1230490843067295D+00
    w(10) =   0.1289057221880822D+00
    w(11) =   0.1324620394046967D+00
    w(12) =   0.1336545721861062D+00
    w(13) =   0.1324620394046967D+00
    w(14) =   0.1289057221880822D+00
    w(15) =   0.1230490843067295D+00
    w(16) =   0.1149966402224114D+00
    w(17) =   0.1048920914645414D+00
    w(18) =   0.9291576606003514D-01
    w(19) =   0.7928141177671895D-01
    w(20) =   0.6423242140852586D-01
    w(21) =   0.4803767173108464D-01
    w(22) =   0.3098800585697944D-01
    w(23) =   0.1341185948714167D-01

  else if ( n == 24 ) then

    x( 1) =  -0.9951872199970213D+00    
    x( 2) =  -0.9747285559713095D+00    
    x( 3) =  -0.9382745520027327D+00    
    x( 4) =  -0.8864155270044011D+00    
    x( 5) =  -0.8200019859739029D+00    
    x( 6) =  -0.7401241915785544D+00    
    x( 7) =  -0.6480936519369755D+00    
    x( 8) =  -0.5454214713888396D+00    
    x( 9) =  -0.4337935076260451D+00    
    x(10) =  -0.3150426796961634D+00    
    x(11) =  -0.1911188674736163D+00    
    x(12) =  -0.6405689286260562D-01
    x(13) =   0.6405689286260562D-01
    x(14) =   0.1911188674736163D+00    
    x(15) =   0.3150426796961634D+00    
    x(16) =   0.4337935076260451D+00    
    x(17) =   0.5454214713888396D+00    
    x(18) =   0.6480936519369755D+00    
    x(19) =   0.7401241915785544D+00    
    x(20) =   0.8200019859739029D+00    
    x(21) =   0.8864155270044011D+00    
    x(22) =   0.9382745520027327D+00    
    x(23) =   0.9747285559713095D+00    
    x(24) =   0.9951872199970213D+00    
 
    w( 1) =   0.1234122979998730D-01
    w( 2) =   0.2853138862893375D-01
    w( 3) =   0.4427743881741982D-01
    w( 4) =   0.5929858491543672D-01
    w( 5) =   0.7334648141108031D-01
    w( 6) =   0.8619016153195320D-01
    w( 7) =   0.9761865210411380D-01
    w( 8) =   0.1074442701159656D+00    
    w( 9) =   0.1155056680537256D+00    
    w(10) =   0.1216704729278035D+00    
    w(11) =   0.1258374563468283D+00    
    w(12) =   0.1279381953467521D+00    
    w(13) =   0.1279381953467521D+00    
    w(14) =   0.1258374563468283D+00    
    w(15) =   0.1216704729278035D+00    
    w(16) =   0.1155056680537256D+00   
    w(17) =   0.1074442701159656D+00    
    w(18) =   0.9761865210411380D-01
    w(19) =   0.8619016153195320D-01
    w(20) =   0.7334648141108031D-01
    w(21) =   0.5929858491543672D-01
    w(22) =   0.4427743881741982D-01
    w(23) =   0.2853138862893375D-01
    w(24) =   0.1234122979998730D-01

  else if ( n == 25 ) then

    x( 1) =  -0.9955569697904981D+00    
    x( 2) =  -0.9766639214595175D+00    
    x( 3) =  -0.9429745712289743D+00    
    x( 4) =  -0.8949919978782754D+00    
    x( 5) =  -0.8334426287608340D+00    
    x( 6) =  -0.7592592630373577D+00    
    x( 7) =  -0.6735663684734684D+00    
    x( 8) =  -0.5776629302412229D+00    
    x( 9) =  -0.4730027314457150D+00    
    x(10) =  -0.3611723058093879D+00    
    x(11) =  -0.2438668837209884D+00    
    x(12) =  -0.1228646926107104D+00    
    x(13) =   0.0000000000000000D+00    
    x(14) =   0.1228646926107104D+00  
    x(15) =   0.2438668837209884D+00    
    x(16) =   0.3611723058093879D+00    
    x(17) =   0.4730027314457150D+00    
    x(18) =   0.5776629302412229D+00    
    x(19) =   0.6735663684734684D+00    
    x(20) =   0.7592592630373577D+00    
    x(21) =   0.8334426287608340D+00    
    x(22) =   0.8949919978782754D+00    
    x(23) =   0.9429745712289743D+00    
    x(24) =   0.9766639214595175D+00    
    x(25) =   0.9955569697904981D+00    
 
    w( 1) =   0.1139379850102617D-01
    w( 2) =   0.2635498661503214D-01
    w( 3) =   0.4093915670130639D-01
    w( 4) =   0.5490469597583517D-01
    w( 5) =   0.6803833381235694D-01
    w( 6) =   0.8014070033500101D-01
    w( 7) =   0.9102826198296370D-01
    w( 8) =   0.1005359490670506D+00    
    w( 9) =   0.1085196244742637D+00    
    w(10) =   0.1148582591457116D+00    
    w(11) =   0.1194557635357847D+00    
    w(12) =   0.1222424429903101D+00    
    w(13) =   0.1231760537267154D+00    
    w(14) =   0.1222424429903101D+00    
    w(15) =   0.1194557635357847D+00    
    w(16) =   0.1148582591457116D+00    
    w(17) =   0.1085196244742637D+00    
    w(18) =   0.1005359490670506D+00    
    w(19) =   0.9102826198296370D-01
    w(20) =   0.8014070033500101D-01
    w(21) =   0.6803833381235694D-01
    w(22) =   0.5490469597583517D-01
    w(23) =   0.4093915670130639D-01
    w(24) =   0.2635498661503214D-01
    w(25) =   0.1139379850102617D-01

  else if ( n == 26 ) then

    x( 1) =  -0.9958857011456169D+00    
    x( 2) =  -0.9783854459564710D+00    
    x( 3) =  -0.9471590666617142D+00    
    x( 4) =  -0.9026378619843071D+00    
    x( 5) =  -0.8454459427884981D+00    
    x( 6) =  -0.7763859488206789D+00    
    x( 7) =  -0.6964272604199573D+00    
    x( 8) =  -0.6066922930176181D+00    
    x( 9) =  -0.5084407148245057D+00    
    x(10) =  -0.4030517551234863D+00    
    x(11) =  -0.2920048394859569D+00    
    x(12) =  -0.1768588203568902D+00    
    x(13) =  -0.5923009342931320D-01
    x(14) =   0.5923009342931320D-01
    x(15) =   0.1768588203568902D+00    
    x(16) =   0.2920048394859569D+00    
    x(17) =   0.4030517551234863D+00    
    x(18) =   0.5084407148245057D+00    
    x(19) =   0.6066922930176181D+00    
    x(20) =   0.6964272604199573D+00    
    x(21) =   0.7763859488206789D+00    
    x(22) =   0.8454459427884981D+00    
    x(23) =   0.9026378619843071D+00    
    x(24) =   0.9471590666617142D+00    
    x(25) =   0.9783854459564710D+00    
    x(26) =   0.9958857011456169D+00    
 
    w( 1) =   0.1055137261734304D-01
    w( 2) =   0.2441785109263173D-01
    w( 3) =   0.3796238329436282D-01
    w( 4) =   0.5097582529714782D-01
    w( 5) =   0.6327404632957484D-01
    w( 6) =   0.7468414976565967D-01
    w( 7) =   0.8504589431348521D-01
    w( 8) =   0.9421380035591416D-01
    w( 9) =   0.1020591610944255D+00    
    w(10) =   0.1084718405285765D+00    
    w(11) =   0.1133618165463197D+00    
    w(12) =   0.1166604434852967D+00    
    w(13) =   0.1183214152792622D+00    
    w(14) =   0.1183214152792622D+00    
    w(15) =   0.1166604434852967D+00    
    w(16) =   0.1133618165463197D+00    
    w(17) =   0.1084718405285765D+00    
    w(18) =   0.1020591610944255D+00    
    w(19) =   0.9421380035591416D-01
    w(20) =   0.8504589431348521D-01
    w(21) =   0.7468414976565967D-01
    w(22) =   0.6327404632957484D-01
    w(23) =   0.5097582529714782D-01
    w(24) =   0.3796238329436282D-01
    w(25) =   0.2441785109263173D-01
    w(26) =   0.1055137261734304D-01

  else if ( n == 27 ) then

    x( 1) =  -0.9961792628889886D+00    
    x( 2) =  -0.9799234759615012D+00    
    x( 3) =  -0.9509005578147051D+00    
    x( 4) =  -0.9094823206774911D+00    
    x( 5) =  -0.8562079080182945D+00    
    x( 6) =  -0.7917716390705082D+00    
    x( 7) =  -0.7170134737394237D+00    
    x( 8) =  -0.6329079719464952D+00    
    x( 9) =  -0.5405515645794569D+00    
    x(10) =  -0.4411482517500269D+00    
    x(11) =  -0.3359939036385089D+00    
    x(12) =  -0.2264593654395369D+00    
    x(13) =  -0.1139725856095300D+00    
    x(14) =   0.0000000000000000D+00    
    x(15) =   0.1139725856095300D+00    
    x(16) =   0.2264593654395369D+00    
    x(17) =   0.3359939036385089D+00    
    x(18) =   0.4411482517500269D+00    
    x(19) =   0.5405515645794569D+00    
    x(20) =   0.6329079719464952D+00    
    x(21) =   0.7170134737394237D+00    
    x(22) =   0.7917716390705082D+00    
    x(23) =   0.8562079080182945D+00    
    x(24) =   0.9094823206774911D+00    
    x(25) =   0.9509005578147051D+00    
    x(26) =   0.9799234759615012D+00    
    x(27) =   0.9961792628889886D+00    
 
    w( 1) =   0.9798996051294232D-02
    w( 2) =   0.2268623159618062D-01
    w( 3) =   0.3529705375741969D-01
    w( 4) =   0.4744941252061504D-01
    w( 5) =   0.5898353685983366D-01
    w( 6) =   0.6974882376624561D-01
    w( 7) =   0.7960486777305781D-01
    w( 8) =   0.8842315854375689D-01
    w( 9) =   0.9608872737002842D-01
    w(10) =   0.1025016378177459D+00    
    w(11) =   0.1075782857885332D+00    
    w(12) =   0.1112524883568452D+00    
    w(13) =   0.1134763461089651D+00    
    w(14) =   0.1142208673789570D+00    
    w(15) =   0.1134763461089651D+00    
    w(16) =   0.1112524883568452D+00    
    w(17) =   0.1075782857885332D+00    
    w(18) =   0.1025016378177459D+00    
    w(19) =   0.9608872737002842D-01
    w(20) =   0.8842315854375689D-01
    w(21) =   0.7960486777305781D-01
    w(22) =   0.6974882376624561D-01
    w(23) =   0.5898353685983366D-01
    w(24) =   0.4744941252061504D-01
    w(25) =   0.3529705375741969D-01
    w(26) =   0.2268623159618062D-01
    w(27) =   0.9798996051294232D-02

  else if ( n == 28 ) then

    x( 1) =  -0.9964424975739544D+00    
    x( 2) =  -0.9813031653708728D+00    
    x( 3) =  -0.9542592806289382D+00    
    x( 4) =  -0.9156330263921321D+00    
    x( 5) =  -0.8658925225743951D+00    
    x( 6) =  -0.8056413709171791D+00    
    x( 7) =  -0.7356108780136318D+00    
    x( 8) =  -0.6566510940388650D+00    
    x( 9) =  -0.5697204718114017D+00    
    x(10) =  -0.4758742249551183D+00    
    x(11) =  -0.3762515160890787D+00    
    x(12) =  -0.2720616276351780D+00    
    x(13) =  -0.1645692821333808D+00    
    x(14) =  -0.5507928988403427D-01
    x(15) =   0.5507928988403427D-01
    x(16) =   0.1645692821333808D+00    
    x(17) =   0.2720616276351780D+00    
    x(18) =   0.3762515160890787D+00    
    x(19) =   0.4758742249551183D+00    
    x(20) =   0.5697204718114017D+00    
    x(21) =   0.6566510940388650D+00    
    x(22) =   0.7356108780136318D+00    
    x(23) =   0.8056413709171791D+00    
    x(24) =   0.8658925225743951D+00    
    x(25) =   0.9156330263921321D+00    
    x(26) =   0.9542592806289382D+00    
    x(27) =   0.9813031653708728D+00    
    x(28) =   0.9964424975739544D+00    
 
    w( 1) =   0.9124282593094672D-02
    w( 2) =   0.2113211259277118D-01
    w( 3) =   0.3290142778230441D-01
    w( 4) =   0.4427293475900429D-01
    w( 5) =   0.5510734567571667D-01
    w( 6) =   0.6527292396699959D-01
    w( 7) =   0.7464621423456877D-01
    w( 8) =   0.8311341722890127D-01
    w( 9) =   0.9057174439303289D-01
    w(10) =   0.9693065799792999D-01
    w(11) =   0.1021129675780608D+00    
    w(12) =   0.1060557659228464D+00    
    w(13) =   0.1087111922582942D+00    
    w(14) =   0.1100470130164752D+00    
    w(15) =   0.1100470130164752D+00    
    w(16) =   0.1087111922582942D+00    
    w(17) =   0.1060557659228464D+00    
    w(18) =   0.1021129675780608D+00   
    w(19) =   0.9693065799792999D-01
    w(20) =   0.9057174439303289D-01
    w(21) =   0.8311341722890127D-01
    w(22) =   0.7464621423456877D-01
    w(23) =   0.6527292396699959D-01
    w(24) =   0.5510734567571667D-01
    w(25) =   0.4427293475900429D-01
    w(26) =   0.3290142778230441D-01
    w(27) =   0.2113211259277118D-01
    w(28) =   0.9124282593094672D-02

  else if ( n == 29 ) then

    x( 1) =  -0.9966794422605966D+00    
    x( 2) =  -0.9825455052614132D+00    
    x( 3) =  -0.9572855957780877D+00    
    x( 4) =  -0.9211802329530588D+00    
    x( 5) =  -0.8746378049201028D+00    
    x( 6) =  -0.8181854876152524D+00    
    x( 7) =  -0.7524628517344771D+00    
    x( 8) =  -0.6782145376026865D+00    
    x( 9) =  -0.5962817971382278D+00    
    x(10) =  -0.5075929551242276D+00    
    x(11) =  -0.4131528881740087D+00    
    x(12) =  -0.3140316378676399D+00    
    x(13) =  -0.2113522861660011D+00    
    x(14) =  -0.1062782301326792D+00    
    x(15) =   0.0000000000000000D+00    
    x(16) =   0.1062782301326792D+00    
    x(17) =   0.2113522861660011D+00    
    x(18) =   0.3140316378676399D+00    
    x(19) =   0.4131528881740087D+00    
    x(20) =   0.5075929551242276D+00    
    x(21) =   0.5962817971382278D+00    
    x(22) =   0.6782145376026865D+00    
    x(23) =   0.7524628517344771D+00    
    x(24) =   0.8181854876152524D+00    
    x(25) =   0.8746378049201028D+00    
    x(26) =   0.9211802329530588D+00    
    x(27) =   0.9572855957780877D+00    
    x(28) =   0.9825455052614132D+00    
    x(29) =   0.9966794422605966D+00    
 
    w( 1) =   0.8516903878746365D-02
    w( 2) =   0.1973208505612276D-01
    w( 3) =   0.3074049220209360D-01
    w( 4) =   0.4140206251868281D-01
    w( 5) =   0.5159482690249799D-01
    w( 6) =   0.6120309065707916D-01
    w( 7) =   0.7011793325505125D-01
    w( 8) =   0.7823832713576385D-01
    w( 9) =   0.8547225736617248D-01
    w(10) =   0.9173775713925882D-01
    w(11) =   0.9696383409440862D-01
    w(12) =   0.1010912737599150D+00    
    w(13) =   0.1040733100777293D+00    
    w(14) =   0.1058761550973210D+00    
    w(15) =   0.1064793817183143D+00    
    w(16) =   0.1058761550973210D+00    
    w(17) =   0.1040733100777293D+00    
    w(18) =   0.1010912737599150D+00    
    w(19) =   0.9696383409440862D-01
    w(20) =   0.9173775713925882D-01
    w(21) =   0.8547225736617248D-01
    w(22) =   0.7823832713576385D-01
    w(23) =   0.7011793325505125D-01
    w(24) =   0.6120309065707916D-01
    w(25) =   0.5159482690249799D-01
    w(26) =   0.4140206251868281D-01
    w(27) =   0.3074049220209360D-01
    w(28) =   0.1973208505612276D-01
    w(29) =   0.8516903878746365D-02

  else if ( n == 30 ) then

    x( 1) =  -0.9968934840746495D+00    
    x( 2) =  -0.9836681232797472D+00    
    x( 3) =  -0.9600218649683075D+00    
    x( 4) =  -0.9262000474292743D+00    
    x( 5) =  -0.8825605357920526D+00    
    x( 6) =  -0.8295657623827684D+00    
    x( 7) =  -0.7677774321048262D+00    
    x( 8) =  -0.6978504947933158D+00    
    x( 9) =  -0.6205261829892429D+00    
    x(10) =  -0.5366241481420199D+00    
    x(11) =  -0.4470337695380892D+00    
    x(12) =  -0.3527047255308781D+00    
    x(13) =  -0.2546369261678899D+00    
    x(14) =  -0.1538699136085835D+00    
    x(15) =  -0.5147184255531770D-01
    x(16) =   0.5147184255531770D-01
    x(17) =   0.1538699136085835D+00    
    x(18) =   0.2546369261678899D+00    
    x(19) =   0.3527047255308781D+00    
    x(20) =   0.4470337695380892D+00    
    x(21) =   0.5366241481420199D+00    
    x(22) =   0.6205261829892429D+00    
    x(23) =   0.6978504947933158D+00    
    x(24) =   0.7677774321048262D+00    
    x(25) =   0.8295657623827684D+00    
    x(26) =   0.8825605357920526D+00    
    x(27) =   0.9262000474292743D+00    
    x(28) =   0.9600218649683075D+00    
    x(29) =   0.9836681232797472D+00    
    x(30) =   0.9968934840746495D+00    
 
    w( 1) =   0.7968192496166648D-02
    w( 2) =   0.1846646831109099D-01
    w( 3) =   0.2878470788332330D-01
    w( 4) =   0.3879919256962704D-01
    w( 5) =   0.4840267283059405D-01
    w( 6) =   0.5749315621761905D-01
    w( 7) =   0.6597422988218052D-01
    w( 8) =   0.7375597473770516D-01
    w( 9) =   0.8075589522942023D-01
    w(10) =   0.8689978720108314D-01
    w(11) =   0.9212252223778619D-01
    w(12) =   0.9636873717464424D-01
    w(13) =   0.9959342058679524D-01
    w(14) =   0.1017623897484056D+00    
    w(15) =   0.1028526528935587D+00    
    w(16) =   0.1028526528935587D+00    
    w(17) =   0.1017623897484056D+00    
    w(18) =   0.9959342058679524D-01
    w(19) =   0.9636873717464424D-01
    w(20) =   0.9212252223778619D-01
    w(21) =   0.8689978720108314D-01
    w(22) =   0.8075589522942023D-01
    w(23) =   0.7375597473770516D-01
    w(24) =   0.6597422988218052D-01
    w(25) =   0.5749315621761905D-01
    w(26) =   0.4840267283059405D-01
    w(27) =   0.3879919256962704D-01
    w(28) =   0.2878470788332330D-01
    w(29) =   0.1846646831109099D-01
    w(30) =   0.7968192496166648D-02

  else if ( n == 31 ) then

    x( 1) =  -0.99708748181947707454263838179654D+00   
    x( 2) =  -0.98468590966515248400211329970113D+00
    x( 3) =  -0.96250392509294966178905249675943D+00
    x( 4) =  -0.93075699789664816495694576311725D+00
    x( 5) =  -0.88976002994827104337419200908023D+00
    x( 6) =  -0.83992032014626734008690453594388D+00
    x( 7) =  -0.78173314841662494040636002019484D+00
    x( 8) =  -0.71577678458685328390597086536649D+00
    x( 9) =  -0.64270672292426034618441820323250D+00
    x(10) =  -0.56324916140714926272094492359516D+00
    x(11) =  -0.47819378204490248044059403935649D+00
    x(12) =  -0.38838590160823294306135146128752D+00
    x(13) =  -0.29471806998170161661790389767170D+00
    x(14) =  -0.19812119933557062877241299603283D+00
    x(15) =  -0.99555312152341520325174790118941D-01
    x(16) =   0.00000000000000000000000000000000D+00  
    x(17) =   0.99555312152341520325174790118941D-01
    x(18) =   0.19812119933557062877241299603283D+00    
    x(19) =   0.29471806998170161661790389767170D+00    
    x(20) =   0.38838590160823294306135146128752D+00    
    x(21) =   0.47819378204490248044059403935649D+00    
    x(22) =   0.56324916140714926272094492359516D+00    
    x(23) =   0.64270672292426034618441820323250D+00    
    x(24) =   0.71577678458685328390597086536649D+00    
    x(25) =   0.78173314841662494040636002019484D+00    
    x(26) =   0.83992032014626734008690453594388D+00    
    x(27) =   0.88976002994827104337419200908023D+00    
    x(28) =   0.93075699789664816495694576311725D+00    
    x(29) =   0.96250392509294966178905249675943D+00    
    x(30) =   0.98468590966515248400211329970113D+00    
    x(31) =   0.99708748181947707454263838179654D+00    
 
    w( 1) =   0.74708315792487746093913218970494D-02
    w( 2) =   0.17318620790310582463552990782414D-01
    w( 3) =   0.27009019184979421800608642617676D-01
    w( 4) =   0.36432273912385464024392008749009D-01
    w( 5) =   0.45493707527201102902315857856518D-01
    w( 6) =   0.54103082424916853711666259085477D-01
    w( 7) =   0.62174786561028426910343543686657D-01
    w( 8) =   0.69628583235410366167756126255124D-01
    w( 9) =   0.76390386598776616426357674901331D-01
    w(10) =   0.82392991761589263903823367431962D-01
    w(11) =   0.87576740608477876126198069695333D-01
    w(12) =   0.91890113893641478215362871607150D-01
    w(13) =   0.95290242912319512807204197487597D-01
    w(14) =   0.97743335386328725093474010978997D-01
    w(15) =   0.99225011226672307874875514428615D-01
    w(16) =   0.99720544793426451427533833734349D-01
    w(17) =   0.99225011226672307874875514428615D-01
    w(18) =   0.97743335386328725093474010978997D-01
    w(19) =   0.95290242912319512807204197487597D-01
    w(20) =   0.91890113893641478215362871607150D-01
    w(21) =   0.87576740608477876126198069695333D-01
    w(22) =   0.82392991761589263903823367431962D-01
    w(23) =   0.76390386598776616426357674901331D-01
    w(24) =   0.69628583235410366167756126255124D-01
    w(25) =   0.62174786561028426910343543686657D-01
    w(26) =   0.54103082424916853711666259085477D-01
    w(27) =   0.45493707527201102902315857856518D-01
    w(28) =   0.36432273912385464024392008749009D-01
    w(29) =   0.27009019184979421800608642617676D-01
    w(30) =   0.17318620790310582463552990782414D-01
    w(31) =   0.74708315792487746093913218970494D-02

  else if ( n == 32 ) then

    x(1) =  - 0.997263861849481563544981128665D+00
    x(2) =  - 0.985611511545268335400175044631D+00
    x(3) =  - 0.964762255587506430773811928118D+00
    x(4) =  - 0.934906075937739689170919134835D+00
    x(5) =  - 0.896321155766052123965307243719D+00
    x(6) =  - 0.849367613732569970133693004968D+00
    x(7) =  - 0.794483795967942406963097298970D+00
    x(8) =  - 0.732182118740289680387426665091D+00
    x(9) =  - 0.663044266930215200975115168663D+00
    x(10) = - 0.587715757240762329040745476402D+00
    x(11) = - 0.506899908932229390023747474378D+00
    x(12) = - 0.421351276130635345364119436172D+00
    x(13) = - 0.331868602282127649779916805730D+00
    x(14) = - 0.239287362252137074544603209166D+00
    x(15) = - 0.144471961582796493485186373599D+00
    x(16) = - 0.483076656877383162348125704405D-01
    x(17) =   0.483076656877383162348125704405D-01
    x(18) =   0.144471961582796493485186373599D+00
    x(19) =   0.239287362252137074544603209166D+00
    x(20) =   0.331868602282127649779916805730D+00
    x(21) =   0.421351276130635345364119436172D+00
    x(22) =   0.506899908932229390023747474378D+00
    x(23) =   0.587715757240762329040745476402D+00
    x(24) =   0.663044266930215200975115168663D+00
    x(25) =   0.732182118740289680387426665091D+00
    x(26) =   0.794483795967942406963097298970D+00
    x(27) =   0.849367613732569970133693004968D+00
    x(28) =   0.896321155766052123965307243719D+00
    x(29) =   0.934906075937739689170919134835D+00
    x(30) =   0.964762255587506430773811928118D+00
    x(31) =   0.985611511545268335400175044631D+00
    x(32) =   0.997263861849481563544981128665D+00

    w(1) =  0.701861000947009660040706373885D-02
    w(2) =  0.162743947309056706051705622064D-01
    w(3) =  0.253920653092620594557525897892D-01
    w(4) =  0.342738629130214331026877322524D-01
    w(5) =  0.428358980222266806568786466061D-01
    w(6) =  0.509980592623761761961632446895D-01
    w(7) =  0.586840934785355471452836373002D-01
    w(8) =  0.658222227763618468376500637069D-01
    w(9) =  0.723457941088485062253993564785D-01
    w(10) = 0.781938957870703064717409188283D-01
    w(11) = 0.833119242269467552221990746043D-01
    w(12) = 0.876520930044038111427714627518D-01
    w(13) = 0.911738786957638847128685771116D-01
    w(14) = 0.938443990808045656391802376681D-01
    w(15) = 0.956387200792748594190820022041D-01
    w(16) = 0.965400885147278005667648300636D-01
    w(17) = 0.965400885147278005667648300636D-01
    w(18) = 0.956387200792748594190820022041D-01
    w(19) = 0.938443990808045656391802376681D-01
    w(20) = 0.911738786957638847128685771116D-01
    w(21) = 0.876520930044038111427714627518D-01
    w(22) = 0.833119242269467552221990746043D-01
    w(23) = 0.781938957870703064717409188283D-01
    w(24) = 0.723457941088485062253993564785D-01
    w(25) = 0.658222227763618468376500637069D-01
    w(26) = 0.586840934785355471452836373002D-01
    w(27) = 0.509980592623761761961632446895D-01
    w(28) = 0.428358980222266806568786466061D-01
    w(29) = 0.342738629130214331026877322524D-01
    w(30) = 0.253920653092620594557525897892D-01
    w(31) = 0.162743947309056706051705622064D-01
    w(32) = 0.701861000947009660040706373885D-02

  else if ( n == 33 ) then

    x( 1) =  -0.9974246942464552D+00    
    x( 2) =  -0.9864557262306425D+00
    x( 3) =  -0.9668229096899927D+00
    x( 4) =  -0.9386943726111684D+00    
    x( 5) =  -0.9023167677434336D+00    
    x( 6) =  -0.8580096526765041D+00    
    x( 7) =  -0.8061623562741665D+00    
    x( 8) =  -0.7472304964495622D+00    
    x( 9) =  -0.6817319599697428D+00    
    x(10) =  -0.6102423458363790D+00    
    x(11) =  -0.5333899047863476D+00    
    x(12) =  -0.4518500172724507D+00    
    x(13) =  -0.3663392577480734D+00    
    x(14) =  -0.2776090971524970D+00    
    x(15) =  -0.1864392988279916D+00    
    x(16) =  -0.09363106585473338D+00
    x(17) =   0.000000000000000D+00
    x(18) =   0.09363106585473338D+00
    x(19) =   0.1864392988279916D+00    
    x(20) =   0.2776090971524970D+00    
    x(21) =   0.3663392577480734D+00    
    x(22) =   0.4518500172724507D+00    
    x(23) =   0.5333899047863476D+00    
    x(24) =   0.6102423458363790D+00    
    x(25) =   0.6817319599697428D+00    
    x(26) =   0.7472304964495622D+00    
    x(27) =   0.8061623562741665D+00    
    x(28) =   0.8580096526765041D+00    
    x(29) =   0.9023167677434336D+00    
    x(30) =   0.9386943726111684D+00    
    x(31) =   0.9668229096899927D+00    
    x(32) =   0.9864557262306425D+00    
    x(33) =   0.9974246942464552D+00    
 
    w( 1) =   0.6606227847587558D-02
    w( 2) =   0.1532170151293465D-01
    w( 3) =   0.2391554810174960D-01
    w( 4) =   0.3230035863232891D-01
    w( 5) =   0.4040154133166965D-01
    w( 6) =   0.4814774281871162D-01
    w( 7) =   0.5547084663166357D-01
    w( 8) =   0.6230648253031755D-01
    w( 9) =   0.6859457281865676D-01
    w(10) =   0.7427985484395420D-01
    w(11) =   0.7931236479488685D-01
    w(12) =   0.8364787606703869D-01
    w(13) =   0.8724828761884425D-01
    w(14) =   0.9008195866063859D-01
    w(15) =   0.9212398664331678D-01
    w(16) =   0.9335642606559612D-01
    w(17) =   0.9376844616020999D-01
    w(18) =   0.9335642606559612D-01
    w(19) =   0.9212398664331678D-01
    w(20) =   0.9008195866063859D-01
    w(21) =   0.8724828761884425D-01
    w(22) =   0.8364787606703869D-01
    w(23) =   0.7931236479488685D-01
    w(24) =   0.7427985484395420D-01
    w(25) =   0.6859457281865676D-01
    w(26) =   0.6230648253031755D-01
    w(27) =   0.5547084663166357D-01
    w(28) =   0.4814774281871162D-01
    w(29) =   0.4040154133166965D-01
    w(30) =   0.3230035863232891D-01
    w(31) =   0.2391554810174960D-01
    w(32) =   0.1532170151293465D-01
    w(33) =   0.6606227847587558D-02

  else if ( n == 63 ) then

    x( 1) =  -0.99928298402912378050701628988630D+00    
    x( 2) =  -0.99622401277797010860209018267357D+00    
    x( 3) =  -0.99072854689218946681089469460884D+00    
    x( 4) =  -0.98280881059372723486251140727639D+00    
    x( 5) =  -0.97248403469757002280196067864927D+00    
    x( 6) =  -0.95977944975894192707035416626398D+00    
    x( 7) =  -0.94472613404100980296637531962798D+00    
    x( 8) =  -0.92736092062184320544703138132518D+00    
    x( 9) =  -0.90772630277853155803695313291596D+00    
    x(10) =  -0.88587032850785342629029845731337D+00    
    x(11) =  -0.86184648236412371953961183943106D+00    
    x(12) =  -0.83571355431950284347180776961571D+00    
    x(13) =  -0.80753549577345676005146598636324D+00    
    x(14) =  -0.77738126299037233556333018991104D+00    
    x(15) =  -0.74532464831784741782932166103759D+00    
    x(16) =  -0.71144409958484580785143153770401D+00    
    x(17) =  -0.67582252811498609013110331596954D+00    
    x(18) =  -0.63854710582136538500030695387338D+00    
    x(19) =  -0.59970905187762523573900892686880D+00    
    x(20) =  -0.55940340948628501326769780007005D+00    
    x(21) =  -0.51772881329003324812447758452632D+00    
    x(22) =  -0.47478724799480439992221230985149D+00    
    x(23) =  -0.43068379879511160066208893391863D+00    
    x(24) =  -0.38552639421224789247761502227440D+00    
    x(25) =  -0.33942554197458440246883443159432D+00    
    x(26) =  -0.29249405858625144003615715555067D+00    
    x(27) =  -0.24484679324595336274840459392483D+00    
    x(28) =  -0.19660034679150668455762745706572D+00    
    x(29) =  -0.14787278635787196856983909655297D+00    
    x(30) =  -0.98783356446945279529703669453922D-01
    x(31) =  -0.49452187116159627234233818051808D-01
    x(32) =    0.0000000000000000000000000000000D+00    
    x(33) =   0.49452187116159627234233818051808D-01
    x(34) =   0.98783356446945279529703669453922D-01
    x(35) =   0.14787278635787196856983909655297D+00    
    x(36) =   0.19660034679150668455762745706572D+00    
    x(37) =   0.24484679324595336274840459392483D+00    
    x(38) =   0.29249405858625144003615715555067D+00    
    x(39) =   0.33942554197458440246883443159432D+00    
    x(40) =   0.38552639421224789247761502227440D+00    
    x(41) =   0.43068379879511160066208893391863D+00    
    x(42) =   0.47478724799480439992221230985149D+00    
    x(43) =   0.51772881329003324812447758452632D+00    
    x(44) =   0.55940340948628501326769780007005D+00    
    x(45) =   0.59970905187762523573900892686880D+00    
    x(46) =   0.63854710582136538500030695387338D+00    
    x(47) =   0.67582252811498609013110331596954D+00    
    x(48) =   0.71144409958484580785143153770401D+00    
    x(49) =   0.74532464831784741782932166103759D+00    
    x(50) =   0.77738126299037233556333018991104D+00    
    x(51) =   0.80753549577345676005146598636324D+00    
    x(52) =   0.83571355431950284347180776961571D+00    
    x(53) =   0.86184648236412371953961183943106D+00    
    x(54) =   0.88587032850785342629029845731337D+00    
    x(55) =   0.90772630277853155803695313291596D+00    
    x(56) =   0.92736092062184320544703138132518D+00    
    x(57) =   0.94472613404100980296637531962798D+00    
    x(58) =   0.95977944975894192707035416626398D+00    
    x(59) =   0.97248403469757002280196067864927D+00    
    x(60) =   0.98280881059372723486251140727639D+00    
    x(61) =   0.99072854689218946681089469460884D+00    
    x(62) =   0.99622401277797010860209018267357D+00    
    x(63) =   0.99928298402912378050701628988630D+00

    w( 1) =   0.18398745955770837880499331680577D-02
    w( 2) =   0.42785083468637618661951422543371D-02
    w( 3) =   0.67102917659601362519069109850892D-02
    w( 4) =   0.91259686763266563540586445877022D-02
    w( 5) =   0.11519376076880041750750606118707D-01
    w( 6) =   0.13884612616115610824866086365937D-01
    w( 7) =   0.16215878410338338882283672974995D-01
    w( 8) =   0.18507464460161270409260545805144D-01
    w( 9) =   0.20753761258039090775341953421471D-01
    w(10) =   0.22949271004889933148942319561770D-01
    w(11) =   0.25088620553344986618630138068443D-01
    w(12) =   0.27166574359097933225189839439413D-01
    w(13) =   0.29178047208280526945551502154029D-01
    w(14) =   0.31118116622219817508215988557189D-01
    w(15) =   0.32982034883779341765683179672459D-01
    w(16) =   0.34765240645355877697180504642788D-01
    w(17) =   0.36463370085457289630452409787542D-01
    w(18) =   0.38072267584349556763638324927889D-01
    w(19) =   0.39587995891544093984807928149202D-01
    w(20) =   0.41006845759666398635110037009072D-01
    w(21) =   0.42325345020815822982505485403028D-01
    w(22) =   0.43540267083027590798964315704401D-01
    w(23) =   0.44648638825941395370332669516813D-01
    w(24) =   0.45647747876292608685885992608542D-01
    w(25) =   0.46535149245383696510395418746953D-01
    w(26) =   0.47308671312268919080604988338844D-01
    w(27) =   0.47966421137995131411052756195132D-01
    w(28) =   0.48506789097883847864090099145802D-01
    w(29) =   0.48928452820511989944709361549215D-01
    w(30) =   0.49230380423747560785043116988145D-01
    w(31) =   0.49411833039918178967039646116705D-01
    w(32) =   0.49472366623931020888669360420926D-01
    w(33) =   0.49411833039918178967039646116705D-01
    w(34) =   0.49230380423747560785043116988145D-01
    w(35) =   0.48928452820511989944709361549215D-01
    w(36) =   0.48506789097883847864090099145802D-01
    w(37) =   0.47966421137995131411052756195132D-01
    w(38) =   0.47308671312268919080604988338844D-01
    w(39) =   0.46535149245383696510395418746953D-01
    w(40) =   0.45647747876292608685885992608542D-01
    w(41) =   0.44648638825941395370332669516813D-01
    w(42) =   0.43540267083027590798964315704401D-01
    w(43) =   0.42325345020815822982505485403028D-01
    w(44) =   0.41006845759666398635110037009072D-01
    w(45) =   0.39587995891544093984807928149202D-01
    w(46) =   0.38072267584349556763638324927889D-01
    w(47) =   0.36463370085457289630452409787542D-01
    w(48) =   0.34765240645355877697180504642788D-01
    w(49) =   0.32982034883779341765683179672459D-01
    w(50) =   0.31118116622219817508215988557189D-01
    w(51) =   0.29178047208280526945551502154029D-01
    w(52) =   0.27166574359097933225189839439413D-01
    w(53) =   0.25088620553344986618630138068443D-01
    w(54) =   0.22949271004889933148942319561770D-01
    w(55) =   0.20753761258039090775341953421471D-01
    w(56) =   0.18507464460161270409260545805144D-01
    w(57) =   0.16215878410338338882283672974995D-01
    w(58) =   0.13884612616115610824866086365937D-01
    w(59) =   0.11519376076880041750750606118707D-01
    w(60) =   0.91259686763266563540586445877022D-02
    w(61) =   0.67102917659601362519069109850892D-02
    w(62) =   0.42785083468637618661951422543371D-02
    w(63) =   0.18398745955770837880499331680577D-02
 
  else if ( n == 64 ) then

    x(1) =  - 0.999305041735772139456905624346D+00
    x(2) =  - 0.996340116771955279346924500676D+00
    x(3) =  - 0.991013371476744320739382383443D+00
    x(4) =  - 0.983336253884625956931299302157D+00
    x(5) =  - 0.973326827789910963741853507352D+00
    x(6) =  - 0.961008799652053718918614121897D+00
    x(7) =  - 0.946411374858402816062481491347D+00
    x(8) =  - 0.929569172131939575821490154559D+00
    x(9) =  - 0.910522137078502805756380668008D+00
    x(10) = - 0.889315445995114105853404038273D+00
    x(11) = - 0.865999398154092819760783385070D+00
    x(12) = - 0.840629296252580362751691544696D+00
    x(13) = - 0.813265315122797559741923338086D+00
    x(14) = - 0.783972358943341407610220525214D+00
    x(15) = - 0.752819907260531896611863774886D+00
    x(16) = - 0.719881850171610826848940217832D+00
    x(17) = - 0.685236313054233242563558371031D+00
    x(18) = - 0.648965471254657339857761231993D+00
    x(19) = - 0.611155355172393250248852971019D+00
    x(20) = - 0.571895646202634034283878116659D+00
    x(21) = - 0.531279464019894545658013903544D+00
    x(22) = - 0.489403145707052957478526307022D+00
    x(23) = - 0.446366017253464087984947714759D+00
    x(24) = - 0.402270157963991603695766771260D+00
    x(25) = - 0.357220158337668115950442615046D+00
    x(26) = - 0.311322871990210956157512698560D+00
    x(27) = - 0.264687162208767416373964172510D+00
    x(28) = - 0.217423643740007084149648748989D+00
    x(29) = - 0.169644420423992818037313629748D+00
    x(30) = - 0.121462819296120554470376463492D+00
    x(31) = - 0.729931217877990394495429419403D-01
    x(32) = - 0.243502926634244325089558428537D-01
    x(33) =   0.243502926634244325089558428537D-01
    x(34) =   0.729931217877990394495429419403D-01
    x(35) =   0.121462819296120554470376463492D+00
    x(36) =   0.169644420423992818037313629748D+00
    x(37) =   0.217423643740007084149648748989D+00
    x(38) =   0.264687162208767416373964172510D+00
    x(39) =   0.311322871990210956157512698560D+00
    x(40) =   0.357220158337668115950442615046D+00
    x(41) =   0.402270157963991603695766771260D+00
    x(42) =   0.446366017253464087984947714759D+00
    x(43) =   0.489403145707052957478526307022D+00
    x(44) =   0.531279464019894545658013903544D+00
    x(45) =   0.571895646202634034283878116659D+00
    x(46) =   0.611155355172393250248852971019D+00
    x(47) =   0.648965471254657339857761231993D+00
    x(48) =   0.685236313054233242563558371031D+00
    x(49) =   0.719881850171610826848940217832D+00
    x(50) =   0.752819907260531896611863774886D+00
    x(51) =   0.783972358943341407610220525214D+00
    x(52) =   0.813265315122797559741923338086D+00
    x(53) =   0.840629296252580362751691544696D+00
    x(54) =   0.865999398154092819760783385070D+00
    x(55) =   0.889315445995114105853404038273D+00
    x(56) =   0.910522137078502805756380668008D+00
    x(57) =   0.929569172131939575821490154559D+00
    x(58) =   0.946411374858402816062481491347D+00
    x(59) =   0.961008799652053718918614121897D+00
    x(60) =   0.973326827789910963741853507352D+00
    x(61) =   0.983336253884625956931299302157D+00
    x(62) =   0.991013371476744320739382383443D+00
    x(63) =   0.996340116771955279346924500676D+00
    x(64) =   0.999305041735772139456905624346D+00

    w(1) =  0.178328072169643294729607914497D-02
    w(2) =  0.414703326056246763528753572855D-02
    w(3) =  0.650445796897836285611736039998D-02
    w(4) =  0.884675982636394772303091465973D-02
    w(5) =  0.111681394601311288185904930192D-01
    w(6) =  0.134630478967186425980607666860D-01
    w(7) =  0.157260304760247193219659952975D-01
    w(8) =  0.179517157756973430850453020011D-01
    w(9) =  0.201348231535302093723403167285D-01
    w(10) = 0.222701738083832541592983303842D-01
    w(11) = 0.243527025687108733381775504091D-01
    w(12) = 0.263774697150546586716917926252D-01
    w(13) = 0.283396726142594832275113052002D-01
    w(14) = 0.302346570724024788679740598195D-01
    w(15) = 0.320579283548515535854675043479D-01
    w(16) = 0.338051618371416093915654821107D-01
    w(17) = 0.354722132568823838106931467152D-01
    w(18) = 0.370551285402400460404151018096D-01
    w(19) = 0.385501531786156291289624969468D-01
    w(20) = 0.399537411327203413866569261283D-01
    w(21) = 0.412625632426235286101562974736D-01
    w(22) = 0.424735151236535890073397679088D-01
    w(23) = 0.435837245293234533768278609737D-01
    w(24) = 0.445905581637565630601347100309D-01
    w(25) = 0.454916279274181444797709969713D-01
    w(26) = 0.462847965813144172959532492323D-01
    w(27) = 0.469681828162100173253262857546D-01
    w(28) = 0.475401657148303086622822069442D-01
    w(29) = 0.479993885964583077281261798713D-01
    w(30) = 0.483447622348029571697695271580D-01
    w(31) = 0.485754674415034269347990667840D-01
    w(32) = 0.486909570091397203833653907347D-01
    w(33) = 0.486909570091397203833653907347D-01
    w(34) = 0.485754674415034269347990667840D-01
    w(35) = 0.483447622348029571697695271580D-01
    w(36) = 0.479993885964583077281261798713D-01
    w(37) = 0.475401657148303086622822069442D-01
    w(38) = 0.469681828162100173253262857546D-01
    w(39) = 0.462847965813144172959532492323D-01
    w(40) = 0.454916279274181444797709969713D-01
    w(41) = 0.445905581637565630601347100309D-01
    w(42) = 0.435837245293234533768278609737D-01
    w(43) = 0.424735151236535890073397679088D-01
    w(44) = 0.412625632426235286101562974736D-01
    w(45) = 0.399537411327203413866569261283D-01
    w(46) = 0.385501531786156291289624969468D-01
    w(47) = 0.370551285402400460404151018096D-01
    w(48) = 0.354722132568823838106931467152D-01
    w(49) = 0.338051618371416093915654821107D-01
    w(50) = 0.320579283548515535854675043479D-01
    w(51) = 0.302346570724024788679740598195D-01
    w(52) = 0.283396726142594832275113052002D-01
    w(53) = 0.263774697150546586716917926252D-01
    w(54) = 0.243527025687108733381775504091D-01
    w(55) = 0.222701738083832541592983303842D-01
    w(56) = 0.201348231535302093723403167285D-01
    w(57) = 0.179517157756973430850453020011D-01
    w(58) = 0.157260304760247193219659952975D-01
    w(59) = 0.134630478967186425980607666860D-01
    w(60) = 0.111681394601311288185904930192D-01
    w(61) = 0.884675982636394772303091465973D-02
    w(62) = 0.650445796897836285611736039998D-02
    w(63) = 0.414703326056246763528753572855D-02
    w(64) = 0.178328072169643294729607914497D-02

  else if ( n == 65 ) then

    x( 1) =  -0.9993260970754129D+00    
    x( 2) =  -0.9964509480618492D+00    
    x( 3) =  -0.9912852761768016D+00    
    x( 4) =  -0.9838398121870350D+00    
    x( 5) =  -0.9741315398335512D+00    
    x( 6) =  -0.9621827547180553D+00    
    x( 7) =  -0.9480209281684076D+00    
    x( 8) =  -0.9316786282287494D+00    
    x( 9) =  -0.9131934405428462D+00    
    x(10) =  -0.8926078805047389D+00    
    x(11) =  -0.8699692949264071D+00    
    x(12) =  -0.8453297528999303D+00    
    x(13) =  -0.8187459259226514D+00    
    x(14) =  -0.7902789574921218D+00    
    x(15) =  -0.7599943224419998D+00    
    x(16) =  -0.7279616763294247D+00    
    x(17) =  -0.6942546952139916D+00    
    x(18) =  -0.6589509061936252D+00    
    x(19) =  -0.6221315090854003D+00    
    x(20) =  -0.5838811896604873D+00    
    x(21) =  -0.5442879248622271D+00    
    x(22) =  -0.5034427804550069D+00    
    x(23) =  -0.4614397015691450D+00    
    x(24) =  -0.4183752966234090D+00    
    x(25) =  -0.3743486151220660D+00    
    x(26) =  -0.3294609198374864D+00    
    x(27) =  -0.2838154539022487D+00    
    x(28) =  -0.2375172033464168D+00    
    x(29) =  -0.1906726556261428D+00    
    x(30) =  -0.1433895546989752D+00    
    x(31) =  -0.9577665320919751D-01
    x(32) =  -0.4794346235317186D-01
    x(33) =    0.000000000000000D+00    
    x(34) =   0.4794346235317186D-01
    x(35) =   0.9577665320919751D-01
    x(36) =   0.1433895546989752D+00    
    x(37) =   0.1906726556261428D+00    
    x(38) =   0.2375172033464168D+00    
    x(39) =   0.2838154539022487D+00    
    x(40) =   0.3294609198374864D+00    
    x(41) =   0.3743486151220660D+00    
    x(42) =   0.4183752966234090D+00    
    x(43) =   0.4614397015691450D+00    
    x(44) =   0.5034427804550069D+00    
    x(45) =   0.5442879248622271D+00    
    x(46) =   0.5838811896604873D+00    
    x(47) =   0.6221315090854003D+00    
    x(48) =   0.6589509061936252D+00    
    x(49) =   0.6942546952139916D+00    
    x(50) =   0.7279616763294247D+00    
    x(51) =   0.7599943224419998D+00    
    x(52) =   0.7902789574921218D+00    
    x(53) =   0.8187459259226514D+00    
    x(54) =   0.8453297528999303D+00    
    x(55) =   0.8699692949264071D+00    
    x(56) =   0.8926078805047389D+00    
    x(57) =   0.9131934405428462D+00    
    x(58) =   0.9316786282287494D+00    
    x(59) =   0.9480209281684076D+00    
    x(60) =   0.9621827547180553D+00    
    x(61) =   0.9741315398335512D+00    
    x(62) =   0.9838398121870350D+00    
    x(63) =   0.9912852761768016D+00    
    x(64) =   0.9964509480618492D+00    
    x(65) =   0.9993260970754129D+00    
 
    w( 1) =   0.1729258251300218D-02
    w( 2) =   0.4021524172003703D-02
    w( 3) =   0.6307942578971821D-02
    w( 4) =   0.8580148266881443D-02
    w( 5) =   0.1083267878959798D-01
    w( 6) =   0.1306031163999490D-01
    w( 7) =   0.1525791214644825D-01
    w( 8) =   0.1742042199767025D-01
    w( 9) =   0.1954286583675005D-01
    w(10) =   0.2162036128493408D-01
    w(11) =   0.2364812969128723D-01
    w(12) =   0.2562150693803776D-01
    w(13) =   0.2753595408845034D-01
    w(14) =   0.2938706778931066D-01
    w(15) =   0.3117059038018911D-01
    w(16) =   0.3288241967636860D-01
    w(17) =   0.3451861839854901D-01
    w(18) =   0.3607542322556527D-01
    w(19) =   0.3754925344825770D-01
    w(20) =   0.3893671920405121D-01
    w(21) =   0.4023462927300549D-01
    w(22) =   0.4143999841724028D-01
    w(23) =   0.4255005424675579D-01
    w(24) =   0.4356224359580051D-01
    w(25) =   0.4447423839508296D-01
    w(26) =   0.4528394102630023D-01
    w(27) =   0.4598948914665173D-01
    w(28) =   0.4658925997223349D-01
    w(29) =   0.4708187401045461D-01
    w(30) =   0.4746619823288551D-01
    w(31) =   0.4774134868124067D-01
    w(32) =   0.4790669250049590D-01
    w(33) =   0.4796184939446662D-01
    w(34) =   0.4790669250049590D-01
    w(35) =   0.4774134868124067D-01
    w(36) =   0.4746619823288551D-01
    w(37) =   0.4708187401045461D-01
    w(38) =   0.4658925997223349D-01
    w(39) =   0.4598948914665173D-01
    w(40) =   0.4528394102630023D-01
    w(41) =   0.4447423839508296D-01
    w(42) =   0.4356224359580051D-01
    w(43) =   0.4255005424675579D-01
    w(44) =   0.4143999841724028D-01
    w(45) =   0.4023462927300549D-01
    w(46) =   0.3893671920405121D-01
    w(47) =   0.3754925344825770D-01
    w(48) =   0.3607542322556527D-01
    w(49) =   0.3451861839854901D-01
    w(50) =   0.3288241967636860D-01
    w(51) =   0.3117059038018911D-01
    w(52) =   0.2938706778931066D-01
    w(53) =   0.2753595408845034D-01
    w(54) =   0.2562150693803776D-01
    w(55) =   0.2364812969128723D-01
    w(56) =   0.2162036128493408D-01
    w(57) =   0.1954286583675005D-01
    w(58) =   0.1742042199767025D-01
    w(59) =   0.1525791214644825D-01
    w(60) =   0.1306031163999490D-01
    w(61) =   0.1083267878959798D-01
    w(62) =   0.8580148266881443D-02
    w(63) =   0.6307942578971821D-02
    w(64) =   0.4021524172003703D-02
    w(65) =   0.1729258251300218D-02
    
  else if ( n == 127 ) then
  
    x(  1) =  -0.99982213041530614629963254927125D+00    
    x(  2) =  -0.99906293435531189513828920479421D+00    
    x(  3) =  -0.99769756618980462107441703193392D+00    
    x(  4) =  -0.99572655135202722663543337085008D+00    
    x(  5) =  -0.99315104925451714736113079489080D+00    
    x(  6) =  -0.98997261459148415760778669967548D+00    
    x(  7) =  -0.98619317401693166671043833175407D+00    
    x(  8) =  -0.98181502080381411003346312451200D+00    
    x(  9) =  -0.97684081234307032681744391886221D+00    
    x( 10) =  -0.97127356816152919228894689830512D+00    
    x( 11) =  -0.96511666794529212109082507703391D+00    
    x( 12) =  -0.95837384942523877114910286998060D+00    
    x( 13) =  -0.95104920607788031054790764659636D+00    
    x( 14) =  -0.94314718462481482734544963026201D+00    
    x( 15) =  -0.93467258232473796857363487794906D+00    
    x( 16) =  -0.92563054405623384912746466814259D+00    
    x( 17) =  -0.91602655919146580931308861741716D+00    
    x( 18) =  -0.90586645826182138280246131760282D+00    
    x( 19) =  -0.89515640941708370896904382642451D+00    
    x( 20) =  -0.88390291468002656994525794802849D+00    
    x( 21) =  -0.87211280599856071141963753428864D+00    
    x( 22) =  -0.85979324109774080981203134414483D+00    
    x( 23) =  -0.84695169913409759845333931085437D+00    
    x( 24) =  -0.83359597615489951437955716480123D+00    
    x( 25) =  -0.81973418036507867415511910167470D+00    
    x( 26) =  -0.80537472720468021466656079404644D+00    
    x( 27) =  -0.79052633423981379994544995252740D+00    
    x( 28) =  -0.77519801587020238244496276354566D+00    
    x( 29) =  -0.75939907785653667155666366659810D+00    
    x( 30) =  -0.74313911167095451292056688997595D+00   
    x( 31) =  -0.72642798867407268553569290153270D+00    
    x( 32) =  -0.70927585412210456099944463906757D+00    
    x( 33) =  -0.69169312100770067015644143286666D+00    
    x( 34) =  -0.67369046373825048534668253831602D+00    
    x( 35) =  -0.65527881165548263027676505156852D+00    
    x( 36) =  -0.63646934240029724134760815684175D+00    
    x( 37) =  -0.61727347512685828385763916340822D+00    
    x( 38) =  -0.59770286357006522938441201887478D+00    
    x( 39) =  -0.57776938897061258000325165713764D+00    
    x( 40) =  -0.55748515286193223292186190687872D+00    
    x( 41) =  -0.53686246972339756745816636353452D+00    
    x( 42) =  -0.51591385950424935727727729906662D+00    
    x( 43) =  -0.49465204002278211739494017368636D+00    
    x( 44) =  -0.47308991924540524164509989939699D+00    
    x( 45) =  -0.45124058745026622733189858020729D+00    
    x( 46) =  -0.42911730928019337626254405355418D+00    
    x( 47) =  -0.40673351568978256340867288124339D+00    
    x( 48) =  -0.38410279579151693577907781452239D+00    
    x( 49) =  -0.36123888860586970607092484346723D+00    
    x( 50) =  -0.33815567472039850137600027657095D+00    
    x( 51) =  -0.31486716786289498148601475374890D+00    
    x( 52) =  -0.29138750639370562079451875284568D+00    
    x( 53) =  -0.26773094472238862088834352027938D+00    
    x( 54) =  -0.24391184465391785797071324453138D+00    
    x( 55) =  -0.21994466666968754245452337866940D+00    
    x( 56) =  -0.19584396114861085150428162519610D+00    
    x( 57) =  -0.17162435953364216500834492248954D+00    
    x( 58) =  -0.14730056544908566938932929319807D+00    
    x( 59) =  -0.12288734577408297172603365288567D+00    
    x( 60) =  -0.98399521677698970751091751509101D-01
    x( 61) =  -0.73851959621048545273440409360569D-01
    x( 62) =  -0.49259562331926630315379321821927D-01
    x( 63) =  -0.24637259757420944614897071846088D-01
    x( 64) =   0.00000000000000000000000000000000D+00    
    x( 65) =   0.24637259757420944614897071846088D-01
    x( 66) =   0.49259562331926630315379321821927D-01
    x( 67) =   0.73851959621048545273440409360569D-01
    x( 68) =   0.98399521677698970751091751509101D-01
    x( 69) =   0.12288734577408297172603365288567D+00    
    x( 70) =   0.14730056544908566938932929319807D+00    
    x( 71) =   0.17162435953364216500834492248954D+00    
    x( 72) =   0.19584396114861085150428162519610D+00    
    x( 73) =   0.21994466666968754245452337866940D+00    
    x( 74) =   0.24391184465391785797071324453138D+00    
    x( 75) =   0.26773094472238862088834352027938D+00    
    x( 76) =   0.29138750639370562079451875284568D+00    
    x( 77) =   0.31486716786289498148601475374890D+00    
    x( 78) =   0.33815567472039850137600027657095D+00    
    x( 79) =   0.36123888860586970607092484346723D+00    
    x( 80) =   0.38410279579151693577907781452239D+00    
    x( 81) =   0.40673351568978256340867288124339D+00    
    x( 82) =   0.42911730928019337626254405355418D+00    
    x( 83) =   0.45124058745026622733189858020729D+00    
    x( 84) =   0.47308991924540524164509989939699D+00    
    x( 85) =   0.49465204002278211739494017368636D+00    
    x( 86) =   0.51591385950424935727727729906662D+00    
    x( 87) =   0.53686246972339756745816636353452D+00    
    x( 88) =   0.55748515286193223292186190687872D+00    
    x( 89) =   0.57776938897061258000325165713764D+00    
    x( 90) =   0.59770286357006522938441201887478D+00    
    x( 91) =   0.61727347512685828385763916340822D+00    
    x( 92) =   0.63646934240029724134760815684175D+00    
    x( 93) =   0.65527881165548263027676505156852D+00    
    x( 94) =   0.67369046373825048534668253831602D+00    
    x( 95) =   0.69169312100770067015644143286666D+00   
    x( 96) =   0.70927585412210456099944463906757D+00    
    x( 97) =   0.72642798867407268553569290153270D+00    
    x( 98) =   0.74313911167095451292056688997595D+00    
    x( 99) =   0.75939907785653667155666366659810D+00    
    x(100) =   0.77519801587020238244496276354566D+00    
    x(101) =   0.79052633423981379994544995252740D+00    
    x(102) =   0.80537472720468021466656079404644D+00    
    x(103) =   0.81973418036507867415511910167470D+00    
    x(104) =   0.83359597615489951437955716480123D+00    
    x(105) =   0.84695169913409759845333931085437D+00    
    x(106) =   0.85979324109774080981203134414483D+00    
    x(107) =   0.87211280599856071141963753428864D+00    
    x(108) =   0.88390291468002656994525794802849D+00    
    x(109) =   0.89515640941708370896904382642451D+00    
    x(110) =   0.90586645826182138280246131760282D+00    
    x(111) =   0.91602655919146580931308861741716D+00    
    x(112) =   0.92563054405623384912746466814259D+00    
    x(113) =   0.93467258232473796857363487794906D+00    
    x(114) =   0.94314718462481482734544963026201D+00    
    x(115) =   0.95104920607788031054790764659636D+00    
    x(116) =   0.95837384942523877114910286998060D+00    
    x(117) =   0.96511666794529212109082507703391D+00    
    x(118) =   0.97127356816152919228894689830512D+00    
    x(119) =   0.97684081234307032681744391886221D+00    
    x(120) =   0.98181502080381411003346312451200D+00    
    x(121) =   0.98619317401693166671043833175407D+00    
    x(122) =   0.98997261459148415760778669967548D+00    
    x(123) =   0.99315104925451714736113079489080D+00    
    x(124) =   0.99572655135202722663543337085008D+00    
    x(125) =   0.99769756618980462107441703193392D+00    
    x(126) =   0.99906293435531189513828920479421D+00    
    x(127) =   0.99982213041530614629963254927125D+00    

    w(  1) =   0.45645726109586654495731936146574D-03
    w(  2) =   0.10622766869538486959954760554099D-02
    w(  3) =   0.16683488125171936761028811985672D-02
    w(  4) =   0.22734860707492547802810838362671D-02
    w(  5) =   0.28772587656289004082883197417581D-02
    w(  6) =   0.34792893810051465908910894094105D-02
    w(  7) =   0.40792095178254605327114733456293D-02
    w(  8) =   0.46766539777779034772638165662478D-02
    w(  9) =   0.52712596565634400891303815906251D-02
    w( 10) =   0.58626653903523901033648343751367D-02
    w( 11) =   0.64505120486899171845442463868748D-02
    w( 12) =   0.70344427036681608755685893032552D-02
    w( 13) =   0.76141028256526859356393930849227D-02
    w( 14) =   0.81891404887415730817235884718726D-02
    w( 15) =   0.87592065795403145773316804234385D-02
    w( 16) =   0.93239550065309714787536985834029D-02
    w( 17) =   0.98830429087554914716648010899606D-02
    w( 18) =   0.10436130863141005225673171997668D-01
    w( 19) =   0.10982883090068975788799657376065D-01
    w( 20) =   0.11522967656921087154811609734510D-01
    w( 21) =   0.12056056679400848183529562144697D-01
    w( 22) =   0.12581826520465013101514365424172D-01
    w( 23) =   0.13099957986718627426172681912499D-01
    w( 24) =   0.13610136522139249906034237533759D-01
    w( 25) =   0.14112052399003395774044161633613D-01
    w( 26) =   0.14605400905893418351737288078952D-01
    w( 27) =   0.15089882532666922992635733981431D-01
    w( 28) =   0.15565203152273955098532590262975D-01
    w( 29) =   0.16031074199309941802254151842763D-01
    w( 30) =   0.16487212845194879399346060358146D-01
    w( 31) =   0.16933342169871654545878815295200D-01
    w( 32) =   0.17369191329918731922164721250350D-01
    w( 33) =   0.17794495722974774231027912900351D-01
    w( 34) =   0.18208997148375106468721469154479D-01
    w( 35) =   0.18612443963902310429440419898958D-01
    w( 36) =   0.19004591238555646611148901044533D-01
    w( 37) =   0.19385200901246454628112623489471D-01
    w( 38) =   0.19754041885329183081815217323169D-01
    w( 39) =   0.20110890268880247225644623956287D-01
    w( 40) =   0.20455529410639508279497065713301D-01
    w( 41) =   0.20787750081531811812652137291250D-01
    w( 42) =   0.21107350591688713643523847921658D-01
    w( 43) =   0.21414136912893259295449693233545D-01
    w( 44) =   0.21707922796373466052301324695331D-01
    w( 45) =   0.21988529885872983756478409758807D-01
    w( 46) =   0.22255787825930280235631416460158D-01
    w( 47) =   0.22509534365300608085694429903050D-01
    w( 48) =   0.22749615455457959852242553240982D-01
    w( 49) =   0.22975885344117206754377437838947D-01
    w( 50) =   0.23188206663719640249922582981729D-01
    w( 51) =   0.23386450514828194170722043496950D-01
    w( 52) =   0.23570496544381716050033676844306D-01
    w( 53) =   0.23740233018760777777714726703424D-01
    w( 54) =   0.23895556891620665983864481754172D-01
    w( 55) =   0.24036373866450369675132086026456D-01
    w( 56) =   0.24162598453819584716522917710986D-01
    w( 57) =   0.24274154023278979833195063936748D-01
    w( 58) =   0.24370972849882214952813561907241D-01
    w( 59) =   0.24452996155301467956140198471529D-01
    w( 60) =   0.24520174143511508275183033290175D-01
    w( 61) =   0.24572466031020653286354137335186D-01
    w( 62) =   0.24609840071630254092545634003360D-01
    w( 63) =   0.24632273575707679066033370218017D-01
    w( 64) =   0.24639752923961094419579417477503D-01
    w( 65) =   0.24632273575707679066033370218017D-01
    w( 66) =   0.24609840071630254092545634003360D-01
    w( 67) =   0.24572466031020653286354137335186D-01
    w( 68) =   0.24520174143511508275183033290175D-01
    w( 69) =   0.24452996155301467956140198471529D-01
    w( 70) =   0.24370972849882214952813561907241D-01
    w( 71) =   0.24274154023278979833195063936748D-01
    w( 72) =   0.24162598453819584716522917710986D-01
    w( 73) =   0.24036373866450369675132086026456D-01
    w( 74) =   0.23895556891620665983864481754172D-01
    w( 75) =   0.23740233018760777777714726703424D-01
    w( 76) =   0.23570496544381716050033676844306D-01
    w( 77) =   0.23386450514828194170722043496950D-01
    w( 78) =   0.23188206663719640249922582981729D-01
    w( 79) =   0.22975885344117206754377437838947D-01
    w( 80) =   0.22749615455457959852242553240982D-01
    w( 81) =   0.22509534365300608085694429903050D-01
    w( 82) =   0.22255787825930280235631416460158D-01
    w( 83) =   0.21988529885872983756478409758807D-01
    w( 84) =   0.21707922796373466052301324695331D-01
    w( 85) =   0.21414136912893259295449693233545D-01
    w( 86) =   0.21107350591688713643523847921658D-01
    w( 87) =   0.20787750081531811812652137291250D-01
    w( 88) =   0.20455529410639508279497065713301D-01
    w( 89) =   0.20110890268880247225644623956287D-01
    w( 90) =   0.19754041885329183081815217323169D-01
    w( 91) =   0.19385200901246454628112623489471D-01
    w( 92) =   0.19004591238555646611148901044533D-01
    w( 93) =   0.18612443963902310429440419898958D-01
    w( 94) =   0.18208997148375106468721469154479D-01
    w( 95) =   0.17794495722974774231027912900351D-01
    w( 96) =   0.17369191329918731922164721250350D-01
    w( 97) =   0.16933342169871654545878815295200D-01
    w( 98) =   0.16487212845194879399346060358146D-01
    w( 99) =   0.16031074199309941802254151842763D-01
    w(100) =   0.15565203152273955098532590262975D-01
    w(101) =   0.15089882532666922992635733981431D-01
    w(102) =   0.14605400905893418351737288078952D-01
    w(103) =   0.14112052399003395774044161633613D-01
    w(104) =   0.13610136522139249906034237533759D-01
    w(105) =   0.13099957986718627426172681912499D-01
    w(106) =   0.12581826520465013101514365424172D-01
    w(107) =   0.12056056679400848183529562144697D-01
    w(108) =   0.11522967656921087154811609734510D-01
    w(109) =   0.10982883090068975788799657376065D-01
    w(110) =   0.10436130863141005225673171997668D-01
    w(111) =   0.98830429087554914716648010899606D-02
    w(112) =   0.93239550065309714787536985834029D-02
    w(113) =   0.87592065795403145773316804234385D-02
    w(114) =   0.81891404887415730817235884718726D-02
    w(115) =   0.76141028256526859356393930849227D-02
    w(116) =   0.70344427036681608755685893032552D-02
    w(117) =   0.64505120486899171845442463868748D-02
    w(118) =   0.58626653903523901033648343751367D-02
    w(119) =   0.52712596565634400891303815906251D-02
    w(120) =   0.46766539777779034772638165662478D-02
    w(121) =   0.40792095178254605327114733456293D-02
    w(122) =   0.34792893810051465908910894094105D-02
    w(123) =   0.28772587656289004082883197417581D-02
    w(124) =   0.22734860707492547802810838362671D-02
    w(125) =   0.16683488125171936761028811985672D-02
    w(126) =   0.10622766869538486959954760554099D-02
    w(127) =   0.45645726109586654495731936146574D-03
 
  else if ( n == 255 ) then

    x(  1) =  -0.9999557053175637D+00
    x(  2) =  -0.9997666213120006D+00
    x(  3) =  -0.9994264746801700D+00
    x(  4) =  -0.9989352412846546D+00
    x(  5) =  -0.9982929861369679D+00
    x(  6) =  -0.9974998041266158D+00
    x(  7) =  -0.9965558144351986D+00
    x(  8) =  -0.9954611594800263D+00
    x(  9) =  -0.9942160046166302D+00
    x( 10) =  -0.9928205380219891D+00
    x( 11) =  -0.9912749706303856D+00
    x( 12) =  -0.9895795360859201D+00
    x( 13) =  -0.9877344906997324D+00
    x( 14) =  -0.9857401134074193D+00
    x( 15) =  -0.9835967057247763D+00
    x( 16) =  -0.9813045917010171D+00
    x( 17) =  -0.9788641178690681D+00
    x( 18) =  -0.9762756531927360D+00
    x( 19) =  -0.9735395890106436D+00
    x( 20) =  -0.9706563389768804D+00
    x( 21) =  -0.9676263389983388D+00
    x( 22) =  -0.9644500471687263D+00
    x( 23) =  -0.9611279436992478D+00
    x( 24) =  -0.9576605308459620D+00
    x( 25) =  -0.9540483328338163D+00
    x( 26) =  -0.9502918957773683D+00
    x( 27) =  -0.9463917875982043D+00
    x( 28) =  -0.9423485979390644D+00
    x( 29) =  -0.9381629380746873D+00
    x( 30) =  -0.9338354408193861D+00
    x( 31) =  -0.9293667604313699D+00
    x( 32) =  -0.9247575725138244D+00
    x( 33) =  -0.9200085739127664D+00
    x( 34) =  -0.9151204826116870D+00
    x( 35) =  -0.9100940376230008D+00
    x( 36) =  -0.9049299988763150D+00
    x( 37) =  -0.8996291471035368D+00
    x( 38) =  -0.8941922837208367D+00
    x( 39) =  -0.8886202307074841D+00
    x( 40) =  -0.8829138304815741D+00
    x( 41) =  -0.8770739457726654D+00
    x( 42) =  -0.8711014594913465D+00
    x( 43) =  -0.8649972745957512D+00
    x( 44) =  -0.8587623139550430D+00
    x( 45) =  -0.8523975202098902D+00
    x( 46) =  -0.8459038556299511D+00
    x( 47) =  -0.8392823019683910D+00
    x( 48) =  -0.8325338603134556D+00
    x( 49) =  -0.8256595509371186D+00
    x( 50) =  -0.8186604131408319D+00
    x( 51) =  -0.8115375050983958D+00
    x( 52) =  -0.8042919036959787D+00
    x( 53) =  -0.7969247043693057D+00
    x( 54) =  -0.7894370209380444D+00
    x( 55) =  -0.7818299854374094D+00
    x( 56) =  -0.7741047479470157D+00
    x( 57) =  -0.7662624764170006D+00
    x( 58) =  -0.7583043564914468D+00
    x( 59) =  -0.7502315913291283D+00
    x( 60) =  -0.7420454014216102D+00
    x( 61) =  -0.7337470244087263D+00
    x( 62) =  -0.7253377148914649D+00
    x( 63) =  -0.7168187442422908D+00
    x( 64) =  -0.7081914004129306D+00
    x( 65) =  -0.6994569877396524D+00
    x( 66) =  -0.6906168267460676D+00
    x( 67) =  -0.6816722539434864D+00
    x( 68) =  -0.6726246216288551D+00
    x( 69) =  -0.6634752976803070D+00
    x( 70) =  -0.6542256653503588D+00
    x( 71) =  -0.6448771230567811D+00
    x( 72) =  -0.6354310841711771D+00
    x( 73) =  -0.6258889768052999D+00
    x( 74) =  -0.6162522435951415D+00
    x( 75) =  -0.6065223414828266D+00
    x( 76) =  -0.5967007414963417D+00
    x( 77) =  -0.5867889285271373D+00
    x( 78) =  -0.5767884011056313D+00
    x( 79) =  -0.5667006711746527D+00
    x( 80) =  -0.5565272638608558D+00
    x( 81) =  -0.5462697172441424D+00
    x( 82) =  -0.5359295821251249D+00
    x( 83) =  -0.5255084217906666D+00
    x( 84) =  -0.5150078117775342D+00
    x( 85) =  -0.5044293396341982D+00
    x( 86) =  -0.4937746046808170D+00
    x( 87) =  -0.4830452177674420D+00
    x( 88) =  -0.4722428010304787D+00
    x( 89) =  -0.4613689876474424D+00
    x( 90) =  -0.4504254215900437D+00
    x( 91) =  -0.4394137573756426D+00
    x( 92) =  -0.4283356598171081D+00
    x( 93) =  -0.4171928037711214D+00
    x( 94) =  -0.4059868738849605D+00
    x( 95) =  -0.3947195643418044D+00
    x( 96) =  -0.3833925786045958D+00
    x( 97) =  -0.3720076291585012D+00
    x( 98) =  -0.3605664372520062D+00
    x( 99) =  -0.3490707326366864D+00
    x(100) =  -0.3375222533056927D+00
    x(101) =  -0.3259227452309905D+00
    x(102) =  -0.3142739620993925D+00
    x(103) =  -0.3025776650474256D+00
    x(104) =  -0.2908356223950708D+00
    x(105) =  -0.2790496093784178D+00
    x(106) =  -0.2672214078812731D+00
    x(107) =  -0.2553528061657641D+00
    x(108) =  -0.2434455986019780D+00
    x(109) =  -0.2315015853966777D+00
    x(110) =  -0.2195225723211354D+00
    x(111) =  -0.2075103704381242D+00
    x(112) =  -0.1954667958281108D+00
    x(113) =  -0.1833936693146885D+00
    x(114) =  -0.1712928161892939D+00
    x(115) =  -0.1591660659352477D+00
    x(116) =  -0.1470152519511620D+00
    x(117) =  -0.1348422112737553D+00
    x(118) =  -0.1226487843001178D+00
    x(119) =  -0.1104368145094688D+00
    x(120) =  -0.9820814818444755D-01
    x(121) =  -0.8596463413198061D-01
    x(122) =  -0.7370812340376778D-01
    x(123) =  -0.6144046901642827D-01
    x(124) =  -0.4916352567134998D-01
    x(125) =  -0.3687914947428402D-01
    x(126) =  -0.2458919765472701D-01
    x(127) =  -0.1229552828513332D-01
    x(128) =    0.000000000000000D+00
    x(129) =   0.1229552828513332D-01
    x(130) =   0.2458919765472701D-01
    x(131) =   0.3687914947428402D-01
    x(132) =   0.4916352567134998D-01
    x(133) =   0.6144046901642827D-01
    x(134) =   0.7370812340376778D-01
    x(135) =   0.8596463413198061D-01
    x(136) =   0.9820814818444755D-01
    x(137) =   0.1104368145094688D+00
    x(138) =   0.1226487843001178D+00
    x(139) =   0.1348422112737553D+00
    x(140) =   0.1470152519511620D+00
    x(141) =   0.1591660659352477D+00
    x(142) =   0.1712928161892939D+00
    x(143) =   0.1833936693146885D+00
    x(144) =   0.1954667958281108D+00
    x(145) =   0.2075103704381242D+00
    x(146) =   0.2195225723211354D+00
    x(147) =   0.2315015853966777D+00
    x(148) =   0.2434455986019780D+00
    x(149) =   0.2553528061657641D+00
    x(150) =   0.2672214078812731D+00
    x(151) =   0.2790496093784178D+00
    x(152) =   0.2908356223950708D+00
    x(153) =   0.3025776650474256D+00
    x(154) =   0.3142739620993925D+00
    x(155) =   0.3259227452309905D+00
    x(156) =   0.3375222533056927D+00
    x(157) =   0.3490707326366864D+00
    x(158) =   0.3605664372520062D+00
    x(159) =   0.3720076291585012D+00
    x(160) =   0.3833925786045958D+00
    x(161) =   0.3947195643418044D+00
    x(162) =   0.4059868738849605D+00
    x(163) =   0.4171928037711214D+00
    x(164) =   0.4283356598171081D+00
    x(165) =   0.4394137573756426D+00
    x(166) =   0.4504254215900437D+00
    x(167) =   0.4613689876474424D+00
    x(168) =   0.4722428010304787D+00
    x(169) =   0.4830452177674420D+00
    x(170) =   0.4937746046808170D+00
    x(171) =   0.5044293396341982D+00
    x(172) =   0.5150078117775342D+00
    x(173) =   0.5255084217906666D+00
    x(174) =   0.5359295821251249D+00
    x(175) =   0.5462697172441424D+00
    x(176) =   0.5565272638608558D+00
    x(177) =   0.5667006711746527D+00
    x(178) =   0.5767884011056313D+00
    x(179) =   0.5867889285271373D+00
    x(180) =   0.5967007414963417D+00
    x(181) =   0.6065223414828266D+00
    x(182) =   0.6162522435951415D+00
    x(183) =   0.6258889768052999D+00
    x(184) =   0.6354310841711771D+00
    x(185) =   0.6448771230567811D+00
    x(186) =   0.6542256653503588D+00
    x(187) =   0.6634752976803070D+00
    x(188) =   0.6726246216288551D+00
    x(189) =   0.6816722539434864D+00
    x(190) =   0.6906168267460676D+00
    x(191) =   0.6994569877396524D+00
    x(192) =   0.7081914004129306D+00
    x(193) =   0.7168187442422908D+00
    x(194) =   0.7253377148914649D+00
    x(195) =   0.7337470244087263D+00
    x(196) =   0.7420454014216102D+00
    x(197) =   0.7502315913291283D+00
    x(198) =   0.7583043564914468D+00
    x(199) =   0.7662624764170006D+00
    x(200) =   0.7741047479470157D+00
    x(201) =   0.7818299854374094D+00
    x(202) =   0.7894370209380444D+00
    x(203) =   0.7969247043693057D+00
    x(204) =   0.8042919036959787D+00
    x(205) =   0.8115375050983958D+00
    x(206) =   0.8186604131408319D+00
    x(207) =   0.8256595509371186D+00
    x(208) =   0.8325338603134556D+00
    x(209) =   0.8392823019683910D+00
    x(210) =   0.8459038556299511D+00
    x(211) =   0.8523975202098902D+00
    x(212) =   0.8587623139550430D+00
    x(213) =   0.8649972745957512D+00
    x(214) =   0.8711014594913465D+00
    x(215) =   0.8770739457726654D+00
    x(216) =   0.8829138304815741D+00
    x(217) =   0.8886202307074841D+00
    x(218) =   0.8941922837208367D+00
    x(219) =   0.8996291471035368D+00
    x(220) =   0.9049299988763150D+00
    x(221) =   0.9100940376230008D+00
    x(222) =   0.9151204826116870D+00
    x(223) =   0.9200085739127664D+00
    x(224) =   0.9247575725138244D+00
    x(225) =   0.9293667604313699D+00
    x(226) =   0.9338354408193861D+00
    x(227) =   0.9381629380746873D+00
    x(228) =   0.9423485979390644D+00
    x(229) =   0.9463917875982043D+00
    x(230) =   0.9502918957773683D+00
    x(231) =   0.9540483328338163D+00
    x(232) =   0.9576605308459620D+00
    x(233) =   0.9611279436992478D+00
    x(234) =   0.9644500471687263D+00
    x(235) =   0.9676263389983388D+00
    x(236) =   0.9706563389768804D+00
    x(237) =   0.9735395890106436D+00
    x(238) =   0.9762756531927360D+00
    x(239) =   0.9788641178690681D+00
    x(240) =   0.9813045917010171D+00
    x(241) =   0.9835967057247763D+00
    x(242) =   0.9857401134074193D+00
    x(243) =   0.9877344906997324D+00
    x(244) =   0.9895795360859201D+00
    x(245) =   0.9912749706303856D+00
    x(246) =   0.9928205380219891D+00
    x(247) =   0.9942160046166302D+00
    x(248) =   0.9954611594800263D+00
    x(249) =   0.9965558144351986D+00
    x(250) =   0.9974998041266158D+00
    x(251) =   0.9982929861369679D+00
    x(252) =   0.9989352412846546D+00
    x(253) =   0.9994264746801700D+00
    x(254) =   0.9997666213120006D+00
    x(255) =   0.9999557053175637D+00

    w(  1) =   0.1136736199914808D-03
    w(  2) =   0.2645938711908564D-03
    w(  3) =   0.4156976252681932D-03
    w(  4) =   0.5667579456482639D-03
    w(  5) =   0.7177364780061286D-03
    w(  6) =   0.8686076661194581D-03
    w(  7) =   0.1019347976427318D-02
    w(  8) =   0.1169934372938800D-02
    w(  9) =   0.1320343990022177D-02
    w( 10) =   0.1470554042778403D-02
    w( 11) =   0.1620541799041545D-02
    w( 12) =   0.1770284570660304D-02
    w( 13) =   0.1919759711713187D-02
    w( 14) =   0.2068944619501569D-02
    w( 15) =   0.2217816736754017D-02
    w( 16) =   0.2366353554396287D-02
    w( 17) =   0.2514532614599710D-02
    w( 18) =   0.2662331513971696D-02
    w( 19) =   0.2809727906820460D-02
    w( 20) =   0.2956699508457498D-02
    w( 21) =   0.3103224098519095D-02
    w( 22) =   0.3249279524294296D-02
    w( 23) =   0.3394843704053401D-02
    w( 24) =   0.3539894630372244D-02
    w( 25) =   0.3684410373449933D-02
    w( 26) =   0.3828369084417135D-02
    w( 27) =   0.3971748998634907D-02
    w( 28) =   0.4114528438981242D-02
    w( 29) =   0.4256685819126112D-02
    w( 30) =   0.4398199646792759D-02
    w( 31) =   0.4539048527006180D-02
    w( 32) =   0.4679211165326077D-02
    w( 33) =   0.4818666371065699D-02
    w( 34) =   0.4957393060495050D-02
    w( 35) =   0.5095370260027839D-02
    w( 36) =   0.5232577109391968D-02
    w( 37) =   0.5368992864783177D-02
    w( 38) =   0.5504596902000804D-02
    w( 39) =   0.5639368719565862D-02
    w( 40) =   0.5773287941820301D-02
    w( 41) =   0.5906334322007422D-02
    w( 42) =   0.6038487745332765D-02
    w( 43) =   0.6169728232005295D-02
    w( 44) =   0.6300035940257733D-02
    w( 45) =   0.6429391169346602D-02
    w( 46) =   0.6557774362530328D-02
    w( 47) =   0.6685166110026254D-02
    w( 48) =   0.6811547151944815D-02
    w( 49) =   0.6936898381201466D-02
    w( 50) =   0.7061200846405536D-02
    w( 51) =   0.7184435754724984D-02
    w( 52) =   0.7306584474728122D-02
    w( 53) =   0.7427628539199977D-02
    w( 54) =   0.7547549647934514D-02
    w( 55) =   0.7666329670501377D-02
    w( 56) =   0.7783950648986801D-02
    w( 57) =   0.7900394800708624D-02
    w( 58) =   0.8015644520904983D-02
    w( 59) =   0.8129682385395602D-02
    w( 60) =   0.8242491153216323D-02
    w( 61) =   0.8354053769225508D-02
    w( 62) =   0.8464353366682819D-02
    w( 63) =   0.8573373269798925D-02
    w( 64) =   0.8681096996256795D-02
    w( 65) =   0.8787508259703609D-02
    w( 66) =   0.8892590972213036D-02
    w( 67) =   0.8996329246717397D-02
    w( 68) =   0.9098707399409718D-02
    w( 69) =   0.9199709952114802D-02
    w( 70) =   0.9299321634629343D-02
    w( 71) =   0.9397527387030594D-02
    w( 72) =   0.9494312361953241D-02
    w( 73) =   0.9589661926834022D-02
    w( 74) =   0.9683561666124043D-02
    w( 75) =   0.9775997383468165D-02
    w( 76) =   0.9866955103851452D-02
    w( 77) =   0.9956421075711706D-02
    w( 78) =   0.1004438177301882D-01
    w( 79) =   0.1013082389731963D-01
    w( 80) =   0.1021573437974821D-01
    w( 81) =   0.1029910038300220D-01
    w( 82) =   0.1038090930328312D-01
    w( 83) =   0.1046114877220228D-01
    w( 84) =   0.1053980665865038D-01
    w( 85) =   0.1061687107063194D-01
    w( 86) =   0.1069233035706287D-01
    w( 87) =   0.1076617310953212D-01
    w( 88) =   0.1083838816402652D-01
    w( 89) =   0.1090896460261843D-01
    w( 90) =   0.1097789175511656D-01
    w( 91) =   0.1104515920067912D-01
    w( 92) =   0.1111075676938929D-01
    w( 93) =   0.1117467454379268D-01
    w( 94) =   0.1123690286039691D-01
    w( 95) =   0.1129743231113249D-01
    w( 96) =   0.1135625374477508D-01
    w( 97) =   0.1141335826832922D-01
    w( 98) =   0.1146873724837283D-01
    w( 99) =   0.1152238231236217D-01
    w(100) =   0.1157428534989815D-01
    w(101) =   0.1162443851395193D-01
    w(102) =   0.1167283422205182D-01
    w(103) =   0.1171946515742932D-01
    w(104) =   0.1176432427012535D-01
    w(105) =   0.1180740477805627D-01
    w(106) =   0.1184870016803913D-01
    w(107) =   0.1188820419677619D-01
    w(108) =   0.1192591089179929D-01
    w(109) =   0.1196181455237226D-01
    w(110) =   0.1199590975035326D-01
    w(111) =   0.1202819133101508D-01
    w(112) =   0.1205865441382472D-01
    w(113) =   0.1208729439318107D-01
    w(114) =   0.1211410693911137D-01
    w(115) =   0.1213908799792579D-01
    w(116) =   0.1216223379283022D-01
    w(117) =   0.1218354082449738D-01
    w(118) =   0.1220300587159574D-01
    w(119) =   0.1222062599127671D-01
    w(120) =   0.1223639851961942D-01
    w(121) =   0.1225032107203351D-01
    w(122) =   0.1226239154361966D-01
    w(123) =   0.1227260810948789D-01
    w(124) =   0.1228096922503318D-01
    w(125) =   0.1228747362616942D-01
    w(126) =   0.1229212032952021D-01
    w(127) =   0.1229490863256759D-01
    w(128) =   0.1229583811375833D-01
    w(129) =   0.1229490863256759D-01
    w(130) =   0.1229212032952021D-01
    w(131) =   0.1228747362616942D-01
    w(132) =   0.1228096922503318D-01
    w(133) =   0.1227260810948789D-01
    w(134) =   0.1226239154361966D-01
    w(135) =   0.1225032107203351D-01
    w(136) =   0.1223639851961942D-01
    w(137) =   0.1222062599127671D-01
    w(138) =   0.1220300587159574D-01
    w(139) =   0.1218354082449738D-01
    w(140) =   0.1216223379283022D-01
    w(141) =   0.1213908799792579D-01
    w(142) =   0.1211410693911137D-01
    w(143) =   0.1208729439318107D-01
    w(144) =   0.1205865441382472D-01
    w(145) =   0.1202819133101508D-01
    w(146) =   0.1199590975035326D-01
    w(147) =   0.1196181455237226D-01
    w(148) =   0.1192591089179929D-01
    w(149) =   0.1188820419677619D-01
    w(150) =   0.1184870016803913D-01
    w(151) =   0.1180740477805627D-01
    w(152) =   0.1176432427012535D-01
    w(153) =   0.1171946515742932D-01
    w(154) =   0.1167283422205182D-01
    w(155) =   0.1162443851395193D-01
    w(156) =   0.1157428534989815D-01
    w(157) =   0.1152238231236217D-01
    w(158) =   0.1146873724837283D-01
    w(159) =   0.1141335826832922D-01
    w(160) =   0.1135625374477508D-01
    w(161) =   0.1129743231113249D-01
    w(162) =   0.1123690286039691D-01
    w(163) =   0.1117467454379268D-01
    w(164) =   0.1111075676938929D-01
    w(165) =   0.1104515920067912D-01
    w(166) =   0.1097789175511656D-01
    w(167) =   0.1090896460261843D-01
    w(168) =   0.1083838816402652D-01
    w(169) =   0.1076617310953212D-01
    w(170) =   0.1069233035706287D-01
    w(171) =   0.1061687107063194D-01
    w(172) =   0.1053980665865038D-01
    w(173) =   0.1046114877220228D-01
    w(174) =   0.1038090930328312D-01
    w(175) =   0.1029910038300220D-01
    w(176) =   0.1021573437974821D-01
    w(177) =   0.1013082389731963D-01
    w(178) =   0.1004438177301882D-01
    w(179) =   0.9956421075711706D-02
    w(180) =   0.9866955103851452D-02
    w(181) =   0.9775997383468165D-02
    w(182) =   0.9683561666124043D-02
    w(183) =   0.9589661926834022D-02
    w(184) =   0.9494312361953241D-02
    w(185) =   0.9397527387030594D-02
    w(186) =   0.9299321634629343D-02
    w(187) =   0.9199709952114802D-02
    w(188) =   0.9098707399409718D-02
    w(189) =   0.8996329246717397D-02
    w(190) =   0.8892590972213036D-02
    w(191) =   0.8787508259703609D-02
    w(192) =   0.8681096996256795D-02
    w(193) =   0.8573373269798925D-02
    w(194) =   0.8464353366682819D-02
    w(195) =   0.8354053769225508D-02
    w(196) =   0.8242491153216323D-02
    w(197) =   0.8129682385395602D-02
    w(198) =   0.8015644520904983D-02
    w(199) =   0.7900394800708624D-02
    w(200) =   0.7783950648986801D-02
    w(201) =   0.7666329670501377D-02
    w(202) =   0.7547549647934514D-02
    w(203) =   0.7427628539199977D-02
    w(204) =   0.7306584474728122D-02
    w(205) =   0.7184435754724984D-02
    w(206) =   0.7061200846405536D-02
    w(207) =   0.6936898381201466D-02
    w(208) =   0.6811547151944815D-02
    w(209) =   0.6685166110026254D-02
    w(210) =   0.6557774362530328D-02
    w(211) =   0.6429391169346602D-02
    w(212) =   0.6300035940257733D-02
    w(213) =   0.6169728232005295D-02
    w(214) =   0.6038487745332765D-02
    w(215) =   0.5906334322007422D-02
    w(216) =   0.5773287941820301D-02
    w(217) =   0.5639368719565862D-02
    w(218) =   0.5504596902000804D-02
    w(219) =   0.5368992864783177D-02
    w(220) =   0.5232577109391968D-02
    w(221) =   0.5095370260027839D-02
    w(222) =   0.4957393060495050D-02
    w(223) =   0.4818666371065699D-02
    w(224) =   0.4679211165326077D-02
    w(225) =   0.4539048527006180D-02
    w(226) =   0.4398199646792759D-02
    w(227) =   0.4256685819126112D-02
    w(228) =   0.4114528438981242D-02
    w(229) =   0.3971748998634907D-02
    w(230) =   0.3828369084417135D-02
    w(231) =   0.3684410373449933D-02
    w(232) =   0.3539894630372244D-02
    w(233) =   0.3394843704053401D-02
    w(234) =   0.3249279524294296D-02
    w(235) =   0.3103224098519095D-02
    w(236) =   0.2956699508457498D-02
    w(237) =   0.2809727906820460D-02
    w(238) =   0.2662331513971696D-02
    w(239) =   0.2514532614599710D-02
    w(240) =   0.2366353554396287D-02
    w(241) =   0.2217816736754017D-02
    w(242) =   0.2068944619501569D-02
    w(243) =   0.1919759711713187D-02
    w(244) =   0.1770284570660304D-02
    w(245) =   0.1620541799041545D-02
    w(246) =   0.1470554042778403D-02
    w(247) =   0.1320343990022177D-02
    w(248) =   0.1169934372938800D-02
    w(249) =   0.1019347976427318D-02
    w(250) =   0.8686076661194581D-03
    w(251) =   0.7177364780061286D-03
    w(252) =   0.5667579456482639D-03
    w(253) =   0.4156976252681932D-03
    w(254) =   0.2645938711908564D-03
    w(255) =   0.1136736199914808D-03

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) &
      '  Legal values are 1 through 33, 63, 64, 65, 127 and 255.'
    stop

  end if

  return
end
subroutine legendre_set_x1 ( order, xtab, weight )

!*****************************************************************************80
!
!! LEGENDRE_SET_X1 sets a Gauss-Legendre rule for ( 1 + X ) * F(X) on [-1,1].
!
!  Integration region:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1 + X
!
!  Integral to approximate:
!
!    integral ( -1 <= X <= 1 ) ( 1 + X ) * F(X) dX
!
!  Approximate integral:
!
!    sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4,G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    ORDER must be between 1 and 9.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) weight(order)

  if ( order == 1 ) then

    xtab(1) =  0.333333333333333333333333333333D+00

    weight(1) = 2.0D+00

  else if ( order == 2 ) then

    xtab(1) = -0.289897948556635619639456814941D+00
    xtab(2) =  0.689897948556635619639456814941D+00

    weight(1) =  0.727834473024091322422523991699D+00
    weight(2) =  1.27216552697590867757747600830D+00

  else if ( order == 3 ) then

    xtab(1) = -0.575318923521694112050483779752D+00
    xtab(2) =  0.181066271118530578270147495862D+00
    xtab(3) =  0.822824080974592105208907712461D+00

    weight(1) =  0.279307919605816490135525088716D+00
    weight(2) =  0.916964425438344986775682378225D+00
    weight(3) =  0.803727654955838523088792533058D+00

  else if ( order == 4 ) then

    xtab(1) = -0.720480271312438895695825837750D+00
    xtab(2) = -0.167180864737833640113395337326D+00
    xtab(3) =  0.446313972723752344639908004629D+00
    xtab(4) =  0.885791607770964635613757614892D+00

    weight(1) =  0.124723883800032328695500588386D+00
    weight(2) =  0.519390190432929763305824811559D+00
    weight(3) =  0.813858272041085443165617903743D+00
    weight(4) =  0.542027653725952464833056696312D+00

  else if ( order == 5 ) then

    xtab(1) = -0.802929828402347147753002204224D+00
    xtab(2) = -0.390928546707272189029229647442D+00
    xtab(3) =  0.124050379505227711989974959990D+00
    xtab(4) =  0.603973164252783654928415726409D+00
    xtab(5) =  0.920380285897062515318386619813D+00

    weight(1) =  0.0629916580867691047411692662740D+00
    weight(2) =  0.295635480290466681402532877367D+00
    weight(3) =  0.585547948338679234792151477424D+00
    weight(4) =  0.668698552377478261966702492391D+00
    weight(5) =  0.387126360906606717097443886545D+00

  else if ( order == 6 ) then

    xtab(1) = -0.853891342639482229703747931639D+00
    xtab(2) = -0.538467724060109001833766720231D+00
    xtab(3) = -0.117343037543100264162786683611D+00
    xtab(4) =  0.326030619437691401805894055838D+00
    xtab(5) =  0.703842800663031416300046295008D+00
    xtab(6) =  0.941367145680430216055899446174D+00

    weight(1) =  0.0349532072544381270240692132496D+00
    weight(2) =  0.175820662202035902032706497222D+00
    weight(3) =  0.394644603562621056482338042193D+00
    weight(4) =  0.563170215152795712476307356284D+00
    weight(5) =  0.542169988926074467362761586552D+00
    weight(6) =  0.289241322902034734621817304499D+00

  else if ( order == 7 ) then

    xtab(1) = -0.887474878926155707068695617935D+00
    xtab(2) = -0.639518616526215270024840114382D+00
    xtab(3) = -0.294750565773660725252184459658D+00
    xtab(4) =  0.0943072526611107660028971153047D+00
    xtab(5) =  0.468420354430821063046421216613D+00
    xtab(6) =  0.770641893678191536180719525865D+00
    xtab(7) =  0.955041227122575003782349000858D+00

    weight(1) =  0.0208574488112296163587654972151D+00
    weight(2) =  0.109633426887493901777324193433D+00
    weight(3) =  0.265538785861965879934591955055D+00
    weight(4) =  0.428500262783494679963649011999D+00
    weight(5) =  0.509563589198353307674937943100D+00
    weight(6) =  0.442037032763498409684482945478D+00
    weight(7) =  0.223869453693964204606248453720D+00

  else if ( order == 8 ) then

    xtab(1) = -0.910732089420060298533757956283D+00
    xtab(2) = -0.711267485915708857029562959544D+00
    xtab(3) = -0.426350485711138962102627520502D+00
    xtab(4) = -0.0903733696068532980645444599064D+00
    xtab(5) =  0.256135670833455395138292079035D+00
    xtab(6) =  0.571383041208738483284917464837D+00
    xtab(7) =  0.817352784200412087992517083851D+00
    xtab(8) =  0.964440169705273096373589797925D+00

    weight(1) =  0.0131807657689951954189692640444D+00
    weight(2) =  0.0713716106239448335742111888042D+00
    weight(3) =  0.181757278018795592332221684383D+00
    weight(4) =  0.316798397969276640481632757440D+00
    weight(5) =  0.424189437743720042818124385645D+00
    weight(6) =  0.450023197883549464687088394417D+00
    weight(7) =  0.364476094545494505382889847132D+00
    weight(8) =  0.178203217446223725304862478136D+00

  else if ( order == 9 ) then

    xtab(1) = -0.927484374233581078117671398464D+00
    xtab(2) = -0.763842042420002599615429776011D+00
    xtab(3) = -0.525646030370079229365386614293D+00
    xtab(4) = -0.236234469390588049278459503207D+00
    xtab(5) =  0.0760591978379781302337137826389D+00
    xtab(6) =  0.380664840144724365880759065541D+00
    xtab(7) =  0.647766687674009436273648507855D+00
    xtab(8) =  0.851225220581607910728163628088D+00
    xtab(9) =  0.971175180702246902734346518378D+00

    weight(1) =  0.00872338834309252349019620448007D+00
    weight(2) =  0.0482400171391415162069086091476D+00
    weight(3) =  0.127219285964216005046760427743D+00
    weight(4) =  0.233604781180660442262926091607D+00
    weight(5) =  0.337433287379681397577000079834D+00
    weight(6) =  0.401235236773473158616600898930D+00
    weight(7) =  0.394134968689382820640692081477D+00
    weight(8) =  0.304297020437232650320317215016D+00
    weight(9) =  0.145112014093119485838598391765D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET_X1 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of ORDER = ', order
    stop

  end if

  return
end
subroutine legendre_set_x2 ( order, xtab, weight )

!*****************************************************************************80
!
!! LEGENDRE_SET_X2 sets a Gauss-Legendre rule for ( 1 + X )^2 * F(X) on [-1,1].
!
!  Integration region:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    ( 1 + X )^2
!
!  Integral to approximate:
!
!    integral ( -1 <= X <= 1 ) ( 1 + X )^2 * F(X) dX
!
!  Approximate integral:
!
!    sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    ORDER must be between 1 and 9.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) weight(order)

  if ( order == 1 ) then

    xtab(1) =  0.5D+00

    weight(1) =  2.66666666666666666666666666666D+00

  else if ( order == 2 ) then

    xtab(1) = -0.0883036880224505775998524725910D+00
    xtab(2) =  0.754970354689117244266519139258D+00

    weight(1) =  0.806287056638603444666851075928D+00
    weight(2) =  1.86037961002806322199981559074D+00

  else if ( order == 3 ) then

    xtab(1) = -0.410004419776996766244796955168D+00
    xtab(2) =  0.305992467923296230556472913192D+00
    xtab(3) =  0.854011951853700535688324041976D+00

    weight(1) =  0.239605624068645584091811926047D+00
    weight(2) =  1.16997015407892817602809616291D+00
    weight(3) =  1.25709088851909290654675857771D+00

  else if ( order == 4 ) then

    xtab(1) = -0.591702835793545726606755921586D+00
    xtab(2) = -0.0340945902087350046811467387661D+00
    xtab(3) =  0.522798524896275389882037174551D+00
    xtab(4) =  0.902998901106005341405865485802D+00

    weight(1) =  0.0828179259993445222751812523731D+00
    weight(2) =  0.549071097383384602539010760334D+00
    weight(3) =  1.14767031839371367238662411421D+00
    weight(4) =  0.887107324890223869465850539752D+00

  else if ( order == 5 ) then

    xtab(1) = -0.702108425894032836232448374820D+00
    xtab(2) = -0.268666945261773544694327777841D+00
    xtab(3) =  0.220227225868961343518209179230D+00
    xtab(4) =  0.653039358456608553790815164028D+00
    xtab(5) =  0.930842120163569816951085142737D+00

    weight(1) =  0.0329106016247920636689299329544D+00
    weight(2) =  0.256444805783695354037991444453D+00
    weight(3) =  0.713601289772720001490035944563D+00
    weight(4) =  1.00959169519929190423066348132D+00
    weight(5) =  0.654118274286167343239045863379D+00

  else if ( order == 6 ) then

    xtab(1) = -0.773611232355123732602532012021D+00
    xtab(2) = -0.431362254623427837535325249187D+00
    xtab(3) = -0.0180728263295041680220798103354D+00
    xtab(4) =  0.395126163954217534500188844163D+00
    xtab(5) =  0.736872116684029732026178298518D+00
    xtab(6) =  0.948190889812665614490712786006D+00

    weight(1) =  0.0146486064549543818622276447204D+00
    weight(2) =  0.125762377479560410622810097040D+00
    weight(3) =  0.410316569036929681761034600615D+00
    weight(4) =  0.756617493988329628546336413760D+00
    weight(5) =  0.859011997894245060846045458784D+00
    weight(6) =  0.500309621812647503028212451747D+00

  else if ( order == 7 ) then

    xtab(1) = -0.822366333126005527278634734418D+00
    xtab(2) = -0.547034493182875002223997992852D+00
    xtab(3) = -0.200043026557985860387937545780D+00
    xtab(4) =  0.171995710805880507163425502299D+00
    xtab(5) =  0.518891747903884926692601716998D+00
    xtab(6) =  0.793821941703901970495546427988D+00
    xtab(7) =  0.959734452453198985538996625765D+00

    weight(1) =  0.00714150426951365443207221475404D+00
    weight(2) =  0.0653034050584375560578544725498D+00
    weight(3) =  0.235377690316228918725962815880D+00
    weight(4) =  0.505171029671130381676271523850D+00
    weight(5) =  0.733870426238362032891332767175D+00
    weight(6) =  0.725590596901489156295739839779D+00
    weight(7) =  0.394212014211504966587433032679D+00

  else if ( order == 8 ) then

    xtab(1) = -0.857017929919813794402037235698D+00
    xtab(2) = -0.631543407166567521509503573952D+00
    xtab(3) = -0.339104543648722903660229021109D+00
    xtab(4) = -0.0111941563689783438801237300122D+00
    xtab(5) =  0.316696017045595559454075475675D+00
    xtab(6) =  0.609049663022520165351466780939D+00
    xtab(7) =  0.834198765028697794599267293239D+00
    xtab(8) =  0.967804480896157932935972899807D+00

    weight(1) =  0.00374814227227757804631954025851D+00
    weight(2) =  0.0357961737041152639660521680263D+00
    weight(3) =  0.137974910241879862433949246199D+00
    weight(4) =  0.326515411108352185491692769217D+00
    weight(5) =  0.547577467373226177976217604887D+00
    weight(6) =  0.682278153375510121675529810121D+00
    weight(7) =  0.614544746137780998436053880546D+00
    weight(8) =  0.318231662453524478640851647411D+00

  else if ( order == 9 ) then

    xtab(1) = -0.882491728426548422828684254270D+00
    xtab(2) = -0.694873684026474640346360850039D+00
    xtab(3) = -0.446537143480670863635920316400D+00
    xtab(4) = -0.159388112702326252531544826624D+00
    xtab(5) =  0.141092709224374414981503995427D+00
    xtab(6) =  0.428217823321559204544020866175D+00
    xtab(7) =  0.676480966471850715860378175342D+00
    xtab(8) =  0.863830940812464825046988286026D+00
    xtab(9) =  0.973668228805771018909618924364D+00

    weight(1) =  0.00209009877215570354392734918986D+00
    weight(2) =  0.0205951891648697848186537272448D+00
    weight(3) =  0.0832489326348178964194106978875D+00
    weight(4) =  0.210746247220398685903797568021D+00
    weight(5) =  0.388325022916052063676224499399D+00
    weight(6) =  0.554275165518437673725822282791D+00
    weight(7) =  0.621388553284444032628761363828D+00
    weight(8) =  0.523916296267173054255512857631D+00
    weight(9) =  0.262081160888317771694556320674D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET_X2 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of ORDER = ', order
    stop

  end if

  return
end
function lens_half_2d ( func, center, r, theta1, theta2, order )

!*****************************************************************************80
!
!! LENS_HALF_2D approximates an integral in a circular half lens in 2D.
!
!  Discussion:
!
!    A circular half lens is formed by drawing a circular arc,
!    and joining its endpoints.
!
!    This rule for a circular half lens simply views the region as 
!    a product region, with a coordinate "S" that varies along the
!    radial direction, and a coordinate "T" that varies in the perpendicular
!    direction, and whose extent varies as a function of S.
!
!    A Gauss-Legendre rule is used to construct a product rule that is
!    applied to the region.  The accuracy of the Gauss-Legendre rule,
!    which is valid for a rectangular product rule region, does not
!    apply straightforwardly to this region, since the limits in the
!    "T" coordinate are being handled implicitly.
!
!    This is simply an application of the QMULT_2D algorithm of Stroud.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function of two
!    variables which is to be integrated, of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, real ( kind = 8 ) CENTER(2), the center of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angles of the rays
!    that begin and end the arc.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the Gauss-Legendre rule
!    to be used.  Legal values include 1 through 20, 32 or 64.
!
!    Output, real ( kind = 8 ) LENS_HALF_2D, the approximate value
!    of the integral of the function over the half lens.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) order

  real ( kind = 8 ) ax
  real ( kind = 8 ) ay
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) cx
  real ( kind = 8 ) cy
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lens_half_2d
  real ( kind = 8 ) quad
  real ( kind = 8 ) r
  real ( kind = 8 ) s_length
  real ( kind = 8 ) sx
  real ( kind = 8 ) sy
  real ( kind = 8 ) t_length
  real ( kind = 8 ) tdirx
  real ( kind = 8 ) tdiry
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) thi
  real ( kind = 8 ) tx
  real ( kind = 8 ) ty
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) xtab(order)
!
!  Retrieve the Legendre rule of the given order.
!
  call legendre_set ( order, xtab, weight )
!
!  Determine the points A (on the secant) and B (on the circumference)
!  that will form the "S" direction.
!
  ax = center(1) + r * 0.5D+00 * ( cos ( theta1 ) + cos ( theta2 ) )
  ay = center(2) + r * 0.5D+00 * ( sin ( theta1 ) + sin ( theta2 ) )

  bx = center(1) + r * cos ( 0.5D+00 * ( theta1 + theta2 ) )
  by = center(2) + r * sin ( 0.5D+00 * ( theta1 + theta2 ) )
!
!  Find the length of the line between A and B.
!
  s_length = sqrt ( ( ax - bx ) * ( ax - bx ) & 
                  + ( ay - by ) * ( ay - by ) )

  if ( s_length == 0.0D+00 ) then
    lens_half_2d = 0.0D+00
    return
  end if
!
!  Determine the unit vector in the T direction.
!
  tdirx = ( ay - by ) / s_length
  tdiry = ( bx - ax ) / s_length

  quad = 0.0D+00

  do i = 1, order

    w1 = 0.5D+00 * s_length * weight(i)
!
!  Map the quadrature point to an S coordinate.
!
    sx = ( ( 1.0D+00 - xtab(i) ) * ax   &
         + ( 1.0D+00 + xtab(i) ) * bx ) &
         /   2.0D+00
    sy = ( ( 1.0D+00 - xtab(i) ) * ay   &
         + ( 1.0D+00 + xtab(i) ) * by ) &
         /   2.0D+00
!
!  Determine the length of the line in the T direction, from the
!  S axis to the circle circumference.
!
    thi = sqrt ( ( r - 0.25D+00 * ( 1.0D+00 - xtab(i) ) * s_length ) &
                         *        ( 1.0D+00 - xtab(i) ) * s_length )
! 
!  Determine the maximum and minimum T coordinates by going
!  up and down in the T direction from the S axis.
!
    cx = sx + tdirx * thi
    cy = sy + tdiry * thi
    dx = sx - tdirx * thi
    dy = sy - tdiry * thi
!
!  Find the length of the T direction.
!
    t_length = sqrt ( ( cx - dx ) * ( cx - dx ) &
                    + ( cy - dy ) * ( cy - dy ) )

    do j = 1, order

      w2 = 0.5D+00 * t_length * weight(j)
!
!  Map the quadrature point to a T coordinate.
!
      tx = ( ( 1.0D+00 - xtab(j) ) * cx   &
           + ( 1.0D+00 + xtab(j) ) * dx ) &
           /   2.0D+00
      ty = ( ( 1.0D+00 - xtab(j) ) * cy   &
           + ( 1.0D+00 + xtab(j) ) * dy ) &
           /   2.0D+00

      quad = quad + w1 * w2 * func ( tx, ty )

    end do

  end do

  lens_half_2d = quad

  return
end
function lens_half_area_2d ( r, theta1, theta2 )

!*****************************************************************************80
!
!! LENS_HALF_AREA_2D returns the area of a circular half lens in 2D.
!
!  Discussion:
!
!    A circular half lens is formed by drawing a circular arc, 
!    and joining its endpoints.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angles of the rays
!    that begin and end the arc.
!
!    Output, real ( kind = 8 ) LENS_HALF_AREA_2D, the area of the half lens.
!
  implicit none

  real ( kind = 8 ) circle_sector_area_2d
  real ( kind = 8 ) circle_triangle_area_2d
  real ( kind = 8 ) lens_half_area_2d
  real ( kind = 8 ) r
  real ( kind = 8 ) sector
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) triangle

  sector = circle_sector_area_2d ( r, theta1, theta2 )
  triangle = circle_triangle_area_2d ( r, theta1, theta2 )
  lens_half_area_2d = sector - triangle

  return
end
function lens_half_h_area_2d ( r, h )

!*****************************************************************************80
!
!! LENS_HALF_H_AREA_2D returns the area of a circular half lens in 2D.
!
!  Discussion:
!
!    A circular half lens is formed by drawing a circular arc, and joining 
!    its endpoints.
!
!    This particular half lens is described by the "height" of the region.  
!    In other words, the half lens is the region that would be submerged 
!    if a circle of radius R were standing in water of depth H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) H, the height of the half lens region.
!
!    Output, real ( kind = 8 ) LENS_HALF_H_AREA_2D, the area of the half lens.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) area
  real ( kind = 8 ) h
  real ( kind = 8 ) half_width
  real ( kind = 8 ) lens_half_h_area_2d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) sector
  real ( kind = 8 ) triangle

  if ( h <= 0.0D+00 ) then

    area = 0.0D+00

  else if ( 2.0D+00 * r <= h ) then

    area = pi * r * r

  else

    half_width = sqrt ( h * ( 2.0D+00 * r - h ) )
    angle = 2.0D+00 * atan2 ( half_width, r - h )
    sector = r * r * angle / 2.0D+00
    triangle = ( r - h ) * half_width
    area = sector - triangle

  end if

  lens_half_h_area_2d = area

  return
end
function lens_half_w_area_2d ( r, w )

!*****************************************************************************80
!
!! LENS_HALF_W_AREA_2D returns the area of a circular half lens in 2D.
!
!  Discussion:
!
!    A half lens is formed by drawing a circular arc, and joining its endpoints.
!    This half lens is described by the "width" of the region.  In other words,
!    it is the portion of the circle under water if the width
!    of the water surface is W.  There are two possible values for this
!    area, A and (PI*R*R-A).  The routine returns the smaller of the 
!    two values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) W, the width of the half lens region.
!
!    Output, real ( kind = 8 ) LENS_HALF_W_AREA_2D, the area of the half lens.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) area
  real ( kind = 8 ) h
  real ( kind = 8 ) half_width
  real ( kind = 8 ) lens_half_w_area_2d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) sector
  real ( kind = 8 ) triangle
  real ( kind = 8 ) w

  if ( w <= 0.0D+00 ) then

    area = 0.0D+00

  else if ( 2.0D+00 * r <= w ) then

    area = 0.5D+00 * pi * r * r

  else

    half_width = 0.5D+00 * w
    h = r - sqrt ( r * r - half_width * half_width )
    angle = 2.0D+00 * atan2 ( half_width, r - h )
    sector = r * r * angle / 2.0D+00
    triangle = ( r - h ) * half_width
    area = sector - triangle

  end if

  lens_half_w_area_2d = area

  return
end
subroutine monomial_value ( dim_num, point_num, x, expon, value )

!*****************************************************************************80
!
!! MONOMIAL_VALUE evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points at which the
!    monomial is to be evaluated.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the point coordinates.
!
!    Input, integer ( kind = 4 ) EXPON(DIM_NUM), the exponents.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the value of the monomial.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) expon(dim_num)
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  value(1:point_num) = 1.0D+00

  do dim = 1, dim_num
    if ( 0 /= expon(dim) ) then
      value(1:point_num) = value(1:point_num) * x(dim,1:point_num)**expon(dim)
    end if
  end do

  return
end
subroutine octahedron_unit_nd ( func, n, result )

!*****************************************************************************80
!
!! OCTAHEDRON_UNIT_ND approximates integrals in the unit octahedron in ND.
!
!  Integration region:
!
!    sum ( abs ( X(1:N) ) ) <= 1.
!
!  Discussion:
!
!    A 2*N point 3rd degree formula is used, Stroud number GN:3-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which is to be integrated, of the form:
!
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the octahedron.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) octahedron_unit_volume_nd
  real ( kind = 8 ) quad
  real ( kind = 8 ) r
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x(n)

  w = 1.0D+00 / real ( 2 * n, kind = 8 )
  r = sqrt ( real ( 2 * n, kind = 8 ) &
    / real ( ( n + 1 ) * ( n + 2 ), kind = 8 ) )

  x(1:n) = 0.0D+00

  quad = 0.0D+00
  do i = 1, n
    x(i) = r
    do j = 1, 2
      quad = quad + w * func ( n, x )
      x(i) = -x(i)
    end do
    x(i) = 0.0D+00
  end do

  volume = octahedron_unit_volume_nd ( n )
  result = quad * volume

  return
end
function octahedron_unit_volume_nd ( n )

!*****************************************************************************80
!
!! OCTAHEDRON_UNIT_VOLUME_ND returns the volume of the unit octahedron in ND.
!
!  Integration region:
!
!    sum ( abs ( X(1:N) ) ) <= 1.
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
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) OCTAHEDRON_UNIT_VOLUME_ND, the volume of
!    the unit octahedron.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) octahedron_unit_volume_nd
  real ( kind = 8 ) volume

  volume = 1.0D+00
  do i = 1, n
    volume = volume * 2.0D+00 / real ( i, kind = 8 )
  end do

  octahedron_unit_volume_nd = volume

  return
end
function parallelipiped_volume_3d ( x, y, z )

!*****************************************************************************80
!
!! PARALLELIPIPED_VOLUME_3D returns the volume of a parallelipiped in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(4), Y(4), Z(4), the coordinates of one corner
!    of the parallelipiped, and its 3 immediate neighbors.
!
!    Output, real ( kind = 8 ) PARALLELIPIPED_VOLUME_3D, the volume of
!    the parallelipiped.
!
  implicit none

  real ( kind = 8 ) parallelipiped_volume_3d
  real ( kind = 8 ) x(4)
  real ( kind = 8 ) y(4)
  real ( kind = 8 ) z(4)

  parallelipiped_volume_3d = abs ( &
    ( z(2) - z(1) ) * ( y(4) * x(3) - y(3) * x(4) ) + &
    ( z(3) - z(1) ) * ( x(4) * y(2) - x(2) * y(4) ) + &
    ( z(4) - z(1) ) * ( x(2) * y(3) - x(3) * y(2) ) + &
    ( z(3) - z(2) ) * ( y(4) * x(1) - y(1) * x(4) ) + &
    ( z(4) - z(2) ) * ( x(3) * y(1) - x(1) * y(3) ) + &
    ( z(4) - z(3) ) * ( x(1) * y(2) - x(2) * y(1) ) )

  return
end
function parallelipiped_volume_nd ( n, v )

!*****************************************************************************80
!
!! PARALLELIPIPED_VOLUME_ND returns the volume of a parallelipiped in ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) V(N,N+1), the
!    N+1 columns of V contains the N coordinates of one of the
!    "corners" of the parallelipiped.
!
!    Output, real ( kind = 8 ) PARALLELIPIPED_VOLUME_ND, the volume of
!    the parallelipiped.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 8 ) parallelipiped_volume_nd
  integer ( kind = 4 ) pivot(n)
  real ( kind = 8 ) v(n,n+1)
  real ( kind = 8 ) w(n,n)
!
!  Compute the volume of the N-dimensional parallelipiped.
!
  do i = 1, n
    w(i,1:n) = v(i,2:n+1) - v(i,1)
  end do

  call r8ge_fa ( n, w, pivot, info )

  if ( info /= 0 ) then
    parallelipiped_volume_nd = 0.0D+00
    return
  end if

  call r8ge_det ( n, w, pivot, det )

  parallelipiped_volume_nd = abs ( det )

  return
end
subroutine polygon_1_2d ( n, x, y, result )

!*****************************************************************************80
!
!! POLYGON_1_2D integrates the function 1 over a polygon in 2D.
!
!  Integration region:
!
!    The polygon bounded by the points (X(1:N), Y(1:N)).
!
!  Formula:
!
!    INTEGRAL = 0.5 * sum ( 1 <= I <= N )
!      ( X(I) + X(I-1) ) * ( Y(I) - Y(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!  Discussion:
!
!    The integral of 1 over a polygon is the area of the polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    Volume 96, Number 2, February 1989, pages 131-132.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in counter-clockwise order.
!
!    Output, real ( kind = 8 ) RESULT, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  real ( kind = 8 ) result
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYGON_1_2D - Warning!'
    write ( *, '(a)' ) '  The number of vertices must be at least 3.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    return
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result + 0.5D+00 * ( x(i) + x(im1) ) * ( y(i) - y(im1) )

  end do

  return
end
subroutine polygon_x_2d ( n, x, y, result )

!*****************************************************************************80
!
!! POLYGON_X_2D integrates the function X over a polygon in 2D.
!
!  Integration region:
!
!    The polygon bounded by the points (X(1:N), Y(1:N)).
!
!  Formula:
!
!    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
!      ( X(I)^2 + X(I) * X(I-1) + X(I-1)^2 ) * ( Y(I) - Y(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    Volume 96, Number 2, February 1989, pages 131-132.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in counter-clockwise order.
!
!    Output, real ( kind = 8 ) RESULT, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  real ( kind = 8 ) result
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYGON_X_2D - Warning!'
    write ( *, '(a)' ) '  The number of vertices must be at least 3.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    return
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result + ( x(i) * x(i) + x(i) * x(im1) + x(im1) * x(im1) ) &
      * ( y(i) - y(im1) )

  end do

  result = result / 6.0D+00

  return
end
subroutine polygon_xx_2d ( n, x, y, result )

!*****************************************************************************80
!
!! POLYGON_XX_2D integrates the function X*X over a polygon in 2D.
!
!  Integration region:
!
!    The polygon bounded by the points (X(1:N), Y(1:N)).
!
!  Formula:
!
!    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
!      ( X(I)^3 + X(I)^2 * X(I-1) + X(I) * X(I-1)^2 + X(I-1)^3 )
!      * ( Y(I) - Y(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    Volume 96, Number 2, February 1989, pages 131-132.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in
!    counter-clockwise order.
!
!    Output, real ( kind = 8 ) RESULT, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  real ( kind = 8 ) result
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYGON_XX_2D - Warning!'
    write ( *, '(a)' ) '  The number of vertices must be at least 3.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    return
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result + ( x(i)**3 + x(i) * x(i) * x(im1) &
      + x(i) * x(im1) * x(im1) + x(im1)**3 ) * ( y(i) - y(im1) )

  end do

  result = result / 12.0D+00

  return
end
subroutine polygon_xy_2d ( n, x, y, result )

!*****************************************************************************80
!
!! POLYGON_XY_2D integrates the function X*Y over a polygon in 2D.
!
!  Integration region:
!
!    The polygon bounded by the points (X(1:N), Y(1:N)).
!
!  Formula:
!
!    INTEGRAL = (1/24) * sum ( 1 <= I <= N )
!      ( Y(I)   * ( 3 * X(I)^2 + 2 * X(I) * X(I-1) +     X(I-1)^2 )
!      + Y(I-1) * (     X(I)^2 + 2 * X(I) * X(I-1) + 3 * X(I-1)^2 ) )
!      * ( Y(I) - Y(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    Volume 96, Number 2, February 1989, pages 131-132.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in
!    counter-clockwise order.
!
!    Output, real ( kind = 8 ) RESULT, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  real ( kind = 8 ) result
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYGON_XY_2D - Warning!'
    write ( *, '(a)' ) '  The number of vertices must be at least 3.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    return
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result + ( &
      y(i) * ( 3.0D+00 * x(i) * x(i) + 2.0D+00 * x(i) * x(im1) + x(im1)**2 ) &
      + y(im1) * ( x(i) * x(i) + 2.0D+00 * x(i) * x(im1) &
      + 3.0D+00 * x(im1)**2 ) ) * ( y(i) - y(im1) )

  end do

  result = result / 24.0D+00

  return
end
subroutine polygon_y_2d ( n, x, y, result )

!*****************************************************************************80
!
!! POLYGON_Y_2D integrates the function Y over a polygon in 2D.
!
!  Integration region:
!
!    The polygon bounded by the points (X(1:N), Y(1:N)).
!
!  Formula:
!
!    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
!      - ( Y(I)^2 + Y(I) * Y(I-1) + Y(I-1)^2 ) * ( X(I) - X(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    Volume 96, Number 2, February 1989, pages 131-132.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in
!    counter-clockwise order.
!
!    Output, real ( kind = 8 ) RESULT, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  real ( kind = 8 ) result
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYGON_Y_2D - Warning!'
    write ( *, '(a)' ) '  The number of vertices must be at least 3.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    return
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result - ( y(i) * y(i) + y(i) * y(im1) + y(im1) * y(im1) ) &
      * ( x(i) - x(im1) )

  end do

  result = result / 6.0D+00

  return
end
subroutine polygon_yy_2d ( n, x, y, result )

!*****************************************************************************80
!
!! POLYGON_YY_2D integrates the function Y*Y over a polygon in 2D.
!
!  Integration region:
!
!    The polygon bounded by the points (X(1:N), Y(1:N)).
!
!  Formula:
!
!    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
!      - ( Y(I)^3 + Y(I)^2 * Y(I-1) + Y(I) * Y(I-1)^2 + Y(I-1)^3 )
!      * ( X(I) - X(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    Volume 96, Number 2, February 1989, pages 131-132.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in
!    counter-clockwise order.
!
!    Output, real ( kind = 8 ) RESULT, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) im1
  real ( kind = 8 ) result
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYGON_YY_2D - Warning!'
    write ( *, '(a)' ) '  The number of polygonal vertices must be '
    write ( *, '(a,i8)' ) '  at least 3, but the input polygon has N = ', n
    return
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result - ( &
        y(i)**3 & 
      + y(i) * y(i) * y(im1) &
      + y(i) * y(im1) * y(im1) &
      + y(im1)**3 &
    ) * ( x(i) - x(im1) )

  end do

  result = result / 12.0D+00

  return
end
subroutine pyramid_unit_o01_3d ( func, result )

!*****************************************************************************80
!
!! PYRAMID_UNIT_O01_3D approximates an integral inside the unit pyramid in 3D.
!
!  Discussion:
!
!    A 1 point degree 1 formula is used.
!
!    The (X,Y,Z) integration region can be represented as:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y) 
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
!    radius of the square diminishes, and when Z reaches 1, the square has 
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  real ( kind = 8 ), external :: func
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
!  Quadrature.
!
  quad = 0.0D+00

  x = 0.0D+00
  y = 0.0D+00
  z = 1.0D+00 / 4.0D+00
  w = 1.0D+00

  quad = quad + w * func ( x, y, z )
!
!  Volume.
!
  volume = pyramid_unit_volume_3d ( )
!
!  Result.
!
  result = quad * volume

  return
end
subroutine pyramid_unit_o05_3d ( func, result )

!*****************************************************************************80
!
!! PYRAMID_UNIT_O05_3D approximates an integral inside the unit pyramid in 3D.
!
!  Discussion:
!
!    A 5 point formula is used.
!
!    The (X,Y,Z) integration region can be represented as:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y) 
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
!    radius of the square diminishes, and when Z reaches 1, the square has 
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 5

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ), dimension ( order ) :: w = (/ &
   0.21093750000000000000D+00, &
   0.21093750000000000000D+00, &
   0.21093750000000000000D+00, &
   0.21093750000000000000D+00, &
   0.15625000000000000000D+00 /)
  real ( kind = 8 ), dimension ( order ) :: x = (/ &
  -0.48686449556014765641D+00, &
   0.48686449556014765641D+00, &
   0.48686449556014765641D+00, &
  -0.48686449556014765641D+00, &
   0.00000000000000000000D+00 /)
  real ( kind = 8 ), dimension ( order ) :: y = (/ &
  -0.48686449556014765641D+00, &
  -0.48686449556014765641D+00, &
   0.48686449556014765641D+00, &
   0.48686449556014765641D+00, &
   0.00000000000000000000D+00 /)
  real ( kind = 8 ), dimension ( order ) :: z = (/ &
   0.16666666666666666667D+00, &
   0.16666666666666666667D+00, &
   0.16666666666666666667D+00, &
   0.16666666666666666667D+00, &
   0.70000000000000000000D+00 /)
!
!  Quadrature.
!
  quad = 0.0D+00
  do i = 1, order
    quad = quad + w(i) * func ( x(i), y(i), z(i) )
  end do
!
!  Volume.
!
  volume = pyramid_unit_volume_3d ( )
!
!  Result.
!
  result = quad * volume

  return
end
subroutine pyramid_unit_o06_3d ( func, result )

!*****************************************************************************80
!
!! PYRAMID_UNIT_O06_3D approximates an integral inside the unit pyramid in 3D.
!
!  Discussion:
!
!    A 6 point formula is used.
!
!    The (X,Y,Z) integration region can be represented as:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y) 
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
!    radius of the square diminishes, and when Z reaches 1, the square has 
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 6

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ), dimension ( order ) :: w = (/ &
   0.21000000000000000000D+00, &
   0.21000000000000000000D+00, &
   0.21000000000000000000D+00, &
   0.21000000000000000000D+00, &
   0.06000000000000000000D+00, &
   0.10000000000000000000D+00 /)
  real ( kind = 8 ), dimension ( order ) :: x = (/ &
  -0.48795003647426658968D+00, &
   0.48795003647426658968D+00, &
   0.48795003647426658968D+00, &
  -0.48795003647426658968D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00 /)
  real ( kind = 8 ), dimension ( order ) :: y = (/ &
  -0.48795003647426658968D+00, &
  -0.48795003647426658968D+00, &
   0.48795003647426658968D+00, &
   0.48795003647426658968D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00 /)
  real ( kind = 8 ), dimension ( order ) :: z = (/ &
   0.16666666666666666667D+00, &
   0.16666666666666666667D+00, &
   0.16666666666666666667D+00, &
   0.16666666666666666667D+00, &
   0.58333333333333333333D+00, &
   0.75000000000000000000D+00 /)
!
!  Quadrature.
!
  quad = 0.0D+00
  do i = 1, order
    quad = quad + w(i) * func ( x(i), y(i), z(i) )
  end do
!
!  Volume.
!
  volume = pyramid_unit_volume_3d ( )
!
!  Result.
!
  result = quad * volume

  return
end
subroutine pyramid_unit_o08_3d ( func, result )

!*****************************************************************************80
!
!! PYRAMID_UNIT_O08_3D approximates an integral inside the unit pyramid in 3D.
!
!  Discussion:
!
!    An 8 point formula is used.
!
!    The (X,Y,Z) integration region can be represented as:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y) 
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
!    radius of the square diminishes, and when Z reaches 1, the square has 
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 8

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ), dimension ( order ) :: w = (/ &
   0.075589411559869072938D+00, &
   0.075589411559869072938D+00, &
   0.075589411559869072938D+00, &
   0.075589411559869072938D+00, &
   0.17441058844013092706D+00, &
   0.17441058844013092706D+00, &
   0.17441058844013092706D+00, &
   0.17441058844013092706D+00 /)
  real ( kind = 8 ), dimension ( order ) :: x = (/ &
  -0.26318405556971359557D+00, &
   0.26318405556971359557D+00, &
   0.26318405556971359557D+00, &
  -0.26318405556971359557D+00, &
  -0.50661630334978742377D+00, &
   0.50661630334978742377D+00, &
   0.50661630334978742377D+00, &
  -0.50661630334978742377D+00 /)
  real ( kind = 8 ), dimension ( order ) :: y = (/ &
  -0.26318405556971359557D+00, &
  -0.26318405556971359557D+00, &
   0.26318405556971359557D+00, &
   0.26318405556971359557D+00, &
  -0.50661630334978742377D+00, &
  -0.50661630334978742377D+00, &
   0.50661630334978742377D+00, &
   0.50661630334978742377D+00 /)
  real ( kind = 8 ), dimension ( order ) :: z = (/ &
   0.54415184401122528880D+00, &
   0.54415184401122528880D+00, &
   0.54415184401122528880D+00, &
   0.54415184401122528880D+00, &
   0.12251482265544137787D+00, &
   0.12251482265544137787D+00, &
   0.12251482265544137787D+00, &
   0.12251482265544137787D+00 /)
!
!  Quadrature.
!
  quad = 0.0D+00
  do i = 1, order
    quad = quad + w(i) * func ( x(i), y(i), z(i) )
  end do
!
!  Volume.
!
  volume = pyramid_unit_volume_3d ( )
!
!  Result.
!
  result = quad * volume

  return
end
subroutine pyramid_unit_o08b_3d ( func, result )

!*****************************************************************************80
!
!! PYRAMID_UNIT_O08B_3D approximates an integral inside the unit pyramid in 3D.
!
!  Discussion:
!
!    An 8 point formula is used.
!
!    The (X,Y,Z) integration region can be represented as:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y) 
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
!    radius of the square diminishes, and when Z reaches 1, the square has 
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 8

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ), dimension ( order ) :: w = (/ &
   0.16438287736328777572D+00, &
   0.16438287736328777572D+00, &
   0.16438287736328777572D+00, &
   0.16438287736328777572D+00, &
   0.085617122636712224276D+00, &
   0.085617122636712224276D+00, &
   0.085617122636712224276D+00, &
   0.085617122636712224276D+00 /)
  real ( kind = 8 ), dimension ( order ) :: x = (/ &
  -0.51197009372656270107D+00, &
   0.51197009372656270107D+00, &
   0.51197009372656270107D+00, &
  -0.51197009372656270107D+00, &
  -0.28415447557052037456D+00, &
   0.28415447557052037456D+00, &
   0.28415447557052037456D+00, &
  -0.28415447557052037456D+00 /)
  real ( kind = 8 ), dimension ( order ) :: y = (/ &
  -0.51197009372656270107D+00, &
  -0.51197009372656270107D+00, &
   0.51197009372656270107D+00, &
   0.51197009372656270107D+00, &
  -0.28415447557052037456D+00, &
  -0.28415447557052037456D+00, &
   0.28415447557052037456D+00, &
   0.28415447557052037456D+00 /)
  real ( kind = 8 ), dimension ( order ) :: z = (/ &
   0.11024490204163285720D+00, &
   0.11024490204163285720D+00, &
   0.11024490204163285720D+00, &
   0.11024490204163285720D+00, &
   0.518326526529795714229D+00, &
   0.518326526529795714229D+00, &
   0.518326526529795714229D+00, &
   0.518326526529795714229D+00 /)
!
!  Quadrature.
!
  quad = 0.0D+00
  do i = 1, order
    quad = quad + w(i) * func ( x(i), y(i), z(i) )
  end do
!
!  Volume.
!
  volume = pyramid_unit_volume_3d ( )
!
!  Result.
!
  result = quad * volume

  return
end
subroutine pyramid_unit_o09_3d ( func, result )

!*****************************************************************************80
!
!! PYRAMID_UNIT_O09_3D approximates an integral inside the unit pyramid in 3D.
!
!  Discussion:
!
!    A 9 point formula is used.
!
!    The (X,Y,Z) integration region can be represented as:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y) 
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
!    radius of the square diminishes, and when Z reaches 1, the square has 
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 9

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ), dimension ( order ) :: w = (/ &
   0.13073389672275944791D+00, &
   0.13073389672275944791D+00, &
   0.13073389672275944791D+00, &
   0.13073389672275944791D+00, &
   0.10989110327724055209D+00, &
   0.10989110327724055209D+00, &
   0.10989110327724055209D+00, &
   0.10989110327724055209D+00, &
   0.03750000000000000000D+00 /)
  real ( kind = 8 ), dimension ( order ) :: x = (/ &
  -0.52966422253852215131D+00, &
   0.52966422253852215131D+00, &
   0.52966422253852215131D+00, &
  -0.52966422253852215131D+00, &
  -0.34819753825720418039D+00, &
   0.34819753825720418039D+00, &
   0.34819753825720418039D+00, &
  -0.34819753825720418039D+00, &
   0.00000000000000000000D+00 /)
  real ( kind = 8 ), dimension ( order ) :: y = (/ &
  -0.52966422253852215131D+00, &
  -0.52966422253852215131D+00, &
   0.52966422253852215131D+00, &
   0.52966422253852215131D+00, &
  -0.34819753825720418039D+00, &
  -0.34819753825720418039D+00, &
   0.34819753825720418039D+00, &
   0.34819753825720418039D+00, &
   0.00000000000000000000D+00 /)
  real ( kind = 8 ), dimension ( order ) :: z = (/ &
   0.08176876558246862335D+00, &
   0.08176876558246862335D+00, &
   0.08176876558246862335D+00, &
   0.08176876558246862335D+00, &
   0.400374091560388519511D+00, &
   0.400374091560388519511D+00, &
   0.400374091560388519511D+00, &
   0.400374091560388519511D+00, &
   0.83333333333333333333D+00 /)
!
!  Quadrature.
!
  quad = 0.0D+00
  do i = 1, order
    quad = quad + w(i) * func ( x(i), y(i), z(i) )
  end do
!
!  Volume.
!
  volume = pyramid_unit_volume_3d ( )
!
!  Result.
!
  result = quad * volume

  return
end
subroutine pyramid_unit_o13_3d ( func, result )

!*****************************************************************************80
!
!! PYRAMID_UNIT_O13_3D approximates an integral inside the unit pyramid in 3D.
!
!  Discussion:
!
!    A 13 point formula is used.
!
!    The (X,Y,Z) integration region can be represented as:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y) 
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
!    radius of the square diminishes, and when Z reaches 1, the square has 
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 13

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ), dimension ( order ) :: w = (/ &
   0.063061594202898550725D+00, &
   0.063061594202898550725D+00, &
   0.063061594202898550725D+00, &
   0.063061594202898550725D+00, &
   0.042101946815575556199D+00, &
   0.042101946815575556199D+00, &
   0.042101946815575556199D+00, &
   0.042101946815575556199D+00, &
   0.13172030707666776585D+00, &
   0.13172030707666776585D+00, &
   0.13172030707666776585D+00, &
   0.13172030707666776585D+00, &
   0.05246460761943250889D+00 /)
  real ( kind = 8 ), dimension ( order ) :: x = (/ &
  -0.38510399211870384331D+00, &
   0.38510399211870384331D+00, &
   0.38510399211870384331D+00, &
  -0.38510399211870384331D+00, &
  -0.40345831960728204766D+00, &
   0.40345831960728204766D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
  -0.53157877436961973359D+00, &
   0.53157877436961973359D+00, &
   0.53157877436961973359D+00, &
  -0.53157877436961973359D+00, &
   0.00000000000000000000D+00 /)
  real ( kind = 8 ), dimension ( order ) :: y = (/ &
  -0.38510399211870384331D+00, &
  -0.38510399211870384331D+00, &
   0.38510399211870384331D+00, &
   0.38510399211870384331D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
  -0.40345831960728204766D+00, &
   0.40345831960728204766D+00, &
  -0.53157877436961973359D+00, &
  -0.53157877436961973359D+00, &
   0.53157877436961973359D+00, &
   0.53157877436961973359D+00, &
   0.00000000000000000000D+00 /)
  real ( kind = 8 ), dimension ( order ) :: z = (/ &
  0.428571428571428571429D+00, &
  0.428571428571428571429D+00, &
  0.428571428571428571429D+00, &
  0.428571428571428571429D+00, &
  0.33928571428571428571D+00, &
  0.33928571428571428571D+00, &
  0.33928571428571428571D+00, &
  0.33928571428571428571D+00, &
  0.08496732026143790850D+00, &
  0.08496732026143790850D+00, &
  0.08496732026143790850D+00, &
  0.08496732026143790850D+00, &
  0.76219701803768503595D+00 /)
!
!  Quadrature.
!
  quad = 0.0D+00
  do i = 1, order
    quad = quad + w(i) * func ( x(i), y(i), z(i) )
  end do
!
!  Volume.
!
  volume = pyramid_unit_volume_3d ( )
!
!  Result.
!
  result = quad * volume

  return
end
subroutine pyramid_unit_o18_3d ( func, result )

!*****************************************************************************80
!
!! PYRAMID_UNIT_O18_3D approximates an integral inside the unit pyramid in 3D.
!
!  Discussion:
!
!    An 18 point formula is used.
!
!    The (X,Y,Z) integration region can be represented as:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y) 
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
!    radius of the square diminishes, and when Z reaches 1, the square has 
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 18

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ), dimension ( order ) :: w = (/ &
   0.023330065296255886709D+00, &
   0.037328104474009418735D+00, &
   0.023330065296255886709D+00, &
   0.037328104474009418735D+00, &
   0.059724967158415069975D+00, &
   0.037328104474009418735D+00, &
   0.023330065296255886709D+00, &
   0.037328104474009418735D+00, &
   0.023330065296255886709D+00, &
   0.05383042853090460712D+00, &
   0.08612868564944737139D+00, &
   0.05383042853090460712D+00, &
   0.08612868564944737139D+00, &
   0.13780589703911579422D+00, &
   0.08612868564944737139D+00, &
   0.05383042853090460712D+00, &
   0.08612868564944737139D+00, &
   0.05383042853090460712D+00 /)
  real ( kind = 8 ), dimension ( order ) :: x = (/ &
  -0.35309846330877704481D+00, &
   0.00000000000000000000D+00, &
   0.35309846330877704481D+00, &
  -0.35309846330877704481D+00, &
   0.00000000000000000000D+00, &
   0.35309846330877704481D+00, &
  -0.35309846330877704481D+00, &
   0.00000000000000000000D+00, &
   0.35309846330877704481D+00, &
  -0.67969709567986745790D+00, &
   0.00000000000000000000D+00, &
   0.67969709567986745790D+00, &
  -0.67969709567986745790D+00, &
   0.00000000000000000000D+00, &
   0.67969709567986745790D+00, &
  -0.67969709567986745790D+00, &
   0.00000000000000000000D+00, &
   0.67969709567986745790D+00 /)
  real ( kind = 8 ), dimension ( order ) :: y = (/ &
  -0.35309846330877704481D+00, &
  -0.35309846330877704481D+00, &
  -0.35309846330877704481D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.35309846330877704481D+00, &
   0.35309846330877704481D+00, &
   0.35309846330877704481D+00, &
  -0.67969709567986745790D+00, &
  -0.67969709567986745790D+00, &
  -0.67969709567986745790D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.67969709567986745790D+00, &
   0.67969709567986745790D+00, &
   0.67969709567986745790D+00 /)
  real ( kind = 8 ), dimension ( order ) :: z = (/ &
  0.544151844011225288800D+00, &
  0.544151844011225288800D+00, &
  0.544151844011225288800D+00, &
  0.544151844011225288800D+00, &
  0.544151844011225288800D+00, &
  0.544151844011225288800D+00, &
  0.544151844011225288800D+00, &
  0.544151844011225288800D+00, &
  0.544151844011225288800D+00, &
  0.12251482265544137787D+00, &
  0.12251482265544137787D+00, &
  0.12251482265544137787D+00, &
  0.12251482265544137787D+00, &
  0.12251482265544137787D+00, &
  0.12251482265544137787D+00, &
  0.12251482265544137787D+00, &
  0.12251482265544137787D+00, &
  0.12251482265544137787D+00 /)
!
!  Quadrature.
!
  quad = 0.0D+00
  do i = 1, order
    quad = quad + w(i) * func ( x(i), y(i), z(i) )
  end do
!
!  Volume.
!
  volume = pyramid_unit_volume_3d ( )
!
!  Result.
!
  result = quad * volume

  return
end
subroutine pyramid_unit_o27_3d ( func, result )

!*****************************************************************************80
!
!! PYRAMID_UNIT_O27_3D approximates an integral inside the unit pyramid in 3D.
!
!  Discussion:
!
!    A 27 point formula is used.
!
!    The (X,Y,Z) integration region can be represented as:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y) 
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
!    radius of the square diminishes, and when Z reaches 1, the square has 
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 27

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ), dimension ( order ) :: w = (/ &
   0.036374157653908938268D+00, &
   0.05819865224625430123D+00, &
   0.036374157653908938268D+00, &
   0.05819865224625430123D+00, &
   0.09311784359400688197D+00, &
   0.05819865224625430123D+00, &
   0.036374157653908938268D+00, &
   0.05819865224625430123D+00, &
   0.036374157653908938268D+00, &
   0.033853303069413431019D+00, &
   0.054165284911061489631D+00, &
   0.033853303069413431019D+00, &
   0.054165284911061489631D+00, &
   0.08666445585769838341D+00, &
   0.054165284911061489631D+00, &
   0.033853303069413431019D+00, &
   0.054165284911061489631D+00, &
   0.033853303069413431019D+00, &
   0.006933033103838124540D+00, &
   0.011092852966140999264D+00, &
   0.006933033103838124540D+00, &
   0.011092852966140999264D+00, &
   0.017748564745825598822D+00, &
   0.011092852966140999264D+00, &
   0.006933033103838124540D+00, &
   0.011092852966140999264D+00, &
   0.006933033103838124540D+00 /)
  real ( kind = 8 ), dimension ( order ) :: x = (/ &
  -0.7180557413198889387D+00, &
   0.00000000000000000000D+00, &
   0.7180557413198889387D+00, &
  -0.7180557413198889387D+00, &
   0.00000000000000000000D+00, &
   0.7180557413198889387D+00, &
  -0.7180557413198889387D+00, &
   0.00000000000000000000D+00, &
   0.7180557413198889387D+00, &
  -0.50580870785392503961D+00, &
   0.00000000000000000000D+00, &
   0.50580870785392503961D+00, &
  -0.50580870785392503961D+00, &
   0.00000000000000000000D+00, &
   0.50580870785392503961D+00, &
  -0.50580870785392503961D+00, &
   0.00000000000000000000D+00, &
   0.50580870785392503961D+00, &
  -0.22850430565396735360D+00, &
   0.00000000000000000000D+00, &
   0.22850430565396735360D+00, &
  -0.22850430565396735360D+00, &
   0.00000000000000000000D+00, &
   0.22850430565396735360D+00, &
  -0.22850430565396735360D+00, &
   0.00000000000000000000D+00, &
   0.22850430565396735360D+00 /)
  real ( kind = 8 ), dimension ( order ) :: y = (/ &
  -0.7180557413198889387D+00, &
  -0.7180557413198889387D+00, &
  -0.7180557413198889387D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.7180557413198889387D+00, &
   0.7180557413198889387D+00, &
   0.7180557413198889387D+00, &
  -0.50580870785392503961D+00, &
  -0.50580870785392503961D+00, &
  -0.50580870785392503961D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.50580870785392503961D+00, &
   0.50580870785392503961D+00, &
   0.50580870785392503961D+00, &
  -0.22850430565396735360D+00, &
  -0.22850430565396735360D+00, &
  -0.22850430565396735360D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.22850430565396735360D+00, &
   0.22850430565396735360D+00, &
   0.22850430565396735360D+00 /)
  real ( kind = 8 ), dimension ( order ) :: z = (/ &
  0.07299402407314973216D+00, &
  0.07299402407314973216D+00, &
  0.07299402407314973216D+00, &
  0.07299402407314973216D+00, &
  0.07299402407314973216D+00, &
  0.07299402407314973216D+00, &
  0.07299402407314973216D+00, & 
  0.07299402407314973216D+00, &
  0.07299402407314973216D+00, &
  0.34700376603835188472D+00, &
  0.34700376603835188472D+00, &
  0.34700376603835188472D+00, &
  0.34700376603835188472D+00, &
  0.34700376603835188472D+00, &
  0.34700376603835188472D+00, &
  0.34700376603835188472D+00, &
  0.34700376603835188472D+00, &
  0.34700376603835188472D+00, &
  0.70500220988849838312D+00, &
  0.70500220988849838312D+00, &
  0.70500220988849838312D+00, &
  0.70500220988849838312D+00, &
  0.70500220988849838312D+00, &
  0.70500220988849838312D+00, &
  0.70500220988849838312D+00, &
  0.70500220988849838312D+00, &
  0.70500220988849838312D+00 /)
!
!  Quadrature.
!
  quad = 0.0D+00
  do i = 1, order
    quad = quad + w(i) * func ( x(i), y(i), z(i) )
  end do
!
!  Volume.
!
  volume = pyramid_unit_volume_3d ( )
!
!  Result.
!
  result = quad * volume

  return
end
subroutine pyramid_unit_o48_3d ( func, result )

!*****************************************************************************80
!
!! PYRAMID_UNIT_O48_3D approximates an integral inside the unit pyramid in 3D.
!
!  Discussion:
!
!    An 48 point degree 7 formula, Stroud CN:C2:7-1, is used.
!
!    The (X,Y,Z) integration region can be represented as:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y) 
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
!    radius of the square diminishes, and when Z reaches 1, the square has 
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 48

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ), dimension ( order ) :: w = (/ &
  2.01241939442682455D-002, &
  2.01241939442682455D-002, &
  2.01241939442682455D-002, &
  2.01241939442682455D-002, &
  2.60351137043010779D-002, &
  2.60351137043010779D-002, &
  2.60351137043010779D-002, &
  2.60351137043010779D-002, &
  1.24557795239745531D-002, &
  1.24557795239745531D-002, &
  1.24557795239745531D-002, &
  1.24557795239745531D-002, &
  1.87873998794808156D-003, &
  1.87873998794808156D-003, &
  1.87873998794808156D-003, &
  1.87873998794808156D-003, &
  4.32957927807745280D-002, &
  4.32957927807745280D-002, &
  4.32957927807745280D-002, &
  4.32957927807745280D-002, &
  1.97463249834127288D-002, &
  1.97463249834127288D-002, &
  1.97463249834127288D-002, &
  1.97463249834127288D-002, &
  5.60127223523590526D-002, &
  5.60127223523590526D-002, &
  5.60127223523590526D-002, &
  5.60127223523590526D-002, &
  2.55462562927473852D-002, &
  2.55462562927473852D-002, &
  2.55462562927473852D-002, &
  2.55462562927473852D-002, &
  2.67977366291788643D-002, &
  2.67977366291788643D-002, &
  2.67977366291788643D-002, &
  2.67977366291788643D-002, &
  1.22218992265373354D-002, &
  1.22218992265373354D-002, &
  1.22218992265373354D-002, &
  1.22218992265373354D-002, &
  4.04197740453215038D-003, &
  4.04197740453215038D-003, &
  4.04197740453215038D-003, &
  4.04197740453215038D-003, &
  1.84346316995826843D-003, &
  1.84346316995826843D-003, &
  1.84346316995826843D-003, &
  1.84346316995826843D-003 /)
  real ( kind = 8 ), dimension ( order ) :: x = (/ &
  0.88091731624450909D+00, &    
 -0.88091731624450909D+00, &    
   0.0000000000000000D+00, &    
   0.0000000000000000D+00, &    
  0.70491874112648223D+00, &    
 -0.70491874112648223D+00, &    
   0.0000000000000000D+00, &    
   0.0000000000000000D+00, &    
  0.44712732143189760D+00, &    
 -0.44712732143189760D+00, &    
   0.0000000000000000D+00, &    
   0.0000000000000000D+00, &    
  0.18900486065123448D+00, &    
 -0.18900486065123448D+00, &    
   0.0000000000000000D+00, &    
   0.0000000000000000D+00, &    
  0.36209733410322176D+00, &    
 -0.36209733410322176D+00, &    
 -0.36209733410322176D+00, &    
  0.36209733410322176D+00, &    
  0.76688932060387538D+00, &    
 -0.76688932060387538D+00, &    
 -0.76688932060387538D+00, &    
  0.76688932060387538D+00, &    
  0.28975386476618070D+00, &    
 -0.28975386476618070D+00, &    
 -0.28975386476618070D+00, &    
  0.28975386476618070D+00, &    
  0.61367241226233160D+00, &    
 -0.61367241226233160D+00, &    
 -0.61367241226233160D+00, &    
  0.61367241226233160D+00, &    
  0.18378979287798017D+00, &    
 -0.18378979287798017D+00, &    
 -0.18378979287798017D+00, &    
  0.18378979287798017D+00, &    
  0.38925011625173161D+00, &    
 -0.38925011625173161D+00, &    
 -0.38925011625173161D+00, &    
  0.38925011625173161D+00, &    
  7.76896479525748113D-02, &
 -7.76896479525748113D-02, &
 -7.76896479525748113D-02, &
  7.76896479525748113D-02, &
  0.16453962988669860D+00, &    
 -0.16453962988669860D+00, &    
 -0.16453962988669860D+00, &    
  0.16453962988669860D+00 /)  
  real ( kind = 8 ), dimension ( order ) :: y = (/ &
   0.0000000000000000D+00, &     
   0.0000000000000000D+00, &     
  0.88091731624450909D+00, &     
 -0.88091731624450909D+00, &     
   0.0000000000000000D+00, &     
   0.0000000000000000D+00, &     
  0.70491874112648223D+00, &     
 -0.70491874112648223D+00, &    
   0.0000000000000000D+00, &     
   0.0000000000000000D+00, &     
  0.44712732143189760D+00, &     
 -0.44712732143189760D+00, &     
   0.0000000000000000D+00, &     
   0.0000000000000000D+00, &     
  0.18900486065123448D+00, &     
 -0.18900486065123448D+00, &     
  0.36209733410322176D+00, &     
  0.36209733410322176D+00, &     
 -0.36209733410322176D+00, &     
 -0.36209733410322176D+00, &     
  0.76688932060387538D+00, &     
  0.76688932060387538D+00, &     
 -0.76688932060387538D+00, &     
 -0.76688932060387538D+00, &     
  0.28975386476618070D+00, &     
  0.28975386476618070D+00, &     
 -0.28975386476618070D+00, &     
 -0.28975386476618070D+00, &     
  0.61367241226233160D+00, &     
  0.61367241226233160D+00, &     
 -0.61367241226233160D+00, &     
 -0.61367241226233160D+00, &     
  0.18378979287798017D+00, &     
  0.18378979287798017D+00, &     
 -0.18378979287798017D+00, &     
 -0.18378979287798017D+00, &     
  0.38925011625173161D+00, &     
  0.38925011625173161D+00, &     
 -0.38925011625173161D+00, &     
 -0.38925011625173161D+00, &     
  7.76896479525748113D-02, &
  7.76896479525748113D-02, &
 -7.76896479525748113D-02, &
 -7.76896479525748113D-02, &
  0.16453962988669860D+00, &     
  0.16453962988669860D+00, &     
 -0.16453962988669860D+00, &
 -0.16453962988669860D+00 /)
  real ( kind = 8 ), dimension ( order ) :: z = (/ &
  4.85005494469969989D-02, &
  4.85005494469969989D-02, &
  4.85005494469969989D-02, &
  4.85005494469969989D-02, &
  0.23860073755186201D+00, &     
  0.23860073755186201D+00, &     
  0.23860073755186201D+00, &     
  0.23860073755186201D+00, &     
  0.51704729510436798D+00, &     
  0.51704729510436798D+00, &     
  0.51704729510436798D+00, &     
  0.51704729510436798D+00, &     
  0.79585141789677305D+00, &     
  0.79585141789677305D+00, &     
  0.79585141789677305D+00, &     
  0.79585141789677305D+00, &     
  4.85005494469969989D-02, &
  4.85005494469969989D-02, &
  4.85005494469969989D-02, &
  4.85005494469969989D-02, &
  4.85005494469969989D-02, &
  4.85005494469969989D-02, &
  4.85005494469969989D-02, &
  4.85005494469969989D-02, &
  0.23860073755186201D+00, &     
  0.23860073755186201D+00, &     
  0.23860073755186201D+00, &     
  0.23860073755186201D+00, &     
  0.23860073755186201D+00, &     
  0.23860073755186201D+00, &     
  0.23860073755186201D+00, &     
  0.23860073755186201D+00, &     
  0.51704729510436798D+00, &     
  0.51704729510436798D+00, &     
  0.51704729510436798D+00, &     
  0.51704729510436798D+00, &     
  0.51704729510436798D+00, &     
  0.51704729510436798D+00, &     
  0.51704729510436798D+00, &     
  0.51704729510436798D+00, &     
  0.79585141789677305D+00, &     
  0.79585141789677305D+00, &     
  0.79585141789677305D+00, &     
  0.79585141789677305D+00, &     
  0.79585141789677305D+00, &     
  0.79585141789677305D+00, &     
  0.79585141789677305D+00, &
  0.79585141789677305D+00 /)     
!
!  Quadrature.
!
  quad = 0.0D+00
  do i = 1, order
    quad = quad + w(i) * func ( x(i), y(i), z(i) )
  end do
!
!  Volume.
!
  volume = pyramid_unit_volume_3d ( )
!
!  Result.
!
  result = quad * volume

  return
end
function pyramid_unit_monomial_3d ( alpha, beta, gamma )

!*****************************************************************************80
!
!! PYRAMID_UNIT_MONOMIAL_3D: monomial integral in a unit pyramid in 3D.
!
!  Discussion:
!
!    This routine returns the integral of X^ALPHA Y^BETA Z^GAMMA over
!    the unit pyramid.
!
!    The unit pyramid is defined as:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ALPHA, BETA, GAMMA, the exponents of
!    X, Y and Z in the monomial.
!
!    Output, real ( kind = 8 ) PYRAMID_UNIT_MONOMIAL_3D, the volume of 
!    the pyramid.
!
  implicit none

  integer ( kind = 4 ) alpha
  integer ( kind = 4 ) beta
  integer ( kind = 4 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  real ( kind = 8 ) pyramid_unit_monomial_3d
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) value

  value = 0.0D+00

  if ( mod ( alpha, 2 ) == 0 .and. mod ( beta, 2 ) == 0 ) then    

    i_hi = 2 + alpha + beta

    do i = 0, i_hi
      value = value + r8_mop ( i ) * r8_choose ( i_hi, i ) &
      / real ( i + gamma + 1, kind = 8 )
    end do

    value = value &
          * 2.0D+00 / real ( alpha + 1, kind = 8 ) &
          * 2.0D+00 / real ( beta + 1, kind = 8 )

  end if

  pyramid_unit_monomial_3d = value

  return
end
function pyramid_unit_volume_3d ( )

!*****************************************************************************80
!
!! PYRAMID_UNIT_VOLUME_3D: volume of a unit pyramid with square base in 3D.
!
!  Integration region:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!  Discussion:
!
!    The volume of this unit pyramid is 4/3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) PYRAMID_UNIT_VOLUME_3D, the volume of 
!    the pyramid.
!
  implicit none

  real ( kind = 8 ) pyramid_unit_volume_3d

  pyramid_unit_volume_3d = 4.0D+00 / 3.0D+00

  return
end
function pyramid_volume_3d ( r, h )

!*****************************************************************************80
!
!! PYRAMID_VOLUME_3D returns the volume of a pyramid with square base in 3D.
!
!  Integration region:
!
!    - ( H - Z ) * R <= X <= ( H - Z ) * R
!    - ( H - Z ) * R <= Y <= ( H - Z ) * R
!                  0 <= Z <= H.
!
!  Discussion:
!
!    A pyramid with square base can be regarded as the upper half of a
!    3D octahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the "radius" of the pyramid, that is, half the
!    length of one of the sides of the square base.
!
!    Input, real ( kind = 8 ) H, the height of the pyramid.
!
!    Output, real ( kind = 8 ) PYRAMID_VOLUME_3D, the volume of the pyramid.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ) pyramid_volume_3d
  real ( kind = 8 ) r

  pyramid_volume_3d = ( 4.0D+00 / 3.0D+00 ) * h * r * r

  return
end
subroutine qmdpt ( func, n, nsub, result )

!*****************************************************************************80
!
!! QMDPT carries out product midpoint quadrature for the unit cube in ND.
!
!  Integration region:
!
!    -1 <= X(1:N) <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates the function, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the cube.
!
!    Input, integer ( kind = 4 ) NSUB, the number of subdivisions 
!    (in each dimension).
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ix(n)
  logical more
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x(n)

  w = 1.0D+00 / real ( nsub**n, kind = 8 )
  quad = 0.0D+00

  more = .false.
  ihi = nsub**n

  do i = 1, ihi

    call vec_lex_next ( n, nsub, ix, more )

    x(1:n) = real ( 2 * ix(1:n) + 1 - nsub, kind = 8 ) / real ( nsub, kind = 8 )

    quad = quad + w * func ( n, x )

  end do

  volume = 2.0D+00**n
  result = quad * volume

  return
end
function qmult_1d ( func, a, b )

!*****************************************************************************80
!
!! QMULT_1D approximates an integral over an interval in 1D.
!
!  Integration region:
!
!    A <= X <= B.
!
!  Discussion:
!
!    A 16 point 31-st degree Gauss-Legendre formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X), of the form
!      function func ( x )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of integration.
!
!    Output, real ( kind = 8 ) QMULT_1D, the approximate integral of 
!    the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 16

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) qmult_1d
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)

  call legendre_set ( order, xtab, weight )

  quad = 0.0D+00
  do i = 1, order
    x = 0.5D+00 * ( b - a ) * xtab(i) + 0.5D+00 * ( a + b )
    quad = quad + 0.5D+00 * weight(i) * func ( x )
  end do

  volume = b - a
  qmult_1d = quad * volume

  return
end
function qmult_2d ( func, a, b, fup, flo )

!*****************************************************************************80
!
!! QMULT_2D approximates an integral with varying Y dimension in 2D.
!
!  Integration region:
!
!      A <= X <= B
!
!    and
!
!      FLO(X) <= Y <= FHI(X).
!
!  Discussion:
!
!    A 256 point product of two 16 point 31-st degree Gauss-Legendre
!    quadrature formulas is used.
!
!    This routine could easily be modified to use a different
!    order product rule by changing the value of ORDER.
!
!    Another easy change would allow the X and Y directions to
!    use quadrature rules of different orders.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X,Y), of the form
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of X integration.
!
!    Input, external FUP, FLO, the names of the user
!    supplied functions which evaluate the upper and lower
!    limits of the Y integration, of the form
!
!      function fup(x)
!      real ( kind = 8 ) fup
!      real ( kind = 8 ) x
!
!    and
!
!      function flo(x)
!      real ( kind = 8 ) flo
!      real ( kind = 8 ) x
!
!    Output, real ( kind = 8 ) QMULT_2D, the approximate integral of 
!    the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 16

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ), external :: func
  real ( kind = 8 ), external :: flo
  real ( kind = 8 ), external :: fup
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) quad
  real ( kind = 8 ) qmult_2d
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) y

  call legendre_set ( order, xtab, weight )

  quad = 0.0D+00

  do i = 1, order

    w1 = 0.5D+00 * ( b - a ) * weight(i)
    x = 0.5D+00 * ( b - a ) * xtab(i) + 0.5D+00 * ( b + a )
    c = flo ( x )
    d = fup ( x )

    do j = 1, order

      w2 = 0.5D+00 * ( d - c ) * weight(j)
      y = 0.5D+00 * ( d - c ) * xtab(j) + 0.5D+00 * ( d + c )
      quad = quad + w1 * w2 * func ( x, y )

    end do

  end do

  qmult_2d = quad

  return
end
function qmult_3d ( func, a, b, fup1, flo1, fup2, flo2 )

!*****************************************************************************80
!
!! QMULT_3D approximates an integral with varying Y and Z dimension in 3D.
!
!  Integration region:
!
!      A         <= X <= B,
!    and
!      FLO(X)    <= Y <= FHI(X),
!    and
!      FLO2(X,Y) <= Z <= FHI2(X,Y).
!
!  Discussion:
!
!    A 4096 point product of three 16 point 31-st degree Gauss-Legendre
!    quadrature formulas is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of X integration.
!
!    Input, external FUP1, FLO1, the names of the user
!    supplied functions which evaluate the upper and lower
!    limits of the Y integration, of the form
!
!      function fup1(x)
!      real ( kind = 8 ) fup1
!      real ( kind = 8 ) x
!
!    and
!
!      function flo1(x)
!      real ( kind = 8 ) flo1
!      real ( kind = 8 ) x
!
!    Input, external FUP2, FLO2, the names of the user
!    supplied functions which evaluate the upper and lower
!    limits of the Z integration, of the form
!
!      function fup2(x,y)
!      real ( kind = 8 ) fup2
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    and
!
!      function flo2(x,y)
!      real ( kind = 8 ) flo2
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Output, real ( kind = 8 ) QMULT_3D, the approximate integral of 
!    the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 16

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ), external :: func
  real ( kind = 8 ), external :: flo1
  real ( kind = 8 ), external :: flo2
  real ( kind = 8 ), external :: fup1
  real ( kind = 8 ), external :: fup2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) qmult_3d
  real ( kind = 8 ) quad
  real ( kind = 8 ) volume
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  call legendre_set ( order, xtab, weight )

  quad = 0.0D+00

  do i = 1, order

    x = 0.5D+00 * ( b - a ) * xtab(i) + 0.5D+00 * ( b + a )
    w1 = 0.5D+00 * weight(i)
    c = flo1 ( x )
    d = fup1 ( x )

    do j = 1, order

      w2 = 0.5D+00 * ( d - c ) * weight(j)
      y = 0.5D+00 * ( d - c ) * xtab(j) + 0.5D+00 * ( d + c )
      e = flo2 ( x, y )
      f = fup2 ( x, y )

      do k = 1, order

        w3 = 0.5D+00 * ( f - e ) * weight(k)
        z = 0.5D+00 * ( f - e ) * xtab(k) + 0.5D+00 * ( f + e )
        quad = quad + w1 * w2 * w3 * func ( x, y, z )

      end do

    end do

  end do

  volume = b - a
  qmult_3d = quad * volume

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in R8 arithmetic.
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
!    24 March 2008
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
!    Output, real ( kind = 8 ) R8_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0.0D+00

  else if ( mn == 0 ) then

    value = 1.0D+00

  else

    mx = max ( k, n - k )
    value = real ( mx + 1, kind = 8 )

    do i = 2, mn
      value = ( value * real ( mx + i, kind = 8 ) ) / real ( i, kind = 8 )
    end do

  end if

  r8_choose = value

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
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
!    16 January 1999
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
!    Output, real ( kind = 8 ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = 8 )
  end do

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none
!
!  Coefficients for minimax approximation over (12, INF).
!
  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = 8 ), parameter :: eps = 2.22D-16
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 171.624D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.0D+30
  real ( kind = 8 ), parameter :: xminin = 2.23D-308
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00 
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! R8_GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in references by
!    Cody and Hillstrom, and Hillstrom.
!
!    The program uses rational functions that theoretically approximate
!    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
!    approximation for 12 < X is from the Hart reference, while approximations
!    for X < 12.0D+00 are similar to those in the Cody and Hillstrom
!    reference, but are unpublished.
!
!    The accuracy achieved depends on the arithmetic system, the compiler,
!    intrinsic functions, and proper selection of the machine-dependent
!    constants.
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
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Mathematics of Computation,
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
!    Output, real ( kind = 8 ) R8_GAMMA_LOG, the logarithm of the Gamma 
!    function of X.
!
!  Machine-dependent constants:
!
!  BETA   - radix for the floating-point representation.
!
!  MAXEXP - the smallest positive power of BETA that overflows.
!
!  XBIG   - largest argument for which LN(GAMMA(X)) is representable
!           in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!  FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!  Approximate values for some important machines are:
!
!                            BETA      MAXEXP         XBIG
!
!  CRAY-1        (S.P.)        2        8191       9.62D+2461
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)        2         128       4.08D+36
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!  IBM 3033      (D.P.)       16          63       4.29D+73
!
!                           FRTBIG
!
!  CRAY-1        (S.P.)   3.13D+615
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   1.42D+9
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)   2.25D+76
!  IBM 3033      (D.P.)   2.56D+18
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
  real ( kind = 8 ), parameter :: d1 = -5.772156649015328605195174D-01
  real ( kind = 8 ), parameter :: d2 =  4.227843350984671393993777D-01
  real ( kind = 8 ), parameter :: d4 =  1.791759469228055000094023D+00
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: frtbig = 1.42D+09
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
    r8_gamma_log = huge ( r8_gamma_log )
    return
  end if

  if ( x <= epsilon ( x ) ) then

    res = -log ( x )

  else if ( x <= 1.5D+00 ) then

    if ( x < pnt68 ) then
      corr = -log ( x )
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
    xden = -1.0D+00
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

    res = res / x
    corr = log ( x )
    res = res + sqrtpi - 0.5D+00 * corr
    res = res + x * ( corr - 1.0D+00 )

  end if

  r8_gamma_log = res

  return
end
subroutine r8_hyper_2f1 ( a_input, b_input, c_input, x_input, hf )

!*****************************************************************************80
!
!! R8_HYPER_2F1 evaluates the hypergeometric function F(A,B,C,X).
!
!  Discussion:
!
!    A minor bug was corrected.  The HW variable, used in several places as
!    the "old" value of a quantity being iteratively improved, was not
!    being initialized.  JVB, 11 February 2008.
!
!    The original version of this program allowed the input arguments to
!    be modified, although they were restored to their input values before exit.
!    This is unacceptable if the input arguments are allowed to be constants.
!    The code has been modified so that the input arguments are never modified.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2008
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!    The original FORTRAN77 version of this routine is copyrighted by
!    Shanjie Zhang and Jianming Jin.  However, they give permission to
!    incorporate this routine into a user program provided that the copyright
!    is acknowledged.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_INPUT, B_INPUT, C_INPUT, X_INPUT, 
!    the arguments of the function.  The user is allowed to pass these
!    values as constants or variables.
!    C_INPUT must not be equal to a nonpositive integer.
!    X_INPUT < 1.
!
!    Output, real ( kind = 8 ) HF, the value of the function.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a_input
  real ( kind = 8 ) a0
  real ( kind = 8 ) aa
  real ( kind = 8 ) b
  real ( kind = 8 ) b_input
  real ( kind = 8 ) bb
  real ( kind = 8 ) c
  real ( kind = 8 ) c_input
  real ( kind = 8 ) c0
  real ( kind = 8 ) c1
  real ( kind = 8 ), parameter :: el = 0.5772156649015329D+00
  real ( kind = 8 ) eps
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) g0
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) g3
  real ( kind = 8 ) ga
  real ( kind = 8 ) gabc
  real ( kind = 8 ) gam
  real ( kind = 8 ) gb
  real ( kind = 8 ) gbm
  real ( kind = 8 ) gc
  real ( kind = 8 ) gca
  real ( kind = 8 ) gcab
  real ( kind = 8 ) gcb
  real ( kind = 8 ) gm
  real ( kind = 8 ) hf
  real ( kind = 8 ) hw
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical l0
  logical l1
  logical l2
  logical l3
  logical l4
  logical l5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nm
  real ( kind = 8 ) pa
  real ( kind = 8 ) pb
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r0
  real ( kind = 8 ) r1
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) rm
  real ( kind = 8 ) rp
  real ( kind = 8 ) sm
  real ( kind = 8 ) sp
  real ( kind = 8 ) sp0
  real ( kind = 8 ) x
  real ( kind = 8 ) x_input
  real ( kind = 8 ) x1
!
!  Immediately copy the input arguments!
!
  a = a_input
  b = b_input
  c = c_input
  x = x_input

  l0 = ( c == aint ( c ) ) .and. ( c < 0.0D+00 )
  l1 = ( 1.0D+00 - x < 1.0D-15 ) .and. ( c - a - b <= 0.0D+00 )
  l2 = ( a == aint ( a ) ) .and. ( a < 0.0D+00 )
  l3 = ( b == aint ( b ) ) .and. ( b < 0.0D+00 )
  l4 = ( c - a == aint ( c - a ) ) .and. ( c - a <= 0.0D+00 )
  l5 = ( c - b == aint ( c - b ) ) .and. ( c - b <= 0.0D+00 )

  if ( l0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
    write ( *, '(a)' ) '  Integral C < 0.'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    stop
  end if

  if ( l1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    write ( *, '(a)' ) '  1 - X < 0, C - A - B < 0.'
    stop
  end if

  if ( 0.95D+00 < x ) then
    eps = 1.0D-08
  else
    eps = 1.0D-15
  end if

  if ( x == 0.0D+00 .or. a == 0.0D+00 .or. b == 0.0D+00 ) then

    hf = 1.0D+00
    return

  else if ( 1.0D+00 - x == eps .and. 0.0D+00 < c - a - b ) then

    gc = r8_gamma ( c )
    gcab = r8_gamma ( c - a - b )
    gca = r8_gamma ( c - a )
    gcb = r8_gamma ( c - b )
    hf = gc * gcab / ( gca * gcb )
    return

  else if ( 1.0D+00 + x <= eps .and. abs ( c - a + b - 1.0D+00 ) <= eps ) then

    g0 = sqrt ( pi ) * 2.0D+00**( - a )
    g1 = r8_gamma ( c )
    g2 = r8_gamma ( 1.0D+00 + a / 2.0D+00 - b )
    g3 = r8_gamma ( 0.5D+00 + 0.5D+00 * a )
    hf = g0 * g1 / ( g2 * g3 )
    return

  else if ( l2 .or. l3 ) then

    if ( l2 ) then
      nm = int ( abs ( a ) )
    end if

    if ( l3 ) then
      nm = int ( abs ( b ) )
    end if

    hf = 1.0D+00
    r = 1.0D+00

    do k = 1, nm
      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      hf = hf + r
    end do

    return

  else if ( l4 .or. l5 ) then

    if ( l4 ) then
      nm = int ( abs ( c - a ) )
    end if

    if ( l5 ) then
      nm = int ( abs ( c - b ) )
    end if

    hf = 1.0D+00
    r  = 1.0D+00
    do k = 1, nm
      r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      hf = hf + r
    end do
    hf = ( 1.0D+00 - x )**( c - a - b ) * hf
    return

  end if

  aa = a
  bb = b
  x1 = x

  if ( x < 0.0D+00 ) then
    x = x / ( x - 1.0D+00 )
    if ( a < c .and. b < a .and. 0.0D+00 < b ) then
      a = bb
      b = aa
    end if
    b = c - b
  end if

  if ( 0.75D+00 <= x ) then

    gm = 0.0D+00

    if ( abs ( c - a - b - aint ( c - a - b ) ) < 1.0D-15 ) then

      m = int ( c - a - b )
      ga = r8_gamma ( a )
      gb = r8_gamma ( b )
      gc = r8_gamma ( c )
      gam = r8_gamma ( a + m )
      gbm = r8_gamma ( b + m )

      pa = r8_psi ( a )
      pb = r8_psi ( b )

      if ( m /= 0 ) then
        gm = 1.0D+00
      end if

      do j = 1, abs ( m ) - 1
        gm = gm * j
      end do

      rm = 1.0D+00
      do j = 1, abs ( m )
        rm = rm * j
      end do

      f0 = 1.0D+00
      r0 = 1.0D+00
      r1 = 1.0D+00
      sp0 = 0.0D+00
      sp = 0.0D+00

      if ( 0 <= m ) then

        c0 = gm * gc / ( gam * gbm )
        c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )

        do k = 1, m - 1
          r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) &
            + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / real ( k, kind = 8 )
        end do

        f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + ( 1.0D+00 - a ) &
              / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) &
              + 1.0D+00 / ( b + j + k - 1.0D+00 )
          end do

          rp = pa + pb + 2.0D+00 * el + sp + sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + m + k - 1.0D+00 ) * ( b + m + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      else if ( m < 0 ) then

        m = - m
        c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
        c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

        do k = 1, m - 1
          r0 = r0 * ( a - m + k - 1.0D+00 ) * ( b - m + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / real ( k, kind = 8 )
        end do

        f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) &
            / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + 1.0D+00 / real ( j + k, kind = 8 )
          end do

          rp = pa + pb + 2.0D+00 * el + sp - sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      end if

    else

      ga = r8_gamma ( a )
      gb = r8_gamma ( b )
      gc = r8_gamma ( c )
      gca = r8_gamma ( c - a )
      gcb = r8_gamma ( c - b )
      gcab = r8_gamma ( c - a - b )
      gabc = r8_gamma ( a + b - c )
      c0 = gc * gcab / ( gca * gcb )
      c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
      hf = 0.0D+00
      hw = hf
      r0 = c0
      r1 = c1

      do k = 1, 250

        r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
          / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )

        r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
          / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )

        hf = hf + r0 + r1

        if ( abs ( hf - hw ) < abs ( hf ) * eps ) then
          exit
        end if

        hw = hf

      end do

      hf = hf + c0 + c1

    end if

  else

    a0 = 1.0D+00

    if ( a < c .and. c < 2.0D+00 * a .and. b < c .and. c < 2.0D+00 * b ) then

      a0 = ( 1.0D+00 - x )**( c - a - b )
      a = c - a
      b = c - b

    end if

    hf = 1.0D+00
    hw = hf
    r = 1.0D+00

    do k = 1, 250

      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x

      hf = hf + r

      if ( abs ( hf - hw ) <= abs ( hf ) * eps ) then
        exit
      end if

      hw = hf

    end do

    hf = a0 * hf

  end if

  if ( x1 < 0.0D+00 ) then
    x = x1
    c0 = 1.0D+00 / ( 1.0D+00 - x )**aa
    hf = c0 * hf
  end if

  a = aa
  b = bb

  if ( 120 < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Warning!'
    write ( *, '(a)' ) '  A large number of iterations were needed.'
    write ( *, '(a)' ) '  The accuracy of the results should be checked.'
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
function r8_psi ( xx )

!*****************************************************************************80
!
!! R8_PSI evaluates the function Psi(X).
!
!  Discussion:
!
!    This routine evaluates the logarithmic derivative of the
!    Gamma function,
!
!      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
!             = d/dX LN ( GAMMA(X) )
!
!    for real X, where either
!
!      - XMAX1 < X < - XMIN, and X is not a negative integer,
!
!    or
!
!      XMIN < X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Anthony Strecok, Henry Thacher,
!    Chebyshev Approximations for the Psi Function,
!    Mathematics of Computation,
!    Volume 27, Number 121, January 1973, pages 123-127.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XX, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_PSI, the value of the function.
!
  implicit none

  real ( kind = 8 ) aug
  real ( kind = 8 ) den
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nq
  real ( kind = 8 ), dimension ( 9 ) :: p1 = (/ &
   4.5104681245762934160D-03, &
   5.4932855833000385356D+00, &
   3.7646693175929276856D+02, &
   7.9525490849151998065D+03, &
   7.1451595818951933210D+04, &
   3.0655976301987365674D+05, &
   6.3606997788964458797D+05, &
   5.8041312783537569993D+05, &
   1.6585695029761022321D+05 /)
  real ( kind = 8 ), dimension ( 7 ) :: p2 = (/ &
  -2.7103228277757834192D+00, &
  -1.5166271776896121383D+01, &
  -1.9784554148719218667D+01, &
  -8.8100958828312219821D+00, &
  -1.4479614616899842986D+00, &
  -7.3689600332394549911D-02, &
  -6.5135387732718171306D-21 /)
  real ( kind = 8 ), parameter :: piov4 = 0.78539816339744830962D+00
  real ( kind = 8 ), dimension ( 8 ) :: q1 = (/ &
   9.6141654774222358525D+01, &
   2.6287715790581193330D+03, &
   2.9862497022250277920D+04, &
   1.6206566091533671639D+05, &
   4.3487880712768329037D+05, &
   5.4256384537269993733D+05, &
   2.4242185002017985252D+05, &
   6.4155223783576225996D-08 /)
  real ( kind = 8 ), dimension ( 6 ) :: q2 = (/ &
   4.4992760373789365846D+01, &
   2.0240955312679931159D+02, &
   2.4736979003315290057D+02, &
   1.0742543875702278326D+02, &
   1.7463965060678569906D+01, &
   8.8427520398873480342D-01 /)
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) sgn
  real ( kind = 8 ) upper
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: x01 = 187.0D+00
  real ( kind = 8 ), parameter :: x01d = 128.0D+00
  real ( kind = 8 ), parameter :: x02 = 6.9464496836234126266D-04
  real ( kind = 8 ), parameter :: xinf = 1.70D+38
  real ( kind = 8 ), parameter :: xlarge = 2.04D+15
  real ( kind = 8 ), parameter :: xmax1 = 3.60D+16
  real ( kind = 8 ), parameter :: xmin1 = 5.89D-39
  real ( kind = 8 ), parameter :: xsmall = 2.05D-09
  real ( kind = 8 ) xx
  real ( kind = 8 ) z

  x = xx
  w = abs ( x )
  aug = 0.0D+00
!
!  Check for valid arguments, then branch to appropriate algorithm.
!
  if ( xmax1 <= - x .or. w < xmin1 ) then

    if ( 0.0D+00 < x ) then
      r8_psi = - xinf
    else
      r8_psi = xinf
    end if

    return
  end if

  if ( x < 0.5D+00 ) then
!
!  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
!  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
!
    if ( w <= xsmall ) then

      aug = - 1.0D+00 / x
!
!  Argument reduction for cotangent.
!
    else

      if ( x < 0.0D+00 ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - real ( int ( w ), kind = 8 )
      nq = int ( w * 4.0D+00 )
      w = 4.0D+00 * ( w - real ( nq, kind = 8 ) * 0.25D+00 )
!
!  W is now related to the fractional part of 4.0 * X.
!  Adjust argument to correspond to values in the first
!  quadrant and determine the sign.
!
      n = nq / 2

      if ( n + n /= nq ) then
        w = 1.0D+00 - w
      end if

      z = piov4 * w

      if ( mod ( n, 2 ) /= 0 ) then
        sgn = - sgn
      end if
!
!  Determine the final value for  -pi * cotan(pi*x).
!
      n = ( nq + 1 ) / 2
      if ( mod ( n, 2 ) == 0 ) then
!
!  Check for singularity.
!
        if ( z == 0.0D+00 ) then

          if ( 0.0D+00 < x ) then
            r8_psi = - xinf
          else
            r8_psi = xinf
          end if

          return
        end if

        aug = sgn * ( 4.0D+00 / tan ( z ) )

      else

        aug = sgn * ( 4.0D+00 * tan ( z ) )

      end if

    end if

    x = 1.0D+00 - x

  end if
!
!  0.5 <= X <= 3.0.
!
  if ( x <= 3.0D+00 ) then

    den = x
    upper = p1(1) * x
    do i = 1, 7
      den = ( den + q1(i) ) * x
      upper = ( upper + p1(i+1) ) * x
    end do
    den = ( upper + p1(9) ) / ( den + q1(8) )
    x = ( x - x01 / x01d ) - x02
    r8_psi = den * x + aug
    return

  end if
!
!  3.0 < X.
!
  if ( x < xlarge ) then
    w = 1.0D+00 / ( x * x )
    den = w
    upper = p2(1) * w
    do i = 1, 5
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do
    aug = ( upper + p2(7) ) / ( den + q2(6) ) - 0.5D+00 / x + aug
  end if

  r8_psi = aug + log ( x )

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP switches two R8's.
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
subroutine r8_swap3 ( x1, x2, x3 )

!*****************************************************************************80
!
!! R8_SWAP3 swaps three R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X1, X2, X3.
!
!    On output, the values have been shifted so that
!
!      X1 := X3;
!      X2 := X1;
!      X3 := X2;
!
  implicit none

  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  x0 = x1
  x1 = x3
  x3 = x2
  x2 = x0

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
!      seed = 16807 * seed mod ( 2^31 - 1 )
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
subroutine r8ge_det ( n, a_lu, pivot, det )

!*****************************************************************************80
!
!! R8GE_DET: determinant of a matrix factored by R8GE_FA or R8GE_TRF.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!  Modified:
!
!    29 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A_LU(N,N), the LU factors from R8GE_FA 
!    or R8GE_TRF.
!
!    Input, integer ( kind = 4 ) PIVOT(N), as computed by R8GE_FA or R8GE_TRF.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pivot(n)

  det = 1.0D+00

  do i = 1, n
    det = det * a_lu(i,i)
    if ( pivot(i) /= i ) then
      det = - det
    end if
  end do

  return
end
subroutine r8ge_fa ( n, a, pivot, info )

!*****************************************************************************80
!
!! R8GE_FA performs a LINPACK style PLU factorization of a R8GE matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
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
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N,N), the matrix to be factored.
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

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  info = 0

  do k = 1, n-1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k+1, n
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
      write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
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
    do j = k+1, n

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
    write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
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
  character ( len = * ) title

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
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
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
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

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
subroutine r8vec_even_select ( n, xlo, xhi, ival, xval )

!*****************************************************************************80
!
!! R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / real ( N - 1 )
!
!    Unless N = 1, X(1) = XLO and X(N) = XHI.
!
!    If N = 1, then X(1) = 0.5*(XLO+XHI).
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
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) XLO, XHI, the low and high values.
!
!    Input, integer ( kind = 4 ) IVAL, the index of the desired point.
!    IVAL is normally between 1 and N, but may be any integer value.
!
!    Output, real ( kind = 8 ) XVAL, the IVAL-th of N evenly spaced values
!    between XLO and XHI.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ival
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xval

  if ( n == 1 ) then

    xval = 0.5D+00 * ( xlo + xhi )

  else

    xval = ( real ( n - ival,     kind = 8 ) * xlo   &
           + real (     ival - 1, kind = 8 ) * xhi ) &
           / real ( n        - 1, kind = 8 )

  end if

  return
end
subroutine r8vec_mirror_next ( n, a, done )

!*****************************************************************************80
!
!! R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.
!
!  Discussion:
!
!    In normal use, the user would set every element of A to be positive.
!    The routine will take the input value of A, and output a copy in
!    which the signs of one or more entries have been changed.  Repeatedly
!    calling the routine with the output from the previous call will generate
!    every distinct "variation" of A; that is, all possible sign variations.
!
!    When the output variable DONE is TRUE (or equal to 1), then the
!    output value of A_NEW is the last in the series.
!
!    Note that A may have some zero values.  The routine will essentially
!    ignore such entries; more exactly, it will not stupidly assume that -0
!    is a proper "variation" of 0!
!
!    Also, it is possible to call this routine with the signs of A set
!    in any way you like.  The routine will operate properly, but it
!    will nonethess terminate when it reaches the value of A in which
!    every nonzero entry has negative sign.
!
!
!    More efficient algorithms using the Gray code seem to require internal
!    memory in the routine, which is not one of MATLAB's strong points,
!    or the passing back and forth of a "memory array", or the use of
!    global variables, or unnatural demands on the user.  This form of
!    the routine is about as clean as I can make it.
!
!  Example:
!
!      Input         Output
!    ---------    --------------
!    A            A         DONE
!    ---------    --------  ----
!     1  2  3     -1  2  3  false
!    -1  2  3      1 -2  3  false
!     1 -2  3     -1 -2  3  false
!    -1 -2  3      1  2 -3  false
!     1  2 -3     -1  2 -3  false
!    -1  2 -3      1 -2 -3  false
!     1 -2 -3     -1 -2 -3  false
!    -1 -2 -3      1  2  3  true
!
!     1  0  3     -1  0  3  false
!    -1  0  3      1  0 -3  false
!     1  0 -3     -1  0 -3  false
!    -1  0 -3      1  0  3  true
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2004
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
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, real ( kind = 8 ) A(N), a vector of real numbers.  On 
!    output, some signs have been changed.
!
!    Output, logical DONE, is TRUE if the input vector A was the last element
!    in the series (every entry was nonpositive); the output vector is reset
!    so that all entries are nonnegative, but presumably the ride is over!
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) positive
!
!  Seek the first strictly positive entry of A.
!
  positive = 0
  do i = 1, n
    if ( 0.0D+00 < a(i) ) then
      positive = i
      exit
    end if
  end do
!
!  If there is no strictly positive entry of A, there is no successor.
!
  if ( positive == 0 ) then
    a(1:n) = -a(1:n)
    done = .true.
    return
  end if
!
!  Otherwise, negate A up to the positive entry.
!
  a(1:positive) = -a(1:positive)
  done = .false.

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,g16.8)' ) i, a(i)
  end do

  return
end
subroutine rectangle_3d ( func, a, b, result )

!*****************************************************************************80
!
!! RECTANGLE_3D approximates an integral inside a rectangular block in 3D.
!
!  Integration region:
!
!      A(1) <= X <= B(1),
!    and
!      A(2) <= Y <= B(2),
!    and
!      A(3) <= Z <= B(3).
!
!  Discussion:
!
!    An 8 point third degree formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function which
!    evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, real ( kind = 8 ) A(3), B(3), the lower and upper limits
!    for X, Y and Z.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  real ( kind = 8 ) a(3)
  real ( kind = 8 ) b(3)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) sqr3
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  sqr3 = 1.0D+00 / sqrt ( 3.0D+00 )
  w = 1.0D+00 / 8.0D+00

  quad = 0.0D+00

  do i = 1, 2

    x = sqr3 * ( -1 )**i
    x = 0.5D+00 * ( ( 1.0D+00 - x ) * b(1) + ( 1.0D+00 + x ) * a(1) )

    do j = 1, 2

      y = sqr3 * (  -1 )**j
      y = 0.5D+00 * ( ( 1.0D+00 - y ) * b(2) + ( 1.0D+00 + y ) * a(2) )

      do k = 1, 2

        z = sqr3 * ( -1 )**k
        z = 0.5D+00 * ( ( 1.0D+00 - z ) * b(3) + ( 1.0D+00 + z ) * a(3) )

        quad = quad + w * func ( x, y, z )

      end do

    end do

  end do

  volume = ( b(1) - a(1) ) * ( b(2) - a(2) ) * ( b(3) - a(3) )
  result = volume * quad

  return
end
subroutine rectangle_sub_2d ( func, xval, yval, nsub, order, xtab, ytab, &
  weight, result )

!*****************************************************************************80
!
!! RECTANGLE_SUB_2D carries out a composite quadrature over a rectangle in 2D.
!
!  Integration region:
!
!      XVAL(1) <= X <= XVAL(2),
!    and
!      YVAL(1) <= Y <= YVAL(2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the function to be
!    integrated.  The user must declare the name an EXTERNAL
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!      function func ( x, y )
!    which evaluates the function at the point (X,Y).
!
!    Input, real ( kind = 8 ) XVAL(2), the left and right X coordinates.
!
!    Input, real ( kind = 8 ) YVAL(2), the lower and upper Y coordinates.
!
!    Input, integer ( kind = 4 ) NSUB(2).
!    NSUB(1) is the number of subintervals to use in the X direction,
!    and NSUB(2) is the same thing for Y.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b(2)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nsub(2)
  real ( kind = 8 ) quad_sub
  real ( kind = 8 ) result
  real ( kind = 8 ) result_sub
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume_sub
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) xval(2)
  real ( kind = 8 ) y
  real ( kind = 8 ) yhi
  real ( kind = 8 ) ylo
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) yval(2)

  a(1) = xval(1)
  a(2) = yval(1)
  b(1) = xval(2)
  b(2) = yval(2)

  do i = 1, 2
    if ( a(i) == b(i) ) then
      result = 0.0D+00
      return
    end if
  end do

  do i = 1, 2
    if ( nsub(i) < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RECTANGLE_SUB_2D - Fatal error!'
      write ( *, '(a,i8)' ) '  Nonpositive value of NSUB(I) = ', nsub(i)
      write ( *, '(a,i8)' ) '  for index I = ', i
      stop
    end if
  end do
!
!  Break up the X interval into NSUB(1) subintervals.
!
  volume = 0.0D+00
  result = 0.0D+00

  do i = 1, nsub(1)

    call r8vec_even_select ( nsub(1)+1, a(1), b(1), i, xlo )
    call r8vec_even_select ( nsub(1)+1, a(1), b(1), i+1, xhi )
!
!  Break up the Y interval into NSUB(2) subintervals.
!
    do j = 1, nsub(2)

      call r8vec_even_select ( nsub(2)+1, a(2), b(2), j,   ylo )
      call r8vec_even_select ( nsub(2)+1, a(2), b(2), j+1, yhi )

      quad_sub = 0.0D+00
      do k = 1, order

        x = xlo + 0.5D+00 * ( xtab(k) + 1.0D+00 ) * ( xhi - xlo )
        y = ylo + 0.5D+00 * ( ytab(k) + 1.0D+00 ) * ( yhi - ylo )

        quad_sub = quad_sub + weight(k) * func ( x, y ) / 4.0D+00

      end do

      volume_sub = ( xhi - xlo ) * ( yhi - ylo )
      result_sub = quad_sub * volume_sub

      volume = volume + volume_sub
      result = result + result_sub

    end do

  end do

  return
end
subroutine rule_adjust ( a, b, c, d, order, x, w )

!*****************************************************************************80
!
!! RULE_ADJUST maps a quadrature rule from [A,B] to [C,D].
!
!  Discussion:
!
!    Most quadrature rules are defined on a special interval, like
!    [-1,1] or [0,1].  To integrate over an interval, the abscissas
!    and weights must be adjusted.  This can be done on the fly,
!    or by calling this routine.
!
!    If the weight function W(X) is not 1, then the W vector will
!    require further adjustment by the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the definition interval.
!
!    Input, real ( kind = 8 ) C, D, the endpoints of the integration interval.
!
!    Input, integer ( kind = 4 ) ORDER, the number of abscissas and weights.
!
!    Input/output, real ( kind = 8 ) X(ORDER), W(ORDER), the abscissas
!    and weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) w(order)
  real ( kind = 8 ) x(order)

  x(1:order) = ( ( b - x(1:order)     ) * c   &
               + (     x(1:order) - a ) * d ) &
               / ( b              - a )

  w(1:order) = ( ( d - c ) / ( b - a ) ) * w(1:order)

  return
end
subroutine setsim ( n, v )

!*****************************************************************************80
!
!! SETSIM defines a unit simplex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(n,n+1)

  v(1:n,1:n+1) = 0.0D+00

  do i = 1, n
    v(i,i+1) = 1.0D+00
  end do
 
  return
end
subroutine simplex_nd ( func, n, v, result )

!*****************************************************************************80
!
!! SIMPLEX_ND approximates an integral inside a simplex in ND.
!
!  Discussion:
!
!    An N+1 point second degree formula is used.
!
!    The integration region is the simplex bounded by the origin and a 
!    convex combination of N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X) at the N-dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input/output, real ( kind = 8 ) V(N,N+1).  On input, each of the
!    N+1 columns of V contains the N coordinates of one of the
!    "corners" of the simplex in entries 1 through N, with
!    the last column being left free.
!    On output, V has been overwritten in the process of
!    computing the volume of the simplex.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) j
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) simplex_volume_nd
  real ( kind = 8 ) v(n,n+1)
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x(n)

  c = 1.0D+00 / sqrt ( real ( n + 2, kind = 8 ) )
  w = 1.0D+00 / real ( n + 1, kind = 8 )

  do j = 1, n
    x(j) = w * ( 1.0D+00 - c ) * sum ( v(1:n+1,j) )
  end do

  quad = 0.0D+00

  do j = 1, n + 1

    x(1:n) = x(1:n) + c * v(1:n,j)

    quad = quad + w * func ( n, x )

    x(1:n) = x(1:n) - c * v(1:n,j)

  end do

  volume = simplex_volume_nd ( n, v )
  result = quad * volume

  return
end
subroutine simplex_unit_01_nd ( func, n, result )

!*****************************************************************************80
!
!! SIMPLEX_UNIT_01_ND approximates an integral inside the unit simplex in ND.
!
!  Integration region:
!
!      0 <= X(1:N),
!    and
!      sum ( X(1:N) ) <= 1.
!
!  Discussion:
!
!    An 1 point formula of degree 1 is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Axel Grundmann, Michael Moeller,
!    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
!    SIAM Journal on Numerical Analysis,
!    Volume 15, Number 2, April 1978, pages 282-290.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X) at the N-dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.  
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: coef = 1.0D+00
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) simplex_unit_volume_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) x(n)

  quad = 0.0D+00

  x(1:n) = 1.0D+00 / real ( n, kind = 8 )
  quad = quad + coef * func ( n, x )

  volume = simplex_unit_volume_nd ( n )

  result = quad * volume

  return
end
subroutine simplex_unit_03_nd ( func, n, result )

!*****************************************************************************80
!
!! SIMPLEX_UNIT_03_ND approximates an integral inside the unit simplex in ND.
!
!  Integration region:
!
!      0 <= X(1:N),
!    and
!      sum ( X(1:N) ) <= 1.
!
!  Discussion:
!
!    An N+2 point formula of degree 3 is used.  This is Stroud TN:3-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Axel Grundmann, Michael Moeller,
!    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
!    SIAM Journal on Numerical Analysis,
!    Volume 15, Number 2, April 1978, pages 282-290.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X) at the N-dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.  
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) coef
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) simplex_unit_volume_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) x(n)

  quad = 0.0D+00

  x(1:n) = 1.0D+00 / real ( n + 1, kind = 8 )
  coef = -0.25D+00 * real ( ( n + 1 ) * ( n + 1 ), kind = 8 ) &
       / real ( n + 2, kind = 8 )
  quad = quad + coef * func ( n, x )

  a = 1.0D+00 / real ( n + 3, kind = 8 )
  b = 3.0D+00 / real ( n + 3, kind = 8 )

  x(1:n) = a
  coef = 0.25D+00 * real ( ( n + 3 ) * ( n + 3 ), kind = 8 ) &
    / real ( ( n + 1 ) * ( n + 2 ), kind = 8 )
  quad = quad + coef * func ( n, x )

  do i = 1, n

    x(i) = b
    quad = quad + coef * func ( n, x )
    x(i) = a

  end do

  volume = simplex_unit_volume_nd ( n )

  result = quad * volume

  return
end
subroutine simplex_unit_05_nd ( func, n, result )

!*****************************************************************************80
!
!! SIMPLEX_UNIT_05_ND approximates an integral inside the unit simplex in ND.
!
!  Integration region:
!
!      0 <= X(1:N),
!    and
!      sum ( X(1:N) ) <= 1.
!
!  Discussion:
!
!    An N^2 + 3 N + 3 point formula of degree 5 is used.  This is
!    Stroud formula TN:5-1.
!
!    (For N = 2, the number of points is actually only 7, and
!     for N = 3, the number of points is actually only 15.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    A Fifth Degree Integration Formula for the N-Simplex,
!    SIAM Journal on Numerical Analysis,
!    Volume 6, Number 1, March 1969.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X) at the N-dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.  For this 
!    routine, it must be the case that 2 <= N <= 16.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), dimension ( 2 : 16 ) :: coef1 = (/ &
    0.225D+00, &
    0.118518518519D+00, &
    0.0631521898883D+00, &
    0.235714285714D+00, &
    0.791575476992D+00, &
    1.85798728021D+00, &
    3.53666958042D+00, &
    5.90844340844D+00, &
    9.03765432098D+00, &
    12.9758241758D+00, &
    17.7645108738D+00, &
    23.4375030259D+00, &
    30.0224941950D+00, &
    37.5423613501D+00, &
    46.0161454949D+00 /)
  real ( kind = 8 ), dimension ( 2 : 16 ) :: coef21 = (/ &
    0.12593918054483D+00, &
    0.0719370837790D+00, &
    0.0470456145702D+00, &
    0.0333009774677D+00, &
    0.0248633014592D+00, &
    0.0192679696358D+00, &
    0.0153322153879D+00, &
    0.0124316229901D+00, &
    0.0102112988361D+00, &
    0.00845730697460D+00, &
    0.00703433430999D+00, &
    0.00585330520067D+00, &
    0.00485356735291D+00, &
    0.00399261092720D+00, &
    0.00323988713017D+00 /)
  real ( kind = 8 ), dimension ( 2 : 16 ) :: coef22 = (/ &
    0.13239415278851D+00, &
    0.0690682072263D+00, &
    0.0371530185868D+00, &
   -0.0719253160920D+00, &
   -0.264323879461D+00, &
   -0.537926779961D+00, &
   -0.886895605701D+00, &
   -1.30409181465D+00, &
   -1.78227048964D+00, &
   -2.31462336314D+00, &
   -2.89499045158D+00, &
   -3.51790849765D+00, &
   -4.17858310668D+00, &
   -4.87282884913D+00, &
   -5.59699944261D+00 /)
  real ( kind = 8 ), dimension ( 2 : 16 ) :: coef31 = (/ &
    0.0D+00, &
    0.0529100529100D+00, &
    0.0261368740713D+00, &
    0.0499020181331D+00, &
    0.0782233395867D+00, &
    0.109041040862D+00, &
    0.140874828568D+00, &
    0.172735353396D+00, &
    0.203992490408D+00, &
    0.234263814181D+00, &
    0.263332763315D+00, &
    0.291091849264D+00, &
    0.317504208212D+00, &
    0.342577872069D+00, &
    0.366348654344D+00 /)
  real ( kind = 8 ), dimension ( 2 : 16 ) :: coef32 = (/ &
    0.0D+00, &
    0.0D+00, &
    0.0254485903613D+00, &
    0.0165000982690D+00, &
    0.0115218303668D+00, &
    0.00850478779483D+00, &
    0.00655297510968D+00, &
    0.00522372456259D+00, &
    0.00428017828134D+00, &
    0.00358722367033D+00, &
    0.00306362964360D+00, &
    0.00265836687133D+00, &
    0.00233816221525D+00, &
    0.00208061510846D+00, &
    0.00187022027571D+00 /)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) quad
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) result
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) simplex_unit_volume_nd
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) volume
  real ( kind = 8 ) x(n)

  if ( n < 2 .or. 16 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIMPLEX_UNIT_05_ND - Fatal error!'
    write ( *, '(a)' ) '  Input spatial dimension N out of range.'
    write ( *, '(a,i8)' ) '  N = ', n
    result = 0.0D+00
    return
  end if

  quad = 0.0D+00
!
!  S1
!
  x(1:n) = 1.0D+00 / real ( n + 1, kind = 8 )
  quad = quad + coef1(n) * func ( n, x )
!
!  S21
!
  r1 = ( real ( n + 4, kind = 8 ) - sqrt ( 15.0D+00 ) ) &
    / real ( n * n + 8 * n + 1, kind = 8 )
  s1 = 1.0D+00 - real ( n, kind = 8 ) * r1

  x(1:n) = r1

  do i = 1, n + 1

    quad = quad + coef21(n) * func ( n, x )

    if ( 1 < i ) then
      x(i-1) = r1
    end if

    if ( i < n + 1 ) then
      x(i) = s1
    end if

  end do
!
!  S22
!
  r2 = ( real ( n + 4, kind = 8 ) + sqrt ( 15.0D+00 ) ) &
    / real ( n * n + 8 * n + 1, kind = 8 )
  s2 = 1.0D+00 - real ( n, kind = 8 ) * r2

  x(1:n) = r2

  do i = 1, n + 1

    quad = quad + coef22(n) * func ( n, x )

    if ( 1 < i ) then
      x(i-1) = r2
    end if

    if ( i < n + 1 ) then
      x(i) = s2
    end if

  end do
!
!  S31
!
  u1 = ( real ( n + 7, kind = 8 ) + 2.0D+00 * sqrt ( 15.0D+00 ) ) &
    / real ( n * n + 14 * n - 11, kind = 8 )
  v1 = ( real ( 4 * n - 2, kind = 8 ) &
    - real ( n - 1, kind = 8 ) * sqrt ( 15.0D+00 ) ) &
    / real ( n * n + 14 * n - 11, kind = 8 )

  do i = 1, n

    x(1:n) = u1
    x(i) = v1

    do j = i, n

      if ( i < j - 1 ) then
        x(j-1) = u1
      end if

      x(j) = v1

      quad = quad + coef31(n) * func ( n, x )

    end do

  end do
!
!  S32
!
  u2 = ( real ( n + 7, kind = 8 ) - 2.0D+00 * sqrt ( 15.0D+00 ) ) &
    / real ( n * n + 14 * n - 11, kind = 8 )
  v2 = ( real ( 4 * n - 2, kind = 8 ) &
    + real ( n - 1, kind = 8 ) * sqrt ( 15.0D+00 ) ) &
    / real ( n * n + 14 * n - 11, kind = 8 )

  do i = 1, n

    x(1:n) = u2
    x(i) = v2

    do j = i, n

      if ( i < j - 1 ) then
        x(j-1) = u2
      end if

      x(j) = v2

      quad = quad + coef32(n) * func ( n, x )

    end do

  end do

  volume = simplex_unit_volume_nd ( n )

  result = quad * volume

  return
end
subroutine simplex_unit_05_2_nd ( func, n, result )

!*****************************************************************************80
!
!! SIMPLEX_UNIT_05_2_ND approximates an integral inside the unit simplex in ND.
!
!  Integration region:
!
!      0 <= X(1:N),
!    and
!      sum ( X(1:N) ) <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Axel Grundmann, Michael Moeller,
!    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
!    SIAM Journal on Numerical Analysis,
!    Volume 15, Number 2, April 1978, pages 282-290.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X) at the N-dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) coef
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) simplex_unit_volume_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) x(n)

  quad = 0.0D+00
!
!  Group 1
!
  x(1:n) = 1.0D+00 / real ( n + 1, kind = 8 )
  coef = real ( ( n + 1 )**4, kind = 8 ) &
    / real ( 32 * ( n + 2 ) * ( n + 3 ), kind = 8 )
  quad = quad + coef * func ( n, x )
!
!  Group 2
!
  a = 1.0D+00 / real ( n + 3, kind = 8 )
  b = 3.0D+00 / real ( n + 3, kind = 8 )

  x(1:n) = a
  coef = - real ( ( n + 3 )**4, kind = 8 ) &
    / real ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 4 ), kind = 8 )
  quad = quad + coef * func ( n, x )

  do i = 1, n

    x(i) = b
    quad = quad + coef * func ( n, x )
    x(i) = a

  end do
!
!  Group 3
!
  a = 1.0D+00 / real ( n + 5, kind = 8 )
  b = 5.0D+00 / real ( n + 5, kind = 8 )

  x(1:n) = a
  coef = real ( ( n + 5 )**4, kind = 8 ) &
    / real ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 3 ) * ( n + 4 ), kind = 8 )
  quad = quad + coef * func ( n, x )

  do i = 1, n

    x(i) = b
    quad = quad + coef * func ( n, x )
    x(i) = a

  end do
!
!  Group 4
!
  a = 1.0D+00 / real ( n + 5, kind = 8 )
  b = 3.0D+00 / real ( n + 5, kind = 8 )

  coef = real ( ( n + 5 )**4, kind = 8 ) &
    / real ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 3 ) * ( n + 4 ), kind = 8 )

  do i = 1, n

    x(1:n) = a
    x(i) = b
    quad = quad + coef * func ( n, x )

    do j = i+1, n
      x(j) = b
      quad = quad + coef * func ( n, x )
      x(j) = a
    end do

  end do

  volume = simplex_unit_volume_nd ( n )

  result = quad * volume

  return
end
function simplex_unit_volume_nd ( n )

!*****************************************************************************80
!
!! SIMPLEX_UNIT_VOLUME_ND returns the volume of the unit simplex in ND.
!
!  Integration region:
!
!      0 <= X(1:N),
!    and
!      sum ( X(1:N) ) <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) SIMPLEX_UNIT_VOLUME_ND, the volume of the
!    unit simplex.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) i4_factorial
  real ( kind = 8 ) simplex_unit_volume_nd

  simplex_unit_volume_nd = 1.0D+00 / real ( i4_factorial ( n ), kind = 8 )

  return
end
function simplex_volume_nd ( n, v )

!*****************************************************************************80
!
!! SIMPLEX_VOLUME_ND returns the volume of a simplex in ND.
!
!  Integration region:
!
!    The simplex bounded by the origin and a convex combination of N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) V(N,N+1), the coordinates of the
!    vertices.
!
!    Output, real ( kind = 8 ) SIMPLEX_VOLUME_ND, the volume of 
!    the unit simplex.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  real ( kind = 8 ) simplex_unit_volume_nd
  real ( kind = 8 ) simplex_volume_nd
  real ( kind = 8 ) v(n,n+1)
  real ( kind = 8 ) volume
  real ( kind = 8 ) w(n,n)

  do i = 1, n
    w(i,1:n) = v(i,2:n+1) - v(i,1)
  end do

  call r8ge_fa ( n, w, pivot, info )

  call r8ge_det ( n, w, pivot, det )
!
!  Multiply by the volume of the unit simplex, which serves as a
!  conversion factor between a parallelipiped and the simplex.
!
  simplex_volume_nd = abs ( det ) * simplex_unit_volume_nd ( n )

  return
end
function sin_power_int ( a, b, n )

!*****************************************************************************80
!
!! SIN_POWER_INT evaluates the sine power integral.
!
!  Discussion:
!
!    The function is defined by
!
!      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin ( t ))^n dt
!
!    The algorithm uses the following fact:
!
!      Integral sin^n ( t ) = (1/n) * (
!        sin^(n-1)(t) * cos(t) + ( n-1 ) * Integral sin^(n-2) ( t ) dt )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, integer ( kind = 4 ) N, the power of the sine function.
!
!    Output, real ( kind = 8 ) SIN_POWER_INT, the value of the integral.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mlo
  integer ( kind = 4 ) n
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) sin_power_int
  real ( kind = 8 ) value

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIN_POWER_INT - Fatal error!'
    write ( *, '(a)' ) '  Power N < 0.'
    value = 0.0D+00
    stop
  end if

  sa = sin ( a )
  sb = sin ( b )
  ca = cos ( a )
  cb = cos ( b )

  if ( mod ( n, 2 ) == 0 ) then
    value = b - a
    mlo = 2
  else
    value = ca - cb
    mlo = 3
  end if

  do m = mlo, n, 2
    value = ( real ( m - 1, kind = 8 ) * value &
              + sa**(m-1) * ca - sb**(m-1) * cb ) &
      / real ( m, kind = 8 )
  end do

  sin_power_int = value

  return
end
subroutine sphere_05_nd ( func, n, center, r, result )

!*****************************************************************************80
!
!! SPHERE_05_ND approximates an integral on the surface of a sphere in ND.
!
!  Integration region:
!
!    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2
!
!  Discussion:
!
!    A 2*N+2^N points 5-th degree formula is used, Stroud number UN:5-2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X), at the N dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) CENTER(N), the center of the sphere.
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) center(n)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iadd
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ix(n)
  logical more
  integer ( kind = 4 ) ncard
  real ( kind = 8 ) quad
  real ( kind = 8 ) r
  real ( kind = 8 ) result
  real ( kind = 8 ) sphere_area_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  x1 = 1.0D+00
  x2 = 1.0D+00 / sqrt ( real ( n, kind = 8 ) )

  w1 = 1.0D+00 / real ( n * ( n + 2 ), kind = 8 )
  w2 = real ( n, kind = 8 ) / real ( ( n + 2 ) * 2**n, kind = 8 )

  x(1:n) = center(1:n)

  quad = 0.0D+00

  do i = 1, n
    x(i) = center(i) + r * x1
    quad = quad + w1 * func ( n, x )
    x(i) = center(i) - r * x1
    quad = quad + w1 * func ( n, x )
    x(i) = center(i)
  end do

  more = .false.
  ihi = 2**n

  x(1:n) = center(1:n) - r * x2

  do i = 1, ihi

    call subset_gray_next ( n, ix, more, ncard, iadd )

    if ( iadd /= 0 ) then
      x(iadd) = center(iadd) - ( x(iadd) - center(iadd) )
    end if

    quad = quad + w2 * func ( n, x )

  end do

  volume = sphere_area_nd ( n, r )
  result = quad * volume

  return
end
subroutine sphere_07_1_nd ( func, n, center, r, result )

!*****************************************************************************80
!
!! SPHERE_07_1_ND approximates an integral on the surface of a sphere in ND.
!
!  Integration region:
!
!    sum ( X(1:N) - CENTER(1:N) )^2 = R * R.
!
!  Discussion:
!
!    A 2^N + 2*N*N point 7th degree formula is used, Stroud number UN:7-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X), at the N dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) CENTER(N), the center of the sphere.
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) center(n)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iadd
  integer ( kind = 4 ) ix(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  logical more
  integer ( kind = 4 ) ncard
  real ( kind = 8 ) quad
  real ( kind = 8 ) r
  real ( kind = 8 ) result
  real ( kind = 8 ) sphere_area_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  x(1:n) = center(1:n)

  w1 = real ( 8 - n, kind = 8 ) / real ( n * ( n + 2 ) * ( n + 4 ), kind = 8 )
  w2 = real ( n**3, kind = 8 ) &
    / real ( 2**n * n * ( n + 2 ) * ( n + 4 ), kind = 8 )
  w3 = 4.0D+00 / real ( n * ( n + 2 ) * ( n + 4 ), kind = 8 )

  x1 = 1.0D+00
  x2 = 1.0D+00 / sqrt ( real ( n, kind = 8 ) )
  x3 = 1.0D+00 / sqrt ( 2.0D+00 )

  quad = 0.0D+00
!
!  First term.
!
  do i = 1, n
    x(i) = center(i) + r * x1
    quad = quad + w1 * func ( n, x )
    x(i) = center(i) - r * x1
    quad = quad + w1 * func ( n, x )
    x(i) = center(i)
  end do
!
!  Second term.
!
  x(1:n) = center(1:n) - r * x2

  more = .false.
  jhi = 2**n

  do j = 1, jhi

    call subset_gray_next ( n, ix, more, ncard, iadd )

    if ( iadd /= 0 ) then
      x(iadd) = center(iadd) - ( x(iadd) - center(iadd) )
    end if

    quad = quad + w2 * func ( n, x )

  end do
!
!  Third term.
!
  x(1:n) = center(1:n)

  do i = 1, n-1
    do j = i+1, n
      x(i) = center(i) + r * x3
      x(j) = center(j) + r * x3
      quad = quad + w3 * func ( n, x )
      x(i) = center(i) - r * x3
      x(j) = center(j) + r * x3
      quad = quad + w3 * func ( n, x )
      x(i) = center(i) + r * x3
      x(j) = center(j) - r * x3
      quad = quad + w3 * func ( n, x )
      x(i) = center(i) - r * x3
      x(j) = center(j) - r * x3
      quad = quad + w3 * func ( n, x )
      x(i) = center(i)
      x(j) = center(j)
    end do
  end do

  volume = sphere_area_nd ( n, r )
  result = quad * volume

  return
end
function sphere_area_3d ( r )

!*****************************************************************************80
!
!! SPHERE_AREA_3D computes the area of a sphere in 3D.
!
!  Integration region:
!
!    X*X + Y*Y + Z*Z = R * R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Output, real ( kind = 8 ) SPHERE_AREA_3D, the area of the sphere.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) sphere_area_3d

  sphere_area_3d = 4.0D+00 * pi * r * r

  return
end
function sphere_area_nd ( n, r )

!*****************************************************************************80
!
!! SPHERE_AREA_ND computes the area of a sphere in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) = R * R
!
!  Discussion:
!
!    N   Area
!
!    2   2       * PI   * R
!    3   4       * PI   * R^2
!    4   2       * PI^2 * R^3
!    5   (8/3)   * PI^2 * R^4
!    6             PI^3 * R^5
!    7   (16/15) * PI^3 * R^6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Output, real ( kind = 8 ) SPHERE_AREA_ND, the area of the sphere.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  real ( kind = 8 ) sphere_area_nd
  real ( kind = 8 ) sphere_unit_area_nd

  sphere_area_nd = sphere_unit_area_nd ( n ) * r**( n - 1 )

  return
end
subroutine sphere_cap_area_2d ( r, h, area )

!*****************************************************************************80
!
!! SPHERE_CAP_AREA_2D computes the surface area of a spherical cap in 2D.
!
!  Discussion:
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Consider the point on the radius line which is
!    H units from P.  Draw the circle that lies in the plane perpendicular to
!    the radius, and which intersects the sphere.  The circle divides the sphere
!    into two pieces, and the corresponding disk divides the solid sphere into
!    two pieces.  The spherical cap is the part of the solid sphere that
!    includes the point P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) H, the "height" of the spherical cap. 
!
!    Output, real ( kind = 8 ) AREA, the area of the spherical cap.
!
  implicit none

  real ( kind = 8 ) arc_sine
  real ( kind = 8 ) area
  real ( kind = 8 ) h
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) theta

  if ( h <= 0.0D+00 ) then
    area = 0.0D+00
  else if ( 2.0D+00 * r <= h ) then
    area = 2.0D+00 * pi * r
  else

    theta = 2.0D+00 * arc_sine ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r )

    area = r * theta

    if ( r <= h ) then
      area = 2.0D+00 * pi * r - area
    end if

  end if

  return
end
subroutine sphere_cap_area_3d ( r, h, area )

!*****************************************************************************80
!
!! SPHERE_CAP_AREA_3D computes the surface area of a spherical cap in 3D.
!
!  Discussion:
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Consider the point on the radius line which is
!    H units from P.  Draw the circle that lies in the plane perpendicular to
!    the radius, and which intersects the sphere.  The circle divides the sphere
!    into two pieces, and the corresponding disk divides the solid sphere into
!    two pieces.  The spherical cap is the part of the solid sphere that
!    includes the point P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) H, the "height" of the spherical cap. 
!
!    Output, real ( kind = 8 ) AREA, the area of the spherical cap.
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) h
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r

  if ( h <= 0.0D+00 ) then
    area = 0.0D+00
  else if ( 2.0D+00 * r <= h ) then
    area = 4.0D+00 * pi * r * r
  else
    area = 2.0D+00 * pi * r * h
  end if

  return
end
subroutine sphere_cap_area_nd ( dim_num, r, h, area )

!*****************************************************************************80
!
!! SPHERE_CAP_AREA_ND computes the area of a spherical cap in ND.
!
!  Discussion:
!
!    The spherical cap is a portion of the surface of the sphere:
!
!      sum ( X(1:N)^2 ) = R * R
!
!    which is no more than H units from the uppermost point on the sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Ericson, Victor Zinoviev,
!    Codes on Euclidean Spheres,
!    Elsevier, 2001,
!    LC: QA166.7 E75
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) H, the "thickness" of the spherical cap,
!    which is normally between 0 and 2 * R.
!
!    Output, real ( kind = 8 ) AREA, the area of the spherical cap.
!
  implicit none

  real ( kind = 8 ) arc_sine
  real ( kind = 8 ) area
  real ( kind = 8 ) h
  real ( kind = 8 ) haver_sine
  integer ( kind = 4 ) i
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) r
  real ( kind = 8 ) sphere_area_nd
  real ( kind = 8 ) sphere_k
  real ( kind = 8 ) theta
  real ( kind = 8 ) ti
  real ( kind = 8 ) tj
  real ( kind = 8 ) tk

  if ( h <= 0.0D+00 ) then
    area = 0.0D+00
    return
  end if

  if ( 2.0D+00 * r <= h ) then
    area = sphere_area_nd ( dim_num, r )
    return
  end if
!
!  For cases where R < H < 2 * R, work with the complementary region.
!
  haver_sine = sqrt ( ( 2.0D+00 * r - h ) * h )

  theta = arc_sine ( haver_sine / r )

  if ( dim_num < 1 ) then

    area = -1.0D+00

  else if ( dim_num == 1 ) then

    area = 0.0D+00

  else if ( dim_num == 2 ) then

    area = 2.0D+00 * theta * r

  else

    ti = theta

    tj = ti
    ti = 1.0 - cos ( theta )

    do i = 2, dim_num-2
      tk = tj
      tj = ti
      ti = ( real ( i - 1, kind = 8 ) * tk &
        - cos ( theta ) * sin ( theta )**( i - 1 ) ) &
        / real ( i, kind = 8 )
    end do

    area = sphere_k ( dim_num-1 ) * ti * r**( dim_num - 1 )

  end if
!
!  Adjust for cases where R < H < 2R.
!
  if ( r < h ) then
    area = sphere_area_nd ( dim_num, r ) - area
  end if

  return
end
subroutine sphere_cap_volume_2d ( r, h, volume )

!*****************************************************************************80
!
!! SPHERE_CAP_VOLUME_2D computes the volume of a spherical cap in 2D.
!
!  Discussion:
!
!    Draw any radius R of the circle and denote as P the point where the
!    radius intersects the circle.  Now consider the point Q which lies
!    on the radius and which is H units from P.  The line which is
!    perpendicular to the radius R and passes through Q divides the
!    circle into two pieces.  The piece including the point P is the
!    spherical (circular) cap of height (or thickness) H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) H, the "height" of the spherical cap.
!
!    Output, real ( kind = 8 ) VOLUME, the volume (area) of the spherical cap.
!
  implicit none

  real ( kind = 8 ) arc_sine
  real ( kind = 8 ) h
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) theta
  real ( kind = 8 ) volume

  if ( h <= 0.0D+00 ) then

    volume = 0.0D+00

  else if ( 2.0D+00 * r <= h ) then

    volume = pi * r * r

  else

    theta = 2.0D+00 * arc_sine ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r )
    volume = r * r * ( theta - sin ( theta ) ) / 2.0D+00

    if ( r < h ) then
      volume = pi * r * r - volume
    end if

  end if

  return
end
subroutine sphere_cap_volume_3d ( r, h, volume )

!*****************************************************************************80
!
!! SPHERE_CAP_VOLUME_3D computes the volume of a spherical cap in 3D.
!
!  Discussion:
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Consider the point on the radius line which is
!    H units from P.  Draw the circle that lies in the plane perpendicular to
!    the radius, and which intersects the sphere.  The circle divides the sphere
!    into two pieces, and the corresponding disk divides the solid sphere into
!    two pieces.  The part of the solid sphere that includes the point P
!    is the spherical cap of height (or thickness) H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) H, the "height" of the spherical cap.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the spherical cap.
!
  implicit none

  real ( kind = 8 ) h
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) volume

  if ( h <= 0.0D+00 ) then
    volume = 0.0D+00
  else if ( 2.0D+00 * r <= h ) then
    volume = ( 4.0D+00 / 3.0D+00 ) * pi * r * r * r
  else
    volume = ( 1.0D+00 / 3.0D+00 ) * pi * h * h * ( 3.0D+00 * r - h )
  end if

  return
end
subroutine sphere_cap_volume_nd ( dim_num, r, h, volume )

!*****************************************************************************80
!
!! SPHERE_CAP_VOLUME_ND computes the volume of a spherical cap in ND.
!
!  Discussion:
!
!    The spherical cap is a portion of the surface and interior of the sphere:
!
!      sum ( X(1:N)^2 ) <= R * R
!
!    which is no more than H units from some point P on the sphere.
!
!
!    The algorithm proceeds from the observation that the N-dimensional
!    sphere can be parameterized by a quantity RC that runs along the
!    radius from the center to the point P.  The value of RC at the
!    base of the spherical cap is (R-H) and at P it is R.  We intend to
!    use RC as our integration parameeter.
!
!    The volume of the spherical cap is then the integral, as RC goes
!    from (R-H) to R, of the N-1 dimensional volume of the sphere
!    of radius RS, where RC * RC + RS * RS = R * R.
!
!    The volume of the N-1 dimensional sphere of radius RS is simply 
!    some constants times RS**(N-1).
! 
!    After factoring out the constant terms, and writing RC = R * cos ( T ),
!    and RS = R * sin ( T ), and letting 
!      T_MAX = arc_sine ( sqrt ( ( 2.0D+00 * r - h ) * h / r ) ),
!    the "interesting part" of our integral becomes
!
!      constants * R^N * Integral ( T = 0 to T_MAX ) sin^N ( T ) dT
!
!    The integral of sin^N ( T ) dT can be handled by recursion.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) H, the "thickness" of the spherical cap,
!    which is normally between 0 and 2 * R.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the spherical cap.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) arc_sine
  real ( kind = 8 ) factor1
  real ( kind = 8 ) factor2
  real ( kind = 8 ) h
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) r
  real ( kind = 8 ) sin_power_int
  real ( kind = 8 ) sphere_unit_volume_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume2

  if ( h <= 0.0D+00 ) then
    volume = 0.0D+00
    return
  end if

  if ( 2.0D+00 * r <= h ) then
    call sphere_volume_nd ( dim_num, r, volume )
    return
  end if

  if ( dim_num < 1 ) then

    volume = -1.0D+00

  else if ( dim_num == 1 ) then

    volume = h

  else

    factor1 = sphere_unit_volume_nd ( dim_num - 1 )

    angle = arc_sine ( sqrt ( ( 2.0D+00 * r - h ) * h / r ) )

    factor2 = sin_power_int ( 0.0D+00, angle, dim_num )

    volume = factor1 * factor2 * r**dim_num

    if ( r < h ) then
      call sphere_volume_nd ( dim_num, r, volume2 )
      volume = volume2 - volume
    end if

  end if

  return
end
function sphere_k ( n )

!*****************************************************************************80
!
!! SPHERE_K computes a factor useful for spherical computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Ericson, Victor Zinoviev,
!    Codes on Euclidean Spheres,
!    Elsevier, 2001
!    LC: QA166.7 E75
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) SPHERE_K, the factor.
!
  implicit none

  integer ( kind = 4 ) i4_factorial2
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sphere_k

  if ( mod ( n, 2 ) == 0 ) then
    sphere_k = ( 2.0D+00 * pi )**( n / 2 )
  else
    sphere_k = 2.0D+00 * ( 2.0D+00 * pi )**((n-1)/2)
  end if

  sphere_k = sphere_k / real ( i4_factorial2 ( n - 2 ), kind = 8 )

  return
end
subroutine sphere_monomial_int_nd ( n, r, e, integral )

!*****************************************************************************80
!
!! SPHERE_MONOMIAL_INT_ND integrates a monomial on surface of a sphere in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) = R * R.
!
!  Discussion:
!
!    The sphere may have nonunit radius, but it must be centered at 0.
!
!    The monomial is F(X) = X(1)^E(1) * X(2)^E(2) * ... * X(N)^E(N).
!
!    This routine is useful for testing the accuracy of quadrature
!    rules on the sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 June 2002
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
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, integer ( kind = 4 ) E(N), the exponents of X, Y and Z in 
!    the monomial.  Each exponent must be nonnegative.
!
!    Output, real ( kind = 8 ) INTEGRAL, the integral.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) e(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) integral
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_gamma

  if ( any ( e(1:n) < 0 ) ) then
    integral = - huge ( integral )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_MONOMIAL_INT_ND - Fatal error!'
    write ( *, '(a)' ) '  All exponents must be nonnegative.'
    stop
  end if

  if ( all ( e(1:n) == 0 ) ) then

    integral = 2.0D+00 * sqrt ( pi**n ) &
      / r8_gamma ( 0.5D+00 * real ( n, kind = 8 ) )

  else if ( any ( mod ( e(1:n), 2 ) == 1 ) ) then

    integral = 0.0D+00

  else

    integral = 2.0D+00

    do i = 1, n
      integral = integral * r8_gamma ( 0.5D+00 * real ( e(i) + 1, kind = 8 ) )
    end do

    integral = integral &
      / r8_gamma ( 0.5D+00 * ( real ( sum ( e(1:n) + 1 ), kind = 8 ) ) )

  end if

  integral = integral * r**( sum ( e(1:n) ) + 2 )

  return
end
subroutine sphere_shell_03_nd ( func, n, center, r1, r2, result )

!*****************************************************************************80
!
!! SPHERE_SHELL_03_ND approximates an integral inside a spherical shell in ND.
!
!  Integration region:
!
!    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2.
!
!  Discussion:
!
!    An 2*N point 3-rd degree formula is used, Stroud number SN-Shell:3-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F at the N-vector X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) CENTER(N), the center of the spheres.
!
!    Input, real ( kind = 8 ) R1, R2, the inner and outer radiuses that
!    define the spherical shell.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) center(n)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) result
  real ( kind = 8 ) rho
  real ( kind = 8 ) sphere_shell_volume_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x(n)

  if ( r1 == r2 ) then
    result = 0.0D+00
    return
  end if

  rho = r1 / r2

  r = real ( n, kind = 8 ) * ( 1.0D+00 - rho**(n+2) ) &
    / ( real ( n + 2, kind = 8 ) * ( 1.0D+00 - rho**n ) )
  r = sqrt ( r )
  w = 1.0D+00 / real ( 2 * n, kind = 8 )

  x(1:n) = center(1:n)

  quad = 0.0D+00
  do i = 1, n
    x(i) = center(i) + r * r2
    quad = quad + w * func ( n, x )
    x(i) = center(i) - r * r2
    quad = quad + w * func ( n, x )
    x(i) = center(i)
  end do

  volume = sphere_shell_volume_nd ( n, r1, r2 )
  result = quad * volume

  return
end
function sphere_shell_volume_nd ( n, r1, r2 )

!*****************************************************************************80
!
!! SPHERE_SHELL_VOLUME_ND computes the volume of a spherical shell in ND.
!
!  Integration region:
!
!    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) R1, R2, the radiuses of the inner and 
!    outer spheres.
!
!    Output, real ( kind = 8 ) SPHERE_SHELL_VOLUME_ND, the volume of the
!    spherical shell.
!
  implicit none

  real ( kind = 8 ) ball_volume_nd
  integer ( kind = 4 ) n
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) sphere_shell_volume_nd

  sphere_shell_volume_nd = ball_volume_nd ( n, r2 ) &
    - ball_volume_nd ( n, r1 )

  return
end
subroutine sphere_unit_03_nd ( func, n, result )

!*****************************************************************************80
!
!! SPHERE_UNIT_03_ND approximates integral on surface of the unit sphere in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) = 1.
!
!  Discussion:
!
!    A 2*N point 3rd degree formula is used, Stroud number UN:3-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X), at the N dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) sphere_unit_area_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x(n)

  x(1:n) = 0.0D+00

  w = 1.0D+00 / real ( 2 * n, kind = 8 )

  quad = 0.0D+00
  do i = 1, n
    x(i) = 1.0D+00
    quad = quad + w * func ( n, x )
    x(i) = -1.0D+00
    quad = quad + w * func ( n, x )
    x(i) = 0.0D+00
  end do

  volume = sphere_unit_area_nd ( n )
  result = quad * volume

  return
end
subroutine sphere_unit_04_nd ( func, n, result )

!*****************************************************************************80
!
!! SPHERE_UNIT_04_ND approximates integral on surface of the unit sphere in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) = 1.
!
!  Discussion:
!
!    A 2*N*N point 5th degree formula is used, Stroud number UN:5-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X), at the N dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) s
  real ( kind = 8 ) sphere_unit_area_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) x(n)

  x(1:n) = 0.0D+00

  w1 = real ( 4 - n, kind = 8 ) / real ( 2 * n * ( n + 2 ), kind = 8 )

  quad = 0.0D+00

  do i = 1, n
    x(i) = 1.0D+00
    quad = quad + w1 * func ( n, x )
    x(i) = -1.0D+00
    quad = quad + w1 * func ( n, x )
    x(i) = 0.0D+00
  end do

  s = 1.0D+00 / sqrt ( 2.0D+00 )
  w2 = 1.0D+00 / real ( n * ( n + 2 ), kind = 8 )

  do i = 1, n

    x(i) = s

    do j = i+1, n
      x(j) = s
      quad = quad + w2 * func ( n, x )
      x(j) = -s
      quad = quad + w2 * func ( n, x )
      x(j) = 0.0D+00
    end do

    x(i) = -s

    do j = i+1, n
      x(j) = s
      quad = quad + w2 * func ( n, x )
      x(j) = -s
      quad = quad + w2 * func ( n, x )
      x(j) = 0.0D+00
    end do

    x(i) = 0.0D+00

  end do

  volume = sphere_unit_area_nd ( n )
  result = quad * volume

  return
end
subroutine sphere_unit_05_nd ( func, n, result )

!*****************************************************************************80
!
!! SPHERE_UNIT_05_ND approximates integral on surface of the unit sphere in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) = 1.
!
!  Discussion:
!
!    A 2*N+2^N points 5-th degree formula is used, Stroud number UN:5-2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X), at the N dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iadd
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ix(n)
  logical more
  integer ( kind = 4 ) ncard
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) sphere_unit_area_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  x1 = 1.0D+00
  x2 = 1.0D+00 / sqrt ( real ( n, kind = 8 ) )

  w1 = 1.0D+00 / real ( n * ( n + 2 ), kind = 8 )
  w2 = real ( n ) / real ( ( n + 2 ) * 2**n, kind = 8 )

  x(1:n) = 0.0D+00

  quad = 0.0D+00

  do i = 1, n
    x(i) = x1
    quad = quad + w1 * func ( n, x )
    x(i) = -x1
    quad = quad + w1 * func ( n, x )
    x(i) = 0.0D+00
  end do

  more = .false.
  ihi = 2**n

  x(1:n) = -x2

  do i = 1, ihi

    call subset_gray_next ( n, ix, more, ncard, iadd )

    if ( iadd /= 0 ) then
      x(iadd) = -x(iadd)
    end if

    quad = quad + w2 * func ( n, x )

  end do

  volume = sphere_unit_area_nd ( n )
  result = quad * volume

  return
end
subroutine sphere_unit_07_3d ( func, result )

!*****************************************************************************80
!
!! SPHERE_UNIT_07_3D approximates integral on surface of the unit sphere in 3D.
!
!  Integration region:
!
!    X*X + Y*Y + Z*Z = 1.
!
!  Discussion:
!
!    A 32 point 7-th degree formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order1 = 2
  integer ( kind = 4 ), parameter :: order2 = 4
  integer ( kind = 4 ), parameter :: order3 = 4

  real ( kind = 8 ) angle
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) sphere_unit_area_3d
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight1(order1)
  real ( kind = 8 ) weight2(order2)
  real ( kind = 8 ) weight3(order3)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab1(order1)
  real ( kind = 8 ) xtab2(order2)
  real ( kind = 8 ) xtab3(order3)
  real ( kind = 8 ) y
  real ( kind = 8 ) z
!
!  Set XTAB1 and WATE1.
!
  xtab1(1) = -1.0D+00
  xtab1(2) =  1.0D+00
  weight1(1) = 1.0D+00
  weight1(2) = 1.0D+00
!
!  Set XTAB2 and WATE2.
!
  do j = 1, order2
    angle = pi * real ( 2 * j - 1, kind = 8 ) &
      / real ( 2 * order2, kind = 8 )
    xtab2(j) = cos ( angle )
  end do

  weight2(1:order2) = 1.0D+00 / real ( 4 * order2, kind = 8 )
!
!  Set XTAB3 and WATE3.
!
  call legendre_set ( order3, xtab3, weight3 )

  quad = 0.0D+00
  do i = 1, order1
    do j = 1, order2
      do k = 1, order3

        x = xtab1(i) * sqrt ( 1.0D+00 - xtab2(j) * xtab2(j) ) &
                     * sqrt ( 1.0D+00 - xtab3(k) * xtab3(k) )
        y = xtab1(i) * xtab2(j) * sqrt ( 1.0D+00 - xtab3(k) * xtab3(k) )
        z = xtab1(i) * xtab3(k)

        quad = quad + weight1(i) * weight2(j) * weight3(k) * func ( x, y, z )

      end do
    end do
  end do

  volume = sphere_unit_area_3d ( )
  result = quad * volume

  return
end
subroutine sphere_unit_07_1_nd ( func, n, result )

!*****************************************************************************80
!
!! SPHERE_UNIT_07_1_ND approximates integral on surface of unit sphere in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) = 1.
!
!  Discussion:
!
!    A 2^N + 2*N*N point 7th degree formula is used, Stroud number UN:7-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X), at the N dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iadd
  integer ( kind = 4 ) ix(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  logical more
  integer ( kind = 4 ) ncard
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) sphere_unit_area_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  w1 = real ( 8 - n, kind = 8 ) / real ( n * ( n + 2 ) * ( n + 4 ), kind = 8 )
  w2 = real ( n**3, kind = 8 ) &
    / real ( 2**n * n * ( n + 2 ) * ( n + 4 ), kind = 8 )
  w3 = 4.0D+00 / real ( n * ( n + 2 ) * ( n + 4 ), kind = 8 )

  x1 = 1.0D+00
  x2 = 1.0D+00 / sqrt ( real ( n, kind = 8 ) )
  x3 = 1.0D+00 / sqrt ( 2.0D+00 )

  x(1:n) = 0.0D+00

  quad = 0.0D+00
!
!  First term.
!
  do i = 1, n
    x(i) = x1
    quad = quad + w1 * func ( n, x )
    x(i) = -x1
    quad = quad + w1 * func ( n, x )
    x(i) = 0.0D+00
  end do
!
!  Second term.
!
  x(1:n) = -x2

  more = .false.
  jhi = 2**n

  do j = 1, jhi

    call subset_gray_next ( n, ix, more, ncard, iadd )

    if ( iadd /= 0 ) then
      x(iadd) = -x(iadd)
    end if

    quad = quad + w2 * func ( n, x )

  end do
!
!  Third term.
!
  x(1:n) = 0.0D+00

  do i = 1, n-1
    do j = i+1, n
      x(i) = x3
      x(j) = x3
      quad = quad + w3 * func ( n, x )
      x(i) = -x3
      x(j) = x3
      quad = quad + w3 * func ( n, x )
      x(i) = x3
      x(j) = -x3
      quad = quad + w3 * func ( n, x )
      x(i) = -x3
      x(j) = -x3
      quad = quad + w3 * func ( n, x )
      x(i) = 0.0D+00
      x(j) = 0.0D+00
    end do
  end do

  volume = sphere_unit_area_nd ( n )
  result = quad * volume

  return
end
subroutine sphere_unit_07_2_nd ( func, n, result )

!*****************************************************************************80
!
!! SPHERE_UNIT_07_2_ND approximates integral on surface of unit sphere in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) = 1.
!
!  Discussion:
!
!    A 2^N * ( N + 1 ) point 7th degree formula is used, Stroud number UN:7-2.
!
!    Some of the weights in this quadrature formula are negative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X), at the N dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) iadd
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ix(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  logical more
  integer ( kind = 4 ) ncard
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) sphere_unit_area_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3

  x(1:n) = 0.0D+00

  w1 = - real ( n * n, kind = 8 ) / real ( 2**(n+3) * ( n + 2 ), kind = 8 )
  w2 = real ( ( n + 4 ) * ( n + 4 ), kind = 8 ) &
    / real ( 2**(n+3) * n * ( n + 2 ), kind = 8 )
  x1 = 1.0D+00 / sqrt ( real ( n, kind = 8 ) )
  x2 = sqrt ( 5.0D+00 / real ( n + 4, kind = 8 ) )
  x3 = 1.0D+00 / sqrt ( real ( n + 4, kind = 8 ) )

  quad = 0.0D+00

  x(1:n) = - x1

  more = .false.
  jhi = 2**n

  do j = 1, jhi

    call subset_gray_next ( n, ix, more, ncard, iadd )

    if ( iadd /= 0 ) then
      x(iadd) = - x(iadd)
    end if

    quad = quad + w1 * func ( n, x )

  end do

  do i = 1, n

    x(1:n) = - x3

    x(i) = - x2
    more = .false.

    do j = 1, jhi

      call subset_gray_next ( n, ix, more, ncard, iadd )

      if ( iadd /= 0 ) then
        x(iadd) = - x(iadd)
      end if

      quad = quad + w2 * func ( n, x )

    end do

  end do

  volume = sphere_unit_area_nd ( n )
  result = quad * volume

  return
end
subroutine sphere_unit_11_3d ( func, result )

!*****************************************************************************80
!
!! SPHERE_UNIT_11_3D approximates integral on surface of unit sphere in 3D.
!
!  Integration region:
!
!    X*X + Y*Y + Z*Z = 1.
!
!  Discussion:
!
!    A 50 point 11-th degree formula is used, Stroud number U3:11-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AD McLaren,
!    Mathematics of Computation,
!    Volume 17, pages 361-383, 1963.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) sphere_unit_area_3d
  real ( kind = 8 ) volume
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) w4
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  quad = 0.0D+00

  w1 = 9216.0D+00 / 725760.0D+00
  x = 1.0D+00
  y = 0.0D+00
  z = 0.0D+00
  do i = 1, 2
    x = -x
    do j = 1, 3
      call r8_swap3 ( x, y, z )
      quad = quad + w1 * func ( x, y, z )
    end do
  end do

  w2 = 16384.0D+00 / 725760.0D+00
  x = sqrt ( 0.5D+00 )
  y = sqrt ( 0.5D+00 )
  z = 0.0D+00
  do i = 1, 2
    x = -x
    do j = 1, 2
      y = -y
      do k = 1, 3
        call r8_swap3 ( x, y, z )
        quad = quad + w2 * func ( x, y, z )
      end do
    end do
  end do

  w3 = 15309.0D+00 / 725760.0D+00
  x = sqrt ( 1.0D+00 / 3.0D+00 )
  y = sqrt ( 1.0D+00 / 3.0D+00 )
  z = sqrt ( 1.0D+00 / 3.0D+00 )
  do i = 1, 2
    x = -x
    do j = 1, 2
      y = -y
      do k = 1, 2
        z = -z
        quad = quad + w3 * func ( x, y, z )
      end do
    end do
  end do

  w4 = 14641.0D+00 / 725760.0D+00
  x = sqrt ( 1.0D+00 / 11.0D+00 )
  y = sqrt ( 1.0D+00 / 11.0D+00 )
  z = 3.0D+00 * sqrt ( 1.0D+00 / 11.0D+00 )
  do i = 1, 2
    x = -x
    do j = 1, 2
      y = -y
      do k = 1, 2
        z = -z
        do l = 1, 3
          call r8_swap3 ( x, y, z )
          quad = quad + w4 * func ( x, y, z )
        end do
      end do
    end do
  end do

  volume = sphere_unit_area_3d ( )
  result = quad * volume

  return
end
subroutine sphere_unit_11_nd ( func, n, result )

!*****************************************************************************80
!
!! SPHERE_UNIT_11_ND approximates integral on surface of unit sphere in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) = 1
!
!  Discussion:
!
!    An 2^N * ( N^2 + N + 1 ) point formula of degree 5 is used.
!
!    (For N = 3, the number of points is actually only 56, and
!     for N = 4, the number of points is actually only 240.)
!
!    One element of COEF31 was changed from
!      0.0236339091329 to
!      0.0236639091329
!    by Stroud, when going from his paper to his later textbook.
!    This correction was pointed out by David Wright, 16 February 2010.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    A Fifth Degree Integration Formula for the N-Simplex,
!    SIAM Journal on Numerical Analysis,
!    Volume 6, Number 1, March 1969.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X) at the N-dimensional point
!    X, of the form
!      function func ( n, x )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x(n)
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.  For this 
!    routine, it must be the case that 3 <= N <= 16.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ), dimension ( 3 : 16 ) :: coef1 = (/ &
    0.128571428571D+00, &
    0.0518518518518D+00, &
    0.0211979378646D+00, &
    0.281250000000D+00, &
    1.11934731935D+00, &
    2.82751322751D+00, &
    5.68266145619D+00, &
    9.93785824515D+00, &
    15.8196616478D+00, &
    23.5285714285D+00, &
    33.2409299392D+00, &
    45.1113811729D+00, &
    59.2754264177D+00, &
    75.8518518518D+00 /)
  real ( kind = 8 ), dimension ( 3 : 16 ) :: coef21 = (/ &
    0.163795782462D+00, &
    0.0967270533860D+00, &
    0.0638253880175D+00, &
    0.0452340041459D+00, &
    0.0336329118818D+00, &
    0.0261275095270D+00, &
    0.0208331595340D+00, &
    0.0169937111647D+00, &
    0.0141147212492D+00, &
    0.0118949128383D+00, &
    0.0101424250926D+00, &
    0.00873046796644D+00, &
    0.00757257014768D+00, &
    0.00660813369775D+00 /)
  real ( kind = 8 ), dimension ( 3 : 16 ) :: coef22 = (/ &
    0.126680408014D+00, &
    0.0514210947621D+00, &
    0.0213579471658D+00, &
   -0.108726067638D+00, &
   -0.371589499738D+00, &
   -0.786048144448D+00, &
   -1.36034060198D+00, &
   -2.09547695631D+00, &
   -2.98784764467D+00, &
   -4.03107480702D+00, &
   -5.21726499521D+00, &
   -6.53783099707D+00, &
   -7.98401677102D+00, &
   -9.54722261180D+00 /)
  real ( kind = 8 ), dimension ( 3 : 16 ) :: coef31 = (/ &
    0.0D+00, &
    0.0592592592592D+00, &
    0.0236639091329D+00, &
    0.0525940190875D+00, &
    0.0925052768546D+00, &
    0.141316953438D+00, &
    0.196818580052D+00, &
    0.257027634179D+00, &
    0.320299222258D+00, &
    0.385326226441D+00, &
    0.451098131789D+00, &
    0.516849445559D+00, &
    0.582010515746D+00, &
    0.646165210110D+00 /)
  real ( kind = 8 ), dimension ( 3 : 16 ) :: coef32 = (/ &
    0.0D+00, &
    0.0D+00, &
    0.0316246294890D+00, &
    0.0207194729760D+00, &
    0.0144303800811D+00, &
    0.0105348984135D+00, &
    0.00798435122193D+00, &
    0.00623845929545D+00, &
    0.00499896882962D+00, &
    0.00409176297655D+00, &
    0.00341037426698D+00, &
    0.00288710646943D+00, &
    0.00247745182907D+00, &
    0.00215128820597D+00 /)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iadd
  integer ( kind = 4 ) ix(n)
  integer ( kind = 4 ) j
  logical more
  integer ( kind = 4 ) ncard
  real ( kind = 8 ) quad
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) result
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) sphere_unit_area_nd
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) x(n)

  if ( n < 3 .or. 16 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_UNIT_11_ND - Fatal error!'
    write ( *, '(a)' ) '  Input spatial dimension N out of range.'
    write ( *, '(a,i8)' ) '  N = ', n
    result = 0.0D+00
    return
  end if

  quad = 0.0D+00
!
!  S1
!
  x(1:n) = 1.0D+00 / sqrt ( real ( n, kind = 8 ) )

  more = .false.

  do

    call subset_gray_next ( n, ix, more, ncard, iadd )

    if ( iadd /= 0 ) then
      x(iadd) = -x(iadd)
    end if

    quad = quad + coef1(n) * func ( n, x )

    if ( .not. more ) then
      exit
    end if

  end do
!
!  S21
!
  r1 = ( real ( n + 6, kind = 8 ) - 4.0D+00 * sqrt ( 3.0D+00 ) ) &
    / real ( n * n + 12 * n - 12, kind = 8 )
  r1 = sqrt ( r1 )

  s1 = ( real ( 7 * n - 6, kind = 8 ) &
    + real ( 4 * ( n - 1 ), kind = 8 ) * sqrt ( 3.0D+00 ) ) &
    / real ( n * n + 12 * n - 12, kind = 8 )
  s1 = sqrt ( s1 )

  do i = 1, n

    x(1:n) = r1
    x(i) = s1

    more = .false.

    do

      call subset_gray_next ( n, ix, more, ncard, iadd )

      if ( iadd /= 0 ) then
        x(iadd) = -x(iadd)
      end if

      quad = quad + coef21(n) * func ( n, x )

      if ( .not. more ) then
        exit
      end if

    end do

  end do
!
!  S22
!
  r2 = ( real ( n + 6, kind = 8 ) + 4.0D+00 * sqrt ( 3.0D+00 ) ) &
    / real ( n * n + 12 * n - 12, kind = 8 )
  r2 = sqrt ( r2 )

  s2 = ( real ( 7 * n - 6, kind = 8 ) &
    - real ( 4 * ( n - 1 ), kind = 8 ) * sqrt ( 3.0D+00 ) ) &
    / real ( n * n + 12 * n - 12, kind = 8 )
  s2 = sqrt ( s2 )

  do i = 1, n

    x(1:n) = r2
    x(i) = s2

    more = .false.

    do

      call subset_gray_next ( n, ix, more, ncard, iadd )

      if ( iadd /= 0 ) then
        x(iadd) = -x(iadd)
      end if

      quad = quad + coef22(n) * func ( n, x )

      if ( .not. more ) then
        exit
      end if

    end do

  end do
!
!  S31
!
  u1 = ( real ( n + 12, kind = 8 ) + 8.0D+00 * sqrt ( 3.0D+00 ) ) &
    / real ( n * n + 24 * n - 48, kind = 8 )
  u1 = sqrt ( u1 )

  v1 = ( real ( 7 * n - 12, kind = 8 ) &
    - real ( 4 * n - 8, kind = 8 ) * sqrt ( 3.0D+00 ) ) &
    / real ( n * n + 24 * n - 48, kind = 8 )
  v1 = sqrt ( v1 )

  do i = 1, n

    do j = i+1, n

      x(1:n) = u1
      x(i) = v1
      x(j) = v1

      more = .false.

      do

        call subset_gray_next ( n, ix, more, ncard, iadd )

        if ( iadd /= 0 ) then
          x(iadd) = -x(iadd)
        end if

        quad = quad + coef31(n) * func ( n, x )

        if ( .not. more ) then
          exit
        end if

      end do

    end do

  end do
!
!  S32
!
  u2 = ( real ( n + 12, kind = 8 ) - 8.0D+00 * sqrt ( 3.0D+00 ) ) &
    / real ( n * n + 24 * n - 48, kind = 8 )
  u2 = sqrt ( u2 )

  v2 = ( real ( 7 * n - 12, kind = 8 ) &
    + real ( 4 * n - 8, kind = 8 ) * sqrt ( 3.0D+00 ) ) &
    / real ( n * n + 24 * n - 48, kind = 8 )
  v2 = sqrt ( v2 )

  do i = 1, n

    do j = i+1, n

      x(1:n) = u2
      x(i) = v2
      x(j) = v2

      more = .false.

      do

        call subset_gray_next ( n, ix, more, ncard, iadd )

        if ( iadd /= 0 ) then
          x(iadd) = -x(iadd)
        end if

        quad = quad + coef32(n) * func ( n, x )

        if ( .not. more ) then
          exit
        end if

      end do

    end do

  end do

  area = sphere_unit_area_nd ( n )

  result = quad * area / 2.0D+00**n

  return
end
subroutine sphere_unit_14_3d ( func, result )

!*****************************************************************************80
!
!! SPHERE_UNIT_14_3D approximates integral on surface of unit sphere in 3D.
!
!  Integration region:
!
!    X*X + Y*Y + Z*Z = 1.
!
!  Discussion:
!
!    A 72 point 14-th degree formula is used, Stroud number U3:14-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AD McLaren,
!    Mathematics of Computation,
!    Volume 17, pages 361-383, 1963.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) sphere_unit_area_3d
  real ( kind = 8 ) temp
  real ( kind = 8 ) volume
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( 5 ) :: xtab = (/ &
    -0.151108275D+00, 0.315838353D+00, 0.346307112D+00, -0.101808787D+00, &
    -0.409228403D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), save, dimension ( 5 ) :: ytab = (/ &
    0.155240600D+00, 0.257049387D+00, 0.666277790D+00,  0.817386065D+00, &
    0.501547712D+00 /)
  real ( kind = 8 ) z
  real ( kind = 8 ), save, dimension ( 5 ) :: ztab = (/ &
    0.976251323D+00, 0.913330032D+00, 0.660412970D+00,  0.567022920D+00, &
    0.762221757D+00 /)

  quad = 0.0D+00

  w1 = 125.0D+00 / 10080.0D+00
  x = 0.525731112D+00
  y = 0.850650808D+00
  z = 0.0D+00

  do i = 1, 2
    x = -x
    do j = 1, 2
      y = -y
      do k = 1, 3
        call r8_swap3 ( x, y, z )
        quad = quad + w1 * func ( x, y, z )
      end do
    end do
  end do

  w2 = 143.0D+00 / 10080.0D+00

  do i = 1, 5

    x = xtab(i)
    y = ytab(i)
    z = ztab(i)

    do j = 1, 3

      temp = x
      x = z
      z = -y
      y = -temp

      do k = 1, 3
        call r8_swap3 ( x, y, z )
        quad = quad + w2 * func ( x, y, z )
      end do

      y = -y
      z = -z
      quad = quad + w2 * func ( x, y, z )

    end do

  end do

  volume = sphere_unit_area_3d ( )
  result = quad * volume

  return
end
subroutine sphere_unit_15_3d ( func, result )

!*****************************************************************************80
!
!! SPHERE_UNIT_15_3D approximates integral on surface of unit sphere in 3D.
!
!  Integration region:
!
!    X*X + Y*Y + Z*Z = 1.
!
!  Discussion:
!
!    A 128 point 15-th degree spherical product Gauss formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function which evaluates F(X,Y,Z), of the form
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 8

  real ( kind = 8 ) angle
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) sphere_unit_area_3d
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  call legendre_set ( order, xtab, weight )

  weight(1:order) = weight(1:order) / 32.0D+00

  quad = 0.0D+00

  do j = 1, order

    do k = 1, 16

      angle = real ( k, kind = 8 ) * pi / 8.0D+00
      x = sqrt ( 1.0D+00 - xtab(j)**2 ) * cos ( angle )
      y = sqrt ( 1.0D+00 - xtab(j)**2 ) * sin ( angle )
      z = xtab(j)

      quad = quad + weight(j) * func ( x, y, z )

    end do
  end do

  volume = sphere_unit_area_3d ( )
  result = quad * volume

  return
end
function sphere_unit_area_3d ( )

!*****************************************************************************80
!
!! SPHERE_UNIT_AREA_3D computes the surface area of the unit sphere in 3D.
!
!  Integration region:
!
!    X*X + Y*Y + Z*Z = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) SPHERE_UNIT_AREA_3D, the area of the sphere.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sphere_unit_area_3d

  sphere_unit_area_3d = 4.0D+00 * pi

  return
end
function sphere_unit_area_nd ( n )

!*****************************************************************************80
!
!! SPHERE_UNIT_AREA_ND computes the surface area of the unit sphere in ND.
!
!  Integration region:
!
!    sum ( ( X(1:N) - CENTER(1:N) )^2 ) = R * R.
!
!  Discussion:
!
!    N   Area
!
!    2   2       * PI
!    3   4       * PI
!    4   2       * PI^2
!    5   (8/3)   * PI^2
!    6             PI^3
!    7   (16/15) * PI^3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) SPHERE_UNIT_AREA_ND, the area of the sphere.
!
  implicit none

  real ( kind = 8 ) area
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sphere_unit_area_nd

  if ( mod ( n, 2 ) == 0 ) then
    m = n / 2
    area = 2.0D+00 * ( pi )**m
    do i = 1, m-1
      area = area / real ( i, kind = 8 )
    end do
  else
    m = ( n - 1 ) / 2
    area = 2.0D+00**n * ( pi )**m
    do i = m+1, 2*m
      area = area / real ( i, kind = 8 )
    end do
  end if

  sphere_unit_area_nd = area

  return
end
subroutine sphere_unit_area_values ( n_data, n, area )

!*****************************************************************************80
!
!! SPHERE_UNIT_AREA_VALUES returns some areas of the unit sphere in ND.
!
!  Discussion:
!
!    The formula for the surface area of the unit sphere in N dimensions is:
!
!      Sphere_Unit_Area ( N ) = 2 * PI^(N/2) / Gamma ( N / 2 )
!
!    Some values of the function include:
!
!       N   Area
!
!       2    2        * PI
!       3  ( 4 /    ) * PI
!       4  ( 2 /   1) * PI^2
!       5  ( 8 /   3) * PI^2
!       6  ( 1 /   1) * PI^3
!       7  (16 /  15) * PI^3
!       8  ( 1 /   3) * PI^4
!       9  (32 / 105) * PI^4
!      10  ( 1 /  12) * PI^5
!
!    For the unit sphere, Area(N) = N * Volume(N)
!
!    In Mathematica, the function can be evaluated by:
!
!      2 * Pi^(n/2) / Gamma[n/2]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and
!    N_DATA is set to the index of the test data.  On each subsequent
!    call, N_DATA is incremented and that test data is returned.  When
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, real ( kind = 8 ) AREA, the area of the unit sphere
!    in that dimension.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) area
  real ( kind = 8 ), save, dimension ( n_max ) :: area_vec = (/ &
    0.2000000000000000D+01, &
    0.6283185307179586D+01, &
    0.1256637061435917D+02, &
    0.1973920880217872D+02, &
    0.2631894506957162D+02, &
    0.3100627668029982D+02, &
    0.3307336179231981D+02, &
    0.3246969701133415D+02, &
    0.2968658012464836D+02, &
    0.2550164039877345D+02, &
    0.2072514267328890D+02, &
    0.1602315322625507D+02, &
    0.1183817381218268D+02, &
    0.8389703410491089D+01, &
    0.5721649212349567D+01, &
    0.3765290085742291D+01, &
    0.2396678817591364D+01, &
    0.1478625959000308D+01, &
    0.8858104195716824D+00, &
    0.5161378278002812D+00 /)
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     1, &
     2, &
     3, &
     4, &
     5, &
     6, &
     7, &
     8, &
     9, &
    10, &
    11, &
    12, &
    13, &
    14, &
    15, &
    16, &
    17, &
    18, &
    19, &
    20 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    area = 0.0D+00
  else
    n = n_vec(n_data)
    area = area_vec(n_data)
  end if

  return
end
function sphere_unit_monomial_nd ( n, p )

!*****************************************************************************80
!
!! SPHERE_UNIT_MONOMIAL_ND integrate monomial on surface of unit sphere in ND.
!
!  Integration region:
!
!    sum ( X(1:N)^2 ) == 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gerald Folland,
!    How to Integrate a Polynomial Over a Sphere,
!    American Mathematical Monthly,
!    Volume 108, May 2001, pages 446-448.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, integer ( kind = 4 ) P(N), the exponents of X(1) through X(N) 
!    in the monomial.  The exponents P(N) must be nonnegative.
!
!    Output, real ( kind = 8 ) SPHERE_UNIT_MONOMIAL_ND, the integral of
!    X1**P(1)*X2**P(2)*...*XN**P(N) over the unit sphere.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) sphere_unit_monomial_nd
  real ( kind = 8 ) temp

  if ( any ( mod ( p(1:n), 2 ) == 1 ) ) then
    sphere_unit_monomial_nd = 0.0D+00
    return
  end if

  temp = 0.0D+00
  arg2 = 0.0D+00

  do i = 1, n
    arg1 = real ( p(i) + 1, kind = 8 ) / 2.0D+00
    temp = temp + r8_gamma_log ( arg1 )
    arg2 = arg2 + arg1
  end do
  temp = temp - r8_gamma_log ( arg2 )
  
  sphere_unit_monomial_nd = 2.0D+00 * exp ( temp )

  return
end
function sphere_unit_volume_nd ( dim_num )

!*****************************************************************************80
!
!! SPHERE_UNIT_VOLUME_ND computes the volume of a unit sphere in ND.
!
!  Discussion:
!
!    The unit sphere in ND satisfies:
!
!      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
!
!    Results for the first few values of DIM_NUM are:
!
!     DIM_NUM  Volume
!
!     1    2
!     2    1        * PI
!     3  ( 4 /   3) * PI
!     4  ( 1 /   2) * PI^2
!     5  ( 8 /  15) * PI^2
!     6  ( 1 /   6) * PI^3
!     7  (16 / 105) * PI^3
!     8  ( 1 /  24) * PI^4
!     9  (32 / 945) * PI^4
!    10  ( 1 / 120) * PI^5
!
!    For the unit sphere, Volume(DIM_NUM) = 2 * PI * Volume(DIM_NUM-2)/ DIM_NUM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Output, real ( kind = 8 ) SPHERE_UNIT_VOLUME_ND, the volume of the sphere.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) sphere_unit_volume_nd
  real ( kind = 8 ) volume

  if ( mod ( dim_num, 2 ) == 0 ) then
    m = dim_num / 2
    volume = ( pi )**m
    do i = 1, m
      volume = volume / real ( i, kind = 8 )
    end do
  else
    m = ( dim_num - 1 ) / 2
    volume = ( pi )**m * 2.0D+00**dim_num
    do i = m+1, 2*m+1
      volume = volume / real ( i, kind = 8 )
    end do
  end if

  sphere_unit_volume_nd = volume

  return
end
subroutine sphere_unit_volume_values ( n_data, n, volume )

!*****************************************************************************80
!
!! SPHERE_UNIT_VOLUME_VALUES returns some volumes of the unit sphere in ND.
!
!  Discussion:
!
!    The formula for the volume of the unit sphere in N dimensions is
!
!      Volume(N) = 2 * PI**(N/2) / ( N * Gamma ( N / 2 ) )
!
!    This function satisfies the relationships:
!
!      Volume(N) = 2 * PI * Volume(N-2) / N
!      Volume(N) = Area(N) / N
!
!    Some values of the function include:
!
!       N  Volume
!
!       1    1
!       2    1        * PI
!       3  ( 4 /   3) * PI
!       4  ( 1 /   2) * PI^2
!       5  ( 8 /  15) * PI^2
!       6  ( 1 /   6) * PI^3
!       7  (16 / 105) * PI^3
!       8  ( 1 /  24) * PI^4
!       9  (32 / 945) * PI^4
!      10  ( 1 / 120) * PI^5
!
!    In Mathematica, the function can be evaluated by:
!
!      2 * Pi^(n/2) / ( n * Gamma[n/2] )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and
!    N_DATA is set to the index of the test data.  On each subsequent
!    call, N_DATA is incremented and that test data is returned.  When
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the unit
!    sphere in that dimension.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     1,  2, &
     3,  4, &
     5,  6, &
     7,  8, &
     9, 10, &
    11, 12, &
    13, 14, &
    15, 16, &
    17, 18, &
    19, 20 /)
  real ( kind = 8 ) volume
  real ( kind = 8 ), save, dimension ( n_max ) :: volume_vec = (/ &
    0.2000000000000000D+01, &
    0.3141592653589793D+01, &
    0.4188790204786391D+01, &
    0.4934802200544679D+01, &
    0.5263789013914325D+01, &
    0.5167712780049970D+01, &
    0.4724765970331401D+01, &
    0.4058712126416768D+01, &
    0.3298508902738707D+01, &
    0.2550164039877345D+01, &
    0.1884103879389900D+01, &
    0.1335262768854589D+01, &
    0.9106287547832831D+00, &
    0.5992645293207921D+00, &
    0.3814432808233045D+00, &
    0.2353306303588932D+00, &
    0.1409811069171390D+00, &
    0.8214588661112823D-01, &
    0.4662160103008855D-01, &
    0.2580689139001406D-01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    volume = 0.0D+00
  else
    n = n_vec(n_data)
    volume = volume_vec(n_data)
  end if

  return
end
subroutine sphere_volume_2d ( r, volume )

!*****************************************************************************80
!
!! SPHERE_VOLUME_2D computes the volume of an implicit sphere in 2D.
!
!  Discussion:
!
!    An implicit sphere in 2D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - CENTER(1:DIM_NUM) )^2 ) = R * R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the sphere.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) volume

  volume = pi * r * r

  return
end
subroutine sphere_volume_3d ( r, volume )

!*****************************************************************************80
!
!! SPHERE_VOLUME_3D computes the volume of an implicit sphere in 3D.
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - CENTER(1:DIM_NUM) )^2 ) = R * R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the sphere.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) volume

  volume = ( 4.0D+00 / 3.0D+00 ) * pi * r * r * r

  return
end
subroutine sphere_volume_nd ( dim_num, r, volume )

!*****************************************************************************80
!
!! SPHERE_VOLUME_ND computes the volume of an implicit sphere in ND.
!
!  Discussion:
!
!    An implicit sphere in ND satisfies the equation:
!
!      sum ( ( X(1:N) - CENTER(1:N) )^2 ) = R * R
!
!    where R is the radius and CENTER is the center.
!
!    Results for the first few values of N are:
!
!    DIM_NUM  Volume
!    -     -----------------------
!    2                PI   * R^2
!    3     (4/3)    * PI   * R^3
!    4     (1/2)    * PI^2 * R^4
!    5     (8/15)   * PI^2 * R^5
!    6     (1/6)    * PI^3 * R^6
!    7     (16/105) * PI^3 * R^7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the sphere.
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) r
  real ( kind = 8 ) sphere_unit_volume_nd
  real ( kind = 8 ) volume

  volume = r**dim_num * sphere_unit_volume_nd ( dim_num )

  return
end
subroutine square_sum ( func, center, r, order, xtab, ytab, weight, result )

!*****************************************************************************80
!
!! SQUARE_SUM carries out a quadrature rule over a square.
!
!  Integration region:
!
!      abs ( X - CENTER(1) ) <= R 
!    and
!      abs ( Y - CENTER(2) ) <= R
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the function to be
!    integrated.  The user must declare the name an EXTERNAL
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!      function func(x,y)
!    which evaluates the function at the point (X,Y).
!
!    Input, real ( kind = 8 ) CENTER(2), the center of the square.
!
!    Input, real ( kind = 8 ) R, the radius of the square.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas of 
!    the rule.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) order

  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) r
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) y
  real ( kind = 8 ) ytab(order)

  quad = 0.0D+00
  do i = 1, order
    x = center(1) + r * xtab(i)
    y = center(2) + r * ytab(i)
    quad = quad + 0.25D+00 * weight(i) * func ( x, y )
  end do

  volume = 4.0D+00 * r * r
  result = quad * volume

  return
end
subroutine square_unit_set ( rule, order, xtab, ytab, weight )

!*****************************************************************************80
!
!! SQUARE_UNIT_SET sets quadrature weights and abscissas in the unit square.
!
!  Discussion;
!
!    To get the value of ORDER associated with a given rule, 
!    call SQUARE_UNIT_SIZE first.
!
!  Integration region:
!
!      -1 <= X <= 1,
!    and
!      -1 <= Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gilbert Strang, George Fix,
!    An Analysis of the Finite Element Method,
!    Cambridge, 1973,
!    ISBN: 096140888X,
!    LC: TA335.S77.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the rule number.
!    1, order 1, degree 1 rule.
!    2, order 4, degree 3, rule.
!    3, order 9, degree 5 rule.
!    4, order 12 degree 7 rule, Stroud number C2:7-1.
!    5, order 13 degree 7 rule, Stroud number C2:7-3.
!    6, order 64 degree 15 product rule.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ), parameter :: order2 = 8

  real ( kind = 8 ) a
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  integer ( kind = 4 ) rule
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) weight2(order2)
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) xtab2(order2)
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) z

  if ( rule == 1 ) then

    weight(1) = 4.0D+00

    xtab(1) = 0.0D+00
    ytab(1) = 0.0D+00

  else if ( rule == 2 ) then

    a = 1.0D+00
    s = 1.0D+00 / sqrt ( 3.0D+00 )

    xtab(1:4) =   (/ -s, +s, -s, +s /)
    ytab(1:4) =   (/ -s, -s, +s, +s /)
    weight(1:4) = (/  a,  a,  a,  a /)

  else if ( rule == 3 ) then

    s = sqrt ( 0.6D+00 )
    z = 0.0D+00
    w1 = 64.0D+00 / 81.0D+00
    w2 = 25.0D+00 / 81.0D+00
    w3 = 40.0D+00 / 81.0D+00

    xtab(1:9) =   (/   z,  -s, +s, -s, +s,   z, -s, +s,  z /)
    ytab(1:9) =   (/   z,  -s, -s, +s, +s,  -s,  z,  z, +s /)
    weight(1:9) = (/  w1,  w2, w2, w2, w2,  w3, w3, w3, w3 /)

  else if ( rule == 4 ) then

    r = sqrt ( 6.0D+00 / 7.0D+00 )
    c = 3.0D+00 * sqrt ( 583.0D+00 )
    s = sqrt ( ( 114.0D+00 - c ) / 287.0D+00 )
    t = sqrt ( ( 114.0D+00 + c ) / 287.0D+00 )
    w1 = 4.0D+00 * 49.0D+00 / 810.0D+00
    w2 = 4.0D+00 * ( 178981.0D+00 + 923.0D+00 * c ) / 1888920.0D+00
    w3 = 4.0D+00 * ( 178981.0D+00 - 923.0D+00 * c ) / 1888920.0D+00
    z = 0.0D+00

    xtab(1:12) =   (/   r,  z, -r,  z,   s, -s, -s,  s,  t, -t, -t,  t /)
    ytab(1:12) =   (/   z,  r,  z,  -r,  s,  s, -s, -s,  t,  t, -t, -t /)
    weight(1:12) = (/  w1, w1,  w1, w1, w2, w2, w2, w2, w3, w3, w3, w3 /)

  else if ( rule == 5 ) then

    r = sqrt ( 12.0D+00 / 35.0D+00 )
    c = 3.0D+00 * sqrt ( 186.0D+00 )
    s = sqrt ( ( 93.0D+00 + c ) / 155.0D+00 )
    t = sqrt ( ( 93.0D+00 - c ) / 155.0D+00 )
    w1 =  8.0D+00 / 162.0D+00
    w2 = 98.0D+00 / 162.0D+00
    w3 = 31.0D+00 / 162.0D+00
    z = 0.0D+00

    xtab(1:13) =   (/  z,  r, -r,  z,  z,  s,  s, -s, -s,  t,  t, -t, -t /)
    ytab(1:13) =   (/  z,  z,  z,  r, -r,  t, -t,  t, -t,  s, -s,  s, -s /)
    weight(1:13) = (/ w1, w2, w2, w2, w2, w3, w3, w3, w3, w3, w3, w3, w3 /)

  else if ( rule == 6 ) then

    call legendre_set ( order2, xtab2, weight2 )

    k = 0

    do i = 1, order2

      do j = 1, order2

        k = k + 1
        xtab(k) = xtab2(i)
        ytab(k) = xtab2(j)
        weight(k) = weight2(i) * weight2(j)

      end do

    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SQUARE_UNIT_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
    stop

  end if

  return
end
subroutine square_unit_size ( rule, order )

!*****************************************************************************80
!
!! SQUARE_UNIT_SIZE sizes a quadrature rule in the unit square.
!
!  Integration region:
!
!      -1 <= X <= 1,
!    and
!      -1 <= Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gilbert Strang, George Fix,
!    An Analysis of the Finite Element Method,
!    Cambridge, 1973,
!    ISBN: 096140888X,
!    LC: TA335.S77.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the rule number.
!    1, a 1 point 1st degree rule.
!    2, a 4 point 3rd degree rule.
!    3, a 9 point 5th degree rule.
!    4, a 12 point 7-th degree rule, Stroud number C2:7-1.
!    5, a 13 point 7-th degree rule, Stroud number C2:7-3.
!    6, a 64 point 15-th degree product rule.
!
!    Output, integer ( kind = 4 ) ORDER, the order of the rule.
!
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ) rule

  if ( rule == 1 ) then

    order = 1

  else if ( rule == 2 ) then

    order = 4

  else if ( rule == 3 ) then

    order = 9

  else if ( rule == 4 ) then

    order = 12

  else if ( rule == 5 ) then

    order = 13

  else if ( rule == 6 ) then

    order = 64

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SQUARE_UNIT_SIZE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
    stop

  end if

  return
end
subroutine square_unit_sum ( func, order, xtab, ytab, weight, result )

!*****************************************************************************80
!
!! SQUARE_UNIT_SUM carries out a quadrature rule over the unit square.
!
!  Integration region:
!
!      -1 <= X <= 1, 
!    and
!      -1 <= Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the function to be
!    integrated.  The user must declare the name an EXTERNAL
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!      function func ( x, y )
!    which evaluates the function at the point (X,Y).
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas of 
!    the rule.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) ytab(order)

  quad = 0.0D+00
  do i = 1, order
    quad = quad + weight(i) * func ( xtab(i), ytab(i) ) / 4.0D+00
  end do

  volume = 1.0D+00
  result = quad * volume

  return
end
subroutine subset_gray_next ( n, a, more, ncard, iadd )

!*****************************************************************************80
!
!! SUBSET_GRAY_NEXT generates all subsets of a set of order N, one at a time.
!
!  Discussion:
!
!    It generates the subsets one at a time, by adding or subtracting
!    exactly one element on each step.
!
!    This uses a Gray code ordering of the subsets.
!
!    The user should set MORE = FALSE and the value of N before
!    the first call.  On return, the user may examine A which contains
!    the definition of the new subset, and must check MORE, because
!    as soon as it is FALSE on return, all the subsets have been
!    generated and the user probably should cease calling.
!
!    The first set returned is the empty set.
!
!  Example:
!
!    N = 4
!
!    0 0 0 0
!    1 0 0 0
!    1 1 0 0
!    0 1 0 0
!    0 1 1 0
!    1 1 1 0
!    1 0 1 0
!    0 0 1 0
!    0 0 1 1
!    1 0 1 1
!    1 1 1 1
!    0 1 1 1
!    0 1 0 1
!    1 1 0 1
!    1 0 0 1
!    0 0 0 1
!
!  Modified:
!
!    02 May 2003
!
!  Author:
!
!    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the total set from which
!    subsets will be drawn.
!
!    Input/output, integer ( kind = 4 ) A(N).  On each return, the Gray code 
!    for the newly generated subset.  A(I) = 0 if element I is in the subset,
!    1 otherwise.
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  Normally, MORE will be returned TRUE but once
!    all the subsets have been generated, MORE will be
!    reset FALSE on return and you should stop calling the program.
!
!    Input/output, integer ( kind = 4 ) NCARD, the cardinality of the set 
!    returned, which may be any value between 0 (the empty set) and N (the
!    whole set).
!
!    Output, integer ( kind = 4 ) IADD, the element which was added or removed
!    to the previous subset to generate the current one.  Exception:
!    the empty set is returned on the first call, and IADD is set to 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), parameter :: i4_2 = 2
  integer ( kind = 4 ) iadd
  logical more
  integer ( kind = 4 ) ncard
!
!  The first set returned is the empty set.
!
  if ( .not. more ) then

    a(1:n) = 0

    iadd = 0
    ncard = 0
    more = .true.

  else

    iadd = 1

    if ( mod ( ncard, i4_2 ) /= 0 ) then

      do

        iadd = iadd + 1
        if ( a(iadd-1) /= 0 ) then
          exit
        end if

      end do

    end if

    a(iadd) = 1 - a(iadd)
    ncard = ncard + 2 * a(iadd) - 1
!
!  The last set returned is the singleton A(N).
!
    if ( ncard == a(n) ) then
      more = .false.
    end if

  end if

  return
end
subroutine tetra_07 ( func, x, y, z, result )

!*****************************************************************************80
!
!! TETRA_07 approximates an integral inside a tetrahedron in 3D.
!
!  Integration region:
!
!    Points inside a tetrahedron whose four corners are given.
!
!  Discussion:
!
!    A 64 point 7-th degree conical product Gauss formula is used,
!    Stroud number T3:7-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966, pages 42-43,
!    LC: QA299.4G3S7
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of three variables which is to be integrated,
!    of the form:
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, real ( kind = 8 ) X(4), Y(4), Z(4), the coordinates of 
!    the vertices.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 4

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) t
  real ( kind = 8 ) tetra_volume
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) weight1(order)
  real ( kind = 8 ), save, dimension ( order ) :: weight2 = (/ &
    0.1355069134D+00, 0.2034645680D+00, 0.1298475476D+00, 0.0311809709D+00 /)
  real ( kind = 8 ), save, dimension ( order ) :: weight3 = (/ &
    0.1108884156D+00, 0.1434587898D+00, 0.0686338872D+00, 0.0103522407D+00 /)
  real ( kind = 8 ) x(4)
  real ( kind = 8 ) xtab1(order)
  real ( kind = 8 ), save, dimension ( order ) :: xtab2 = (/ &
    0.0571041961D+00, 0.2768430136D+00, 0.5835904324D+00, 0.8602401357D+00 /)
  real ( kind = 8 ), save, dimension ( order ) :: xtab3 = (/ &
    0.0485005495D+00, 0.2386007376D+00, 0.5170472951D+00, 0.7958514179D+00 /)
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(4)
  real ( kind = 8 ) yval
  real ( kind = 8 ) z(4)
  real ( kind = 8 ) zval
!
!  Get the Gauss-Legendre weights and abscissas for [-1,1].
!
  call legendre_set ( order, xtab1, weight1 )
!
!  Adjust the rule for the interval [0,1].
!
  a = -1.0D+00
  b = +1.0D+00

  c =  0.0D+00
  d =  1.0D+00

  call rule_adjust ( a, b, c, d, order, xtab1, weight1 )
!
!  Carry out the quadrature.
!
  quad = 0.0D+00

  do i = 1, order
    do j = 1, order
      do k = 1, order
!
!  Compute the barycentric coordinates of the point in the unit triangle.
!
        t =                                                 xtab3(k)
        u =                        xtab2(j)   * ( 1.0D+00 - xtab3(k) )
        v = xtab1(i) * ( 1.0D+00 - xtab2(j) ) * ( 1.0D+00 - xtab3(k) )
        w = 1.0D+00 - t - u - v
!
!  Compute the corresponding point in the triangle.
!
        xval = t * x(1) + u * x(2) + v * x(3) + w * x(4)
        yval = t * y(1) + u * y(2) + v * y(3) + w * y(4)
        zval = t * z(1) + u * z(2) + v * z(3) + w * z(4)

        quad = quad + 6.0D+00 * weight1(i) * weight2(j) * weight3(k) &
          * func ( xval, yval, zval )

      end do
    end do
  end do

  volume = tetra_volume ( x, y, z )
  result = quad * volume

  return
end
subroutine tetra_sum ( func, x, y, z, order, xtab, ytab, ztab, weight, result )

!*****************************************************************************80
!
!! TETRA_SUM carries out a quadrature rule in a tetrahedron in 3D.
!
!  Integration region:
!
!    A tetrahedron whose vertices are specified.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, name of the function, of the form:
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, real ( kind = 8 ) X(4), Y(4), Z(4), the vertices.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), ZTAB(ORDER), the
!    abscissas.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) tetra_volume
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x(4)
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(4)
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) yval
  real ( kind = 8 ) z(4)
  real ( kind = 8 ) ztab(order)
  real ( kind = 8 ) zval

  quad = 0.0D+00

  do i = 1, order

    xval =             xtab(i)                       * x(1) &
                               + ytab(i)             * x(2) &
                                         + ztab(i)   * x(3) &
         + ( 1.0D+00 - xtab(i) - ytab(i) - ztab(i) ) * x(4)

    yval =             xtab(i)                       * y(1) &
                               + ytab(i)             * y(2) &
                                         + ztab(i)   * y(3) &
         + ( 1.0D+00 - xtab(i) - ytab(i) - ztab(i) ) * y(4)

    zval =             xtab(i)                       * z(1) &
                               + ytab(i)             * z(2) &
                                         + ztab(i)   * z(3) &
         + ( 1.0D+00 - xtab(i) - ytab(i) - ztab(i) ) * z(4)

    quad = quad + weight(i) * func ( xval, yval, zval )

  end do

  volume = tetra_volume ( x, y, z )
  result = quad * volume

  return
end
subroutine tetra_tproduct ( func, order, x, y, z, result )

!*****************************************************************************80
!
!! TETRA_TPRODUCT approximates an integral in a tetrahedron in 3D.
!
!  Discussion:
!
!    Integration is carried out over the points inside an arbitrary
!    tetrahedron whose four vertices are given.
!
!    An ORDER**3 point (2*ORDER-1)-th degree triangular product
!    Gauss-Legendre rule is used.
!
!    With ORDER = 8, this routine is equivalent to the routine TETR15
!    in the reference, page 367.
!
!    Thanks to Joerg Behrens, jbehren@gwdg.de, for numerous suggestions
!    and corrections.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of three variables which is to be integrated,
!    of the form:
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, integer ( kind = 4 ) ORDER, the order of the basic quadrature rules.
!    ORDER should be between 1 and 9.
!
!    Input, real ( kind = 8 ) X(4), Y(4), Z(4), the vertices
!    of the tetrahedron.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) tetra_volume
  real ( kind = 8 ) volume
  real ( kind = 8 ), dimension ( order ) :: weight0
  real ( kind = 8 ), dimension ( order ) :: weight1
  real ( kind = 8 ), dimension ( order ) :: weight2
  real ( kind = 8 ) x(4)
  real ( kind = 8 ), dimension ( order ) :: xtab0
  real ( kind = 8 ), dimension ( order ) :: xtab1
  real ( kind = 8 ), dimension ( order ) :: xtab2
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(4)
  real ( kind = 8 ) yval
  real ( kind = 8 ) z(4)
  real ( kind = 8 ) zval

  if ( order < 1 .or. 9 < order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRA_TPRODUCT - Fatal error!'
    write ( *, '(a)' ) '  The quadrature rule orders must be between 1 and 9.'
    write ( *, '(a,i8)' ) '  The input value was ORDER = ', order
    stop
  end if
!
!  Get the Gauss-Legendre ORDER point rules on [-1,1] for integrating
!    F(X),
!    X * F(X),
!    X * X * F(X).
!
  call legendre_set ( order, xtab0, weight0 )
  call legendre_set_x1 ( order, xtab1, weight1 )
  call legendre_set_x2 ( order, xtab2, weight2 )
!
!  Adjust the rules from [-1,1] to [0,1].
!
  a = -1.0D+00
  b = +1.0D+00
  c =  0.0D+00
  d =  1.0D+00

  call rule_adjust ( a, b, c, d, order, xtab0, weight0 )

  call rule_adjust ( a, b, c, d, order, xtab1, weight1 )

  call rule_adjust ( a, b, c, d, order, xtab2, weight2 )
!
!  For rules with a weight function that is not 1, the weight vectors
!  require further adjustment.
!
  weight1(1:order) = weight1(1:order) / 2.0D+00
  weight2(1:order) = weight2(1:order) / 4.0D+00
!
!  Carry out the quadrature.
!
  quad = 0.0D+00

  do k = 1, order
    do j = 1, order
      do i = 1, order

        xval = x(1) + ( ( ( x(4) - x(3) )   * xtab0(i) &
                        + ( x(3) - x(2) ) ) * xtab1(j) &
                        + ( x(2) - x(1) ) ) * xtab2(k)

        yval = y(1) + ( ( ( y(4) - y(3) )   * xtab0(i) &
                        + ( y(3) - y(2) ) ) * xtab1(j) &
                        + ( y(2) - y(1) ) ) * xtab2(k)

        zval = z(1) + ( ( ( z(4) - z(3) )   * xtab0(i) &
                        + ( z(3) - z(2) ) ) * xtab1(j) &
                        + ( z(2) - z(1) ) ) * xtab2(k)

        quad = quad + 6.0D+00 * weight0(i) * weight1(j) * weight2(k) &
          * func ( xval, yval, zval )

      end do

    end do

  end do
!
!  Compute the volume of the tetrahedron.
!
  volume = tetra_volume ( x, y, z )
  result = quad * volume

  return
end
subroutine tetra_unit_set ( rule, order, xtab, ytab, ztab, weight )

!*****************************************************************************80
!
!! TETRA_UNIT_SET sets quadrature weights and abscissas in the unit tetrahedron.
!
!  Integration region:
!
!      0 <= X,
!    and
!      0 <= Y,
!    and
!      0 <= Z, 
!    and
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hermann Engels,
!    Numerical Quadrature and Cubature,
!    Academic Press, 1980,
!    ISBN: 012238850X,
!    LC: QA299.3E5.
!
!    Patrick Keast,
!    Moderate Degree Tetrahedral Quadrature Formulas,
!    Computer Methods in Applied Mechanics and Engineering,
!    Volume 55, Number 3, May 1986, pages 339-348.
!
!    Olgierd Zienkiewicz,
!    The Finite Element Method,
!    Sixth Edition,
!    Butterworth-Heinemann, 2005,
!    ISBN: 0750663200,
!    LC: TA640.2.Z54
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!     1, order 1, precision 0, Newton Cotes formula #0, Zienkiewicz #1.
!     2, order 4, precision 1, Newton Cotes formula #1.
!     3, order 4, precision 2, Zienkiewicz #2.
!     4, order 10, precision 2, Newton Cotes formula #2
!     5, order 5, precision 3, Zienkiewicz #3.
!     6, order 8, precision 3, Newton Cotes formula #3.
!     7, order 35, precision 4, Newton Cotes formula #4.
!     8, order 11, precision 4, a Keast rule.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), ZTAB(ORDER),
!    the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) rule
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) z
  real ( kind = 8 ) ztab(order)
!
!  Newton Cotes #0.
!
  if ( rule == 1 ) then

    xtab(1) = 1.0D+00 / 4.0D+00
    ytab(1) = 1.0D+00 / 4.0D+00
    ztab(1) = 1.0D+00 / 4.0D+00
    weight(1) = 1.0D+00
!
!  Newton Cotes #1.
!
  else if ( rule == 2 ) then

    a = 1.0D+00
    b = 1.0D+00 / 4.0D+00
    z = 0.0D+00

    xtab(1:4) =   (/ z, a, z, z /)
    ytab(1:4) =   (/ z, z, a, z /)
    ztab(1:4) =   (/ z, z, z, a /)
    weight(1:4) = (/ b, b, b, b /)
!
!  Zienkiewicz #2.
!
  else if ( rule == 3 ) then

    a =  0.5854101966249685D+00
    b =  0.1381966011250105D+00
    c =  0.25D+00

    xtab(1:4) =   (/ a, b, b, b /)
    ytab(1:4) =   (/ b, a, b, b /)
    ztab(1:4) =   (/ b, b, a, b /)
    weight(1:4) = (/ c, c, c, c /)
!
!  Newton Cotes #2.
!
  else if ( rule == 4 ) then

    a =  1.0D+00
    b =  0.5D+00
    c = -1.0D+00 / 20.0D+00
    d =  4.0D+00 / 20.0D+00
    z =  0.0D+00

    xtab(1:10) =   (/ z, a, z, z, b, z, z, b, b, z /)
    ytab(1:10) =   (/ z, z, a, z, z, b, z, b, z, b /)
    ztab(1:10) =   (/ z, z, z, a, z, z, b, z, b, b /)
    weight(1:10) = (/ c, c, c, c, d, d, d, d, d, d /)
!
!  Zienkiewicz #3.
!
  else if ( rule == 5 ) then

    a =  1.0D+00 / 6.0D+00
    b =  0.25D+00
    c =  0.5D+00
    d = -0.8D+00
    e =  0.45D+00

    xtab(1:5) =   (/ b, c, a, a, a /)
    ytab(1:5) =   (/ b, a, c, a, a /)
    ztab(1:5) =   (/ b, a, a, c, a /)
    weight(1:5) = (/ d, e, e, e, e /)
!
!  Newton Cotes #3.
!  (This is actually formally a 20 point rule, but with 12 zero coefficients!)
!
  else if ( rule == 6 ) then

    a = 1.0D+00
    b = 1.0D+00 / 40.0D+00
    c = 1.0D+00 /  3.0D+00
    d = 9.0D+00 / 40.0D+00
    z = 0.0D+00

    xtab(1:8) =   (/ z, a, z, z, c, c, z, c /)
    ytab(1:8) =   (/ z, z, a, z, c, z, c, c /)
    ztab(1:8) =   (/ z, z, z, a, z, c, c, c /)
    weight(1:8) = (/ b, b, b, b, d, d, d, d /)
!
!  Newton Cotes #4.
!
  else if ( rule == 7 ) then

    a =   0.25D+00
    b =   0.50D+00
    c =   0.75D+00
    d =   1.00D+00
    e =  -5.0D+00 / 420.0D+00
    f = -12.0D+00 / 420.0D+00
    g =  16.0D+00 / 420.0D+00
    h = 128.0D+00 / 420.0D+00
    z =   0.0D+00

    xtab(1:35) =   (/ z, d, z, z, a, z, z, c, c, c, z, a, z, z, a, z, b, z, z, &
                      b, b, z, a, b, a, a, b, z, b, z, a, a, z, a, a /)
    ytab(1:35) =   (/ z, z, d, z, z, a, z, z, a, z, c, c, c, z, z, a, z, b, z, &
                      b, z, b, a, a, b, z, z, a, a, b, b, z, a, a, a /)
    ztab(1:35) =   (/ z, z, z, d, z, z, a, z, z, a, z, z, a, c, c, c, z, z, b, &
                      z, b, b, z, z, z, a, a, a, a, a, a, b, b, b, a /)

    weight(1:35) = (/ e, e, e, e, g, g, g, g, g, g, g, g, g, g, g, g, f, f, f, &
                      f, f, f, g, g, g, g, g, g, g, g, g, g, g, g, h /)
!
!  Keast Rule of order 11
!
  else if ( rule == 8 ) then

    a =  0.25D+00
    b =  11.0D+00 /    14.0D+00
    c =   1.0D+00 /    14.0D+00
    d =  0.25D+00 * ( 1.0D+00 + sqrt ( 5.0D+00 / 14.0D+00 ) )
    e =  0.25D+00 * ( 1.0D+00 - sqrt ( 5.0D+00 / 14.0D+00 ) )
    f = -74.0D+00 /  5625.0D+00
    g = 343.0D+00 / 45000.0D+00
    h =  56.0D+00 /  2250.0D+00

    xtab(1:11) =   (/ a, b, c, c, c, d, d, d, e, e, e /)
    ytab(1:11) =   (/ a, c, b, c, c, d, e, e, d, d, e /)
    ztab(1:11) =   (/ a, c, c, b, c, e, d, e, d, e, d /)

    weight(1:11) = (/ f, g, g, g, g, h, h, h, h, h, h /)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRA_UNIT_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
    stop

  end if

  return
end
subroutine tetra_unit_size ( rule, order )

!*****************************************************************************80
!
!! TETRA_UNIT_SIZE sizes quadrature rules in the unit tetrahedron.
!
!  Integration region:
!
!      0 <= X,
!    and
!      0 <= Y,
!    and
!      0 <= Z, 
!    and
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hermann Engels,
!    Numerical Quadrature and Cubature,
!    Academic Press, 1980,
!    ISBN: 012238850X,
!    LC: QA299.3E5.
!
!    Patrick Keast,
!    Moderate Degree Tetrahedral Quadrature Formulas,
!    Computer Methods in Applied Mechanics and Engineering,
!    Volume 55, Number 3, May 1986, pages 339-348.
!
!    Olgierd Zienkiewicz,
!    The Finite Element Method,
!    Sixth Edition,
!    Butterworth-Heinemann, 2005,
!    ISBN: 0750663200,
!    LC: TA640.2.Z54
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!     1, order 1, precision 0, Newton Cotes formula #0, Zienkiewicz #1.
!     2, order 4, precision 1, Newton Cotes formula #1.
!     3, order 4, precision 2, Zienkiewicz #2.
!     4, order 10, precision 2, Newton Cotes formula #2
!     5, order 5, precision 3, Zienkiewicz #3.
!     6, order 8, precision 3, Newton Cotes formula #3.
!     7, order 35, precision 4, Newton Cotes formula #4.
!     8, order 11, precision 4, a Keast rule.
!
!    Output, integer ( kind = 4 ) ORDER, the order of the rule.
!
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ) rule
!
!  Newton Cotes #0.
!
  if ( rule == 1 ) then

    order = 1
!
!  Newton Cotes #1.
!
  else if ( rule == 2 ) then

    order = 4
!
!  Zienkiewicz #2.
!
  else if ( rule == 3 ) then

    order = 4
!
!  Newton Cotes #2.
!
  else if ( rule == 4 ) then

    order = 10
!
!  Zienkiewicz #3.
!
  else if ( rule == 5 ) then

    order = 5
!
!  Newton Cotes #3.
!  (This is actually formally a 20 point rule, but with 12 zero coefficients!)
!
  else if ( rule == 6 ) then

    order = 8
!
!  Newton Cotes #4.
!
  else if ( rule == 7 ) then

    order = 35
!
!  Keast Rule of order 11
!
  else if ( rule == 8 ) then

    order = 11

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRA_UNIT_SIZE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
    stop

  end if

  return
end
subroutine tetra_unit_sum ( func, order, xtab, ytab, ztab, weight, result )

!*****************************************************************************80
!
!! TETRA_UNIT_SUM carries out a quadrature rule in the unit tetrahedron in 3D.
!
!  Integration region:
!
!      0 <= X,
!    and
!      0 <= Y,
!    and
!      0 <= Z, 
!    and
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of three variables which is to be integrated,
!    of the form:
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), ZTAB(ORDER), the
!    abscissas.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) tetra_unit_volume
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) ztab(order)

  quad = 0.0D+00

  do i = 1, order
    quad = quad + weight(i) * func ( xtab(i), ytab(i), ztab(i) )
  end do

  volume = tetra_unit_volume ( )
  result = quad * volume

  return
end
function tetra_unit_volume ( )

!*****************************************************************************80
!
!! TETRA_UNIT_VOLUME returns the volume of the unit tetrahedron in 3D.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X,
!      0 <= Y,
!      0 <= Z, 
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) TETRA_UNIT_VOLUME, the volume.
!
  implicit none

  real ( kind = 8 ) tetra_unit_volume

  tetra_unit_volume = 1.0D+00 / 6.0D+00

  return
end
function tetra_volume ( x, y, z )

!*****************************************************************************80
!
!! TETRA_VOLUME computes the volume of a tetrahedron in 3D.
!
!  Integration region:
!
!    Points inside a tetrahedron whose four vertices are given.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(4), Y(4), Z(4), the vertices.
!
!    Output, real ( kind = 8 ) TETRA_VOLUME, the volume of the tetrahedron.
!
  implicit none

  real ( kind = 8 ) parallelipiped_volume_3d
  real ( kind = 8 ) tetra_unit_volume
  real ( kind = 8 ) tetra_volume
  real ( kind = 8 ) volume
  real ( kind = 8 ) x(4)
  real ( kind = 8 ) y(4)
  real ( kind = 8 ) z(4)

  volume = parallelipiped_volume_3d ( x, y, z )

  tetra_volume = volume * tetra_unit_volume ( )

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
subroutine torus_1 ( func, r1, r2, n, result )

!*****************************************************************************80
!
!! TORUS_1 approximates an integral on the surface of a torus in 3D.
!
!  Integration region:
!
!    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z = R2 * R2.
!
!  Discussion:
!
!    An (N+1)*(N+2) point N-th degree formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of three variables which is to be integrated,
!    of the form:
!      function func ( x, y, z )
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, real ( kind = 8 ) R1, R2, the two radii that define the torus.
!
!    Input, integer ( kind = 4 ) N, defines the degree of the formula
!    used to approximate the integral.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) ct1
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) result
  real ( kind = 8 ) st1
  real ( kind = 8 ) torus_area_3d
  real ( kind = 8 ) u
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  w = 1.0D+00 / ( r1 * real ( ( n + 1 ) * ( n + 2 ), kind = 8 ) )
  quad = 0.0D+00

  do i = 1, n+1

    angle = 2.0D+00 * pi * real ( i, kind = 8 ) / real ( n + 1, kind = 8 )
    ct1 = cos ( angle )
    st1 = sin ( angle )

    do j = 1, n+2

      angle = 2.0D+00 * pi * real ( j, kind = 8 ) / real ( n + 2, kind = 8 )
      u = r1 + r2 * cos ( angle )
      x = u * ct1
      y = u * st1
      z = r2 * sin ( angle )

      quad = quad + w * u * func ( x, y, z )

    end do

  end do

  volume = torus_area_3d ( r1, r2 )
  result = quad * volume

  return
end
subroutine torus_14s ( func, r1, r2, result )

!*****************************************************************************80
!
!! TORUS_14S approximates an integral inside a torus in 3D.
!
!  Integration region:
!
!    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
!
!  Discussion:
!
!    A 960 point 14-th degree formula is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of three variables which is to be integrated,
!    of the form:
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, real ( kind = 8 ) R1, R2, the two radii that define the torus.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 4

  real ( kind = 8 ) angle
  real ( kind = 8 ) ct
  real ( kind = 8 ) cth
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ), save, dimension ( order ) :: r = (/ &
    0.263499230D+00, 0.574464514D+00, 0.818529487D+00, 0.964659606D+00 /)
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) result
  real ( kind = 8 ) st
  real ( kind = 8 ) sth
  real ( kind = 8 ) torus_volume_3d
  real ( kind = 8 ) u
  real ( kind = 8 ) volume
  real ( kind = 8 ), save, dimension ( order ) :: weight = (/ &
    0.086963711D+00, 0.163036289D+00, 0.163036289D+00, 0.086963711D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  quad = 0.0D+00

  do n = 1, 15

    angle = 2.0D+00 * pi * real ( n, kind = 8 ) / 15.0D+00
    cth = cos ( angle )
    sth = sin ( angle )

    do i = 1, 16

      angle = 2.0D+00 * pi * real ( i, kind = 8 ) / 16.0D+00
      ct = cos ( angle )
      st = sin ( angle )

      do j = 1, order
        u = r1 + r(j) * ct * r2
        x = u * cth
        y = u * sth
        z = r(j) * st * r2
        quad = quad + u * weight(j) * func ( x, y, z ) / ( 120.0D+00 * r1 )
      end do

    end do

  end do

  volume = torus_volume_3d ( r1, r2 )
  result = quad * volume

  return
end
subroutine torus_5s2 ( func, r1, r2, result )

!*****************************************************************************80
!
!! TORUS_5S2 approximates an integral inside a torus in 3D.
!
!  Integration region:
!
!    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
!
!  Discussion:
!
!    A 24 point, 5-th degree formula is used, Stroud number TOR3-S2:5-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of three variables which is to be integrated,
!    of the form:
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, real ( kind = 8 ) R1, R2, the two radii that define the torus.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) cs
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) result
  real ( kind = 8 ) sn
  real ( kind = 8 ) torus_volume_3d
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) u3
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  w = 1.0D+00 / 24.0D+00

  quad = 0.0D+00

  u1 = sqrt ( r1 * r1 + 0.5D+00 * r2 * r2 )
  u2 = sqrt ( r1 * r1 + sqrt ( 2.0D+00 ) * r1 * r2 + r2 * r2 )
  u3 = sqrt ( r1 * r1 - sqrt ( 2.0D+00 ) * r1 * r2 + r2 * r2 )

  do i = 1, 6

    angle = 2.0D+00 * pi * real ( i, kind = 8 ) / 6.0D+00
    cs = cos ( angle )
    sn = sin ( angle )

    x = u1 * cs
    y = u1 * sn
    z = r2 / sqrt ( 2.0D+00 )
    quad = quad + w * func ( x, y, z )

    x = u1 * cs
    y = u1 * sn
    z = -r2 / sqrt ( 2.0D+00 )
    quad = quad + w * func ( x, y, z )

    x = u2 * cs
    y = u2 * sn
    z = 0.0D+00
    quad = quad + w * func ( x, y, z )

    x = u3 * cs
    y = u3 * sn
    z = 0.0D+00
    quad = quad + w * func ( x, y, z )

  end do

  volume = torus_volume_3d ( r1, r2 )
  result = quad * volume

  return
end
subroutine torus_6s2 ( func, r1, r2, result )

!*****************************************************************************80
!
!! TORUS_6S2 approximates an integral inside a torus in 3D.
!
!  Integration region:
!
!    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
!
!  Discussion:
!
!    An 84 point 6-th degree formula is used, Stroud number TOR3-S2:6-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of three variables which is to be integrated,
!    of the form:
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, real ( kind = 8 ) R1, R2, the two radii that define the torus.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 2

  real ( kind = 8 ) cth
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) result
  real ( kind = 8 ), save, dimension ( order ) :: s = (/ &
    0.322914992D+00, 0.644171310D+00 /)
  real ( kind = 8 ) sth
  real ( kind = 8 ) torus_volume_3d
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ), save, dimension ( order ) :: weight = (/ &
    0.387077796D+00, 0.165609800D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  w = 1.0D+00 / ( 7.0D+00 * r1 * pi )

  quad = 0.0D+00

  do n = 1, 7

    u = 0.5D+00 * sqrt ( 3.0D+00 ) * r2
    cth = cos ( 2.0D+00 * pi * real ( n, kind = 8 ) / 7.0D+00 )
    sth = sin ( 2.0D+00 * pi * real ( n, kind = 8 ) / 7.0D+00 )

    do i = 1, 2

      u = -u

      x = ( r1 + u ) * cth
      y = ( r1 + u ) * sth
      z = 0.0D+00
      quad = quad + 0.232710567D+00 * w * ( r1 + u ) * func ( x, y, z )

      x = r1 * cth
      y = r1 * sth
      z = u
      quad = quad + 0.232710567D+00 * w * r1 * func ( x, y, z )

    end do

    do k = 1, order

      u = s(k) * r2
      v = u

      do i = 1, 2

        u = -u

        do j = 1, 2

          v = -v

          x = ( r1 + u ) * cth
          y = ( r1 + u ) * sth
          z = v
          quad = quad + weight(k) * w * ( r1 + u ) * func ( x, y, z )

        end do
      end do
    end do
  end do

  volume = torus_volume_3d ( r1, r2 )
  result = quad * volume

  return
end
function torus_area_3d ( r1, r2 )

!*****************************************************************************80
!
!! TORUS_AREA_3D returns the area of a torus in 3D.
!
!  Integration region:
!
!    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z = R2 * R2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, the two radii that define the torus.
!
!    Output, real ( kind = 8 ) TORUS_AREA_3D, the area of the torus.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) torus_area_3d

  torus_area_3d = 4.0D+00 * pi * pi * r1 * r2

  return
end
subroutine torus_square_14c ( func, r1, r2, result )

!*****************************************************************************80
!
!! TORUS_SQUARE_14C approximates an integral in a "square" torus in 3D.
!
!  Discussion:
!
!    A 14-th degree 960 point formula is used.
!
!  Integration region:
!
!      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
!    and
!      -R2 <= Z <= R2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of three variables which is to be integrated, of the form:
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, real ( kind = 8 ) R1, R2, the radii that define the torus.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: order = 8

  real ( kind = 8 ) angle
  real ( kind = 8 ) cth
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) result
  real ( kind = 8 ) rtab(order)
  real ( kind = 8 ) sth
  real ( kind = 8 ) torus_square_volume_3d
  real ( kind = 8 ) u
  real ( kind = 8 ) volume
  real ( kind = 8 ) w
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  call legendre_set ( order, rtab, weight )

  w = 1.0D+00 / ( 60.0D+00 * r1 )
  quad = 0.0D+00

  do n = 1, 15

    angle = 2.0D+00 * pi * real ( n, kind = 8 ) / 15.0D+00
    cth = cos ( angle )
    sth = sin ( angle )

    do i = 1, order

      u = r1 + rtab(i) * r2
      x = u * cth
      y = u * sth

      do j = 1, order
        z = rtab(j) * r2
        quad = quad + u * w * weight(i) * weight(j) * func ( x, y, z )
      end do

    end do

  end do

  volume = torus_square_volume_3d ( r1, r2 )
  result = quad * volume

  return
end
subroutine torus_square_5c2 ( func, r1, r2, result )

!*****************************************************************************80
!
!! TORUS_SQUARE_5C2 approximates an integral in a "square" torus in 3D.
!
!  Integration region:
!
!      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
!    and
!      -R2 <= Z <= R2.
!
!  Discussion:
!
!    A 24 point 5-th degree formula is used, Stroud number TOR3-C2:5-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of three variables which is to be integrated,
!    of the form:
!      function func ( x, y, z )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!      real ( kind = 8 ) z
!
!    Input, real ( kind = 8 ) R1, the primary radius of the torus.
!
!    Input, real ( kind = 8 ) R2, one-half the length of a side of the
!    square cross-section.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  real ( kind = 8 ), parameter :: b1 = 5.0D+00 / 108.0D+00
  real ( kind = 8 ), parameter :: b2 = 4.0D+00 / 108.0D+00
  real ( kind = 8 ) cs
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) result
  real ( kind = 8 ) sn
  real ( kind = 8 ) torus_square_volume_3d
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) u3
  real ( kind = 8 ) v
  real ( kind = 8 ) volume
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  quad = 0.0D+00

  u1 = sqrt ( r1 * r1 + r2 * r2 )

  v = r2 * sqrt ( 0.6D+00 )

  u2 = sqrt ( r1 * r1 - sqrt ( 3.0D+00 ) * r1 * r2 + r2 * r2 )

  u3 = sqrt ( r1 * r1 + sqrt ( 3.0D+00 ) * r1 * r2 + r2 * r2 )

  do i = 1, 6

    cs = cos ( real ( i, kind = 8 ) * pi / 3.0D+00 )
    sn = sin ( real ( i, kind = 8 ) * pi / 3.0D+00 )

    x = u1 * cs
    y = u1 * sn
    z = v
    quad = quad + b1 * func ( x, y, z )

    z = -v
    quad = quad + b1 * func ( x, y, z )

    x = u2 * cs
    y = u2 * sn
    z = 0.0D+00
    quad = quad + b2 * func ( x, y, z )

    x = u3 * cs
    y = u3 * sn
    z = 0.0D+00
    quad = quad + b2 * func ( x, y, z )

  end do

  volume = torus_square_volume_3d ( r1, r2 )
  result = quad * volume

  return
end
function torus_square_area_3d ( r1, r2 )

!*****************************************************************************80
!
!! TORUS_SQUARE_AREA_3D returns the area of a square torus in 3D.
!
!  Integration region:
!
!      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
!    and
!      -R2 <= Z <= R2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, the two radii that define the torus.
!
!    Output, real ( kind = 8 ) TORUS_SQUARE_AREA_3D, the area of the torus.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) torus_square_area_3d

  torus_square_area_3d = 16.0D+00 * pi * r1 * r2

  return
end
function torus_square_volume_3d ( r1, r2 )

!*****************************************************************************80
!
!! TORUS_SQUARE_VOLUME_3D returns the volume of a square torus in 3D.
!
!  Integration region:
!
!      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
!    and
!      -R2 <= Z <= R2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, the two radii that define the torus.
!
!    Output, real ( kind = 8 ) TORUS_SQUARE_VOLUME_3D, the volume of the torus.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) torus_square_volume_3d

  torus_square_volume_3d = 8.0D+00 * pi * r1 * r2 * r2

  return
end
function torus_volume_3d ( r1, r2 )

!*****************************************************************************80
!
!! TORUS_VOLUME_3D returns the volume of a torus in 3D.
!
!  Integration region:
!
!    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z = R2 * R2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R1, R2, the two radii that define the torus.
!
!    Output, real ( kind = 8 ) TORUS_VOLUME_3D, the volume of the torus.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) torus_volume_3d

  torus_volume_3d = 2.0D+00 * pi * pi * r1 * r2 * r2

  return
end
subroutine triangle_rule_adjust ( xval, yval, order, xtab, ytab, weight, &
  xtab2, ytab2, weight2 )

!*****************************************************************************80
!
!! TRIANGLE_RULE_ADJUST adjusts a unit quadrature rule to an arbitrary triangle.
!
!  Integration region:
!
!      (X,Y) = ALPHA * (X1,Y1) + BETA * (X2,Y2) + ( 1 - ALPHA - BETA ) * (X3,Y3)
!    and
!      0 <= ALPHA <= 1 - BETA
!    and
!      0 <= BETA <= 1 - ALPHA
!
!  Discussion:
!
!    This routine accepts as input abscissas and weights appropriate for
!    quadrature in the unit triangle, and returns abscissas and weights
!    appropriate for quadrature in a given triangle.
!
!    Once this routine has been called, an integral over the given triangle
!    can be approximated as:
!
!      QUAD = sum ( 1 <= I <= ORDER ) WTAB2(I) * FUNC ( XTAB2(I), YTAB2(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVAL(3), YVAL(3), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas for
!    the unit triangle.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights for the unit triangle.
!
!    Output, real ( kind = 8 ) XTAB2(ORDER), YTAB2(ORDER), the adjusted
!    abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT2(ORDER), the adjusted weights.
!
  implicit none

  integer ( kind = 4 ) order

  integer ( kind = 4 ) i
  real ( kind = 8 ) triangle_volume
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) weight2(order)
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) xtab2(order)
  real ( kind = 8 ) xval(3)
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) ytab2(order)
  real ( kind = 8 ) yval(3)

  volume = triangle_volume ( xval, yval )

  do i = 1, order

    xtab2(i) =             xtab(i)             * xval(1) &
             +                       ytab(i)   * xval(2) &
             + ( 1.0D+00 - xtab(i) - ytab(i) ) * xval(3)

    ytab2(i) =             xtab(i)             * yval(1) &
                                   + ytab(i)   * yval(2) &
             + ( 1.0D+00 - xtab(i) - ytab(i) ) * yval(3)

    weight2(i) = weight(i) * 2.0D+00 * volume

  end do

  return
end
subroutine triangle_sub ( func, xval, yval, nsub, order, xtab, ytab, weight, &
  result )

!*****************************************************************************80
!
!! TRIANGLE_SUB carries out quadrature over subdivisions of a triangular region.
!
!  Integration region:
!
!      (X,Y) =       ALPHA          * ( XVAL(1), YVAL(1) )
!            +               BETA   * ( XVAL(2), YVAL(2) )
!            + ( 1 - ALPHA - BETA ) * ( XVAL(3), YVAL(3) )
!    and
!      0 <= ALPHA <= 1 - BETA
!    and
!      0 <= BETA <= 1 - ALPHA
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function of
!    two variables which is to be integrated, of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, real ( kind = 8 ) XVAL(3), YVAL(3), the coordinates of the 
!    triangle vertices.
!
!    Input, integer ( kind = 4 ) NSUB, the number of subdivisions of each side 
!    of the input triangle to be made.  NSUB = 1 means no subdivisions are made.
!    NSUB = 3 means that each side of the triangle is subdivided into
!    three portions, and that the original triangle is subdivided into
!    NSUB * NSUB triangles.  NSUB must be at least 1.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  real ( kind = 8 ) triangle_volume
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) xval(3)
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) yval(3)
!
!  Initialize RESULT, the approximate integral.
!
  result = 0.0D+00
!
!  NSUB must be positive.
!
  if ( nsub <= 0 ) then
    return
  end if
!
!  Initialize QUAD, the quadrature sum.
!
  quad = 0.0D+00
!
!  The sub-triangles can be grouped into NSUB strips.
!
  do i = 1, nsub

    temp1 = 0.0D+00
    temp2 = real ( i, kind = 8 ) / real ( nsub, kind = 8 )

    x2 = xval(2) + temp1 * ( xval(3) - xval(2) ) &
                 + temp2 * ( xval(1) - xval(2) )

    y2 = yval(2) + temp1 * ( yval(3) - yval(2) ) &
                 + temp2 * ( yval(1) - yval(2) )

    temp1 = 0.0D+00
    temp2 = real ( i - 1, kind = 8 ) / real ( nsub, kind = 8 )

    x3 = xval(2) + temp1 * ( xval(3) - xval(2) ) &
                 + temp2 * ( xval(1) - xval(2) )

    y3 = yval(2) + temp1 * ( yval(3) - yval(2) ) &
                 + temp2 * ( yval(1) - yval(2) )
!
!  There are 2*I-1 triangles in strip number I.
!  The next triangle in the strip shares two nodes with the previous one.
!  Compute its corners, (X1,Y1), (X2,Y2), (X3,Y3).
!
    do j = 1, 2*i-1

      x1 = x2
      y1 = y2
      x2 = x3
      y2 = y3
      temp1 = real ( ( ( j + 1 ) / 2 ), kind = 8 ) / real ( nsub, kind = 8 )
      temp2 = real ( ( i - 1 - ( j / 2 ) ), kind = 8 ) / real ( nsub, kind = 8 )

      x3 = xval(2) + temp1 * ( xval(3) - xval(2) ) &
                   + temp2 * ( xval(1) - xval(2) )

      y3 = yval(2) + temp1 * ( yval(3) - yval(2) ) &
                   + temp2 * ( yval(1) - yval(2) )
!
!  Now integrate over the triangle, mapping the points ( XTAB(K), YTAB(K) )
!  into the triangle.
!
      do k = 1, order

        x = x2 + xtab(k) * ( x3 - x2 ) + ytab(k) * ( x1 - x2 )
        y = y2 + xtab(k) * ( y3 - y2 ) + ytab(k) * ( y1 - y2 )
        quad = quad + weight(k) * func ( x, y )

       end do

    end do

  end do

  volume = triangle_volume ( xval, yval ) / real ( nsub * nsub, kind = 8 )
  result = quad * volume

  return
end
subroutine triangle_sum ( func, xval, yval, order, xtab, ytab, weight, result )

!*****************************************************************************80
!
!! TRIANGLE_SUM carries out a unit quadrature rule in an arbitrary triangle.
!
!  Integration region:
!
!      (X,Y) =       ALPHA          * (X1,Y1) 
!            +               BETA   * (X2,Y2) 
!            + ( 1 - ALPHA - BETA ) * (X3,Y3)
!    and
!      0 <= ALPHA <= 1 - BETA
!    and
!      0 <= BETA <= 1 - ALPHA
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function of
!    two variables which is to be integrated, of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, real ( kind = 8 ) XVAL(3), YVAL(3), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) triangle_volume
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) xval(3)
  real ( kind = 8 ) y
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) yval(3)

  quad = 0.0D+00

  do i = 1, order

    x =             xtab(i)             * xval(1) &
      +                       ytab(i)   * xval(2) &
      + ( 1.0D+00 - xtab(i) - ytab(i) ) * xval(3)

    y =             xtab(i)             * yval(1) &
      +                       ytab(i)   * yval(2) &
      + ( 1.0D+00 - xtab(i) - ytab(i) ) * yval(3)

    quad = quad + weight(i) * func ( x, y )

  end do

  volume = triangle_volume ( xval, yval )
  result = quad * volume

  return
end
subroutine triangle_sum_adjusted ( func, order, xtab, ytab, weight, result )

!*****************************************************************************80
!
!! TRIANGLE_SUM_ADJUSTED carries out an adjusted quadrature rule in a triangle.
!
!  Integration region:
!
!      (X,Y) =       ALPHA          * (X1,Y1) 
!                          + BETA   * (X2,Y2) 
!            + ( 1 - ALPHA - BETA ) * (X3,Y3)
!    and
!      0 <= ALPHA <= 1 - BETA
!    and
!      0 <= BETA <= 1 - ALPHA
!
!  Discussion:
!
!    It is assumed that a quadrature rule approprate for the unit triangle
!    was generated, and then adjusted to a particular triangle by calling
!    TRIANGLE_RULE_ADJUST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied function of
!    two variables which is to be integrated, of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) ytab(order)

  result = 0.0D+00

  do i = 1, order
    result = result + weight(i) * func ( xtab(i), ytab(i) )
  end do

  return
end
subroutine triangle_unit_product_set ( rule, order, xtab, ytab, weight )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_PRODUCT_SET sets a product rule on the unit triangle.
!
!  Discussion:
!
!    For a given order of accuracy, a product rule on a triangle usually
!    uses more points than necessary.  That is, there is usually a rule
!    of the same order that uses fewer points.
!
!    However, one advantage of product rules is that a rule of any
!    desired order can be generated automatically.
!   
!    The integration region is:
!
!      0 <= X,
!    and
!      0 <= Y, 
!    and
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the order of the 1D rule.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ) rule

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order0
  integer ( kind = 4 ) order1
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) weight0(rule)
  real ( kind = 8 ) weight1(rule)
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) xtab0(rule)
  real ( kind = 8 ) xtab1(rule)
  real ( kind = 8 ) ytab(order)

  a = -1.0D+00
  b = +1.0D+00
  c =  0.0D+00
  d = +1.0D+00

  order0 = rule
  call legendre_set ( order0, xtab0, weight0 )
  call rule_adjust ( a, b, c, d, order0, xtab0, weight0 )

  order1 = rule
  call legendre_set_x1 ( order1, xtab1, weight1 )
  call rule_adjust ( a, b, c, d, order1, xtab1, weight1 )

  k = 0
  do j = 1, order1
    do i = 1, order0
      k = k + 1
      xtab(k) = 1.0D+00 - xtab1(j)
      ytab(k) = xtab0(i) * xtab1(j)
      weight(k) = weight0(i) * weight1(j)
    end do
  end do

  return
end
subroutine triangle_unit_product_size ( rule, order )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_PRODUCT_SIZE sizes a product rule on the unit triangle.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X,
!    and
!      0 <= Y, 
!    and
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the order of the 1D rule.
!
!    Output, integer ( kind = 4 ) ORDER, the order of the rule. 
!
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ) rule

  order = rule * rule

  return
end
subroutine triangle_unit_set ( rule, order, xtab, ytab, weight )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_SET sets a quadrature rule in the unit triangle.
!
!  Discussion:
!
!    The user is responsible for determining the value of ORDER,
!    and appropriately dimensioning the arrays XTAB, YTAB and
!    WEIGHT so that they can accommodate the data.
!
!    The value of ORDER for each rule can be found by invoking
!    the function TRIANGLE_RULE_SIZE.
!
!  Integration region:
!
!      0 <= X,
!    and
!      0 <= Y, 
!    and
!      X + Y <= 1.
!
!  Graph:
!
!      ^
!    1 | *
!      | |\
!    Y | | \
!      | |  \
!    0 | *---*
!      +------->
!        0 X 1
!
!   The rules are accessed by an index number, RULE.  The indices,
!   and the descriptions of the corresponding rules, are:
!
!     1, ORDER =  1, precision 1, Zienkiewicz #1.
!     2, ORDER =  2, precision 1, (the "vertex rule").
!     3, ORDER =  3, precision 2, Strang and Fix formula #1.
!     4, ORDER =  3, precision 2, Strang and Fix formula #2,
!                                 Zienkiewicz #2.
!     5, ORDER =  4, precision 3, Strang and Fix formula #3,
!                                 Zienkiewicz #3.
!     6, ORDER =  6, precision 3, Strang and Fix formula #4.
!     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
!     8, ORDER =  6, precision 4, Strang and Fix formula #5.
!     9, ORDER =  7, precision 4, Strang and Fix formula #6.
!    10, ORDER =  7, precision 5, Strang and Fix formula #7,
!                                 Stroud formula T2:5-1, 
!                                 Zienkiewicz #4, 
!                                 Schwarz Table 2.2.
!    11, ORDER =  9, precision 6, Strang and Fix formula #8.
!    12, ORDER = 12, precision 6, Strang and Fix formula #9.
!    13, ORDER = 13, precision 7, Strang and Fix formula #10.
!        Note that there is a typographical error in Strang and Fix
!        which lists the value of the XSI(3) component of the
!        last generator point as 0.4869... when it should be 0.04869...
!    14, ORDER =  7, precision 3.
!    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
!    16, ORDER = 64, precision 15, triangular product Gauss rule.
!    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
!    18, ORDER = 19, precision 9, from TRIEX, ACM TOMS #612.
!    19, ORDER = 28, precision 11, from TRIEX, ACM TOMS #612.
!    20, ORDER = 37, precision 13, from ACM TOMS #706.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jarle Berntsen, Terje Espelid,
!    Algorithm 706,
!    DCUTRI: an algorithm for adaptive cubature over a collection of triangles, 
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, September 1992, pages 329-342.
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!    Dirk Laurie,
!    Algorithm 584,
!    CUBTRI, Automatic Cubature Over a Triangle,
!    ACM Transactions on Mathematical Software,
!    Volume 8, Number 2, 1982, pages 210-218.
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!    Hans Rudolf Schwarz,
!    Finite Element Methods,
!    Academic Press, 1988,
!    ISBN: 0126330107,
!    LC: TA347.F5.S3313.
!
!    Gilbert Strang, George Fix,
!    An Analysis of the Finite Element Method,
!    Cambridge, 1973,
!    ISBN: 096140888X,
!    LC: TA335.S77.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!    Olgierd Zienkiewicz,
!    The Finite Element Method,
!    Sixth Edition,
!    Butterworth-Heinemann, 2005,
!    ISBN: 0750663200,
!    LC: TA640.2.Z54
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order2
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  integer ( kind = 4 ) rule
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) w4
  real ( kind = 8 ) w5
  real ( kind = 8 ) w6
  real ( kind = 8 ) w7
  real ( kind = 8 ) w8
  real ( kind = 8 ) w9
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) weight1(8)
  real ( kind = 8 ) weight2(8)
  real ( kind = 8 ) wx
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) xtab1(8)
  real ( kind = 8 ) xtab2(8)
  real ( kind = 8 ) y
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) z
!
!  1 point, precision 1.
!
  if ( rule == 1 ) then

    xtab(1)   = 0.33333333333333333333D+00

    ytab(1)   = 0.33333333333333333333D+00

    weight(1) = 1.00000000000000000000D+00
!
!  3 points, precision 1, the "vertex rule".
!
  else if ( rule == 2 ) then

    xtab(1) =   1.00000000000000000000D+00
    xtab(2) =   0.00000000000000000000D+00
    xtab(3) =   0.00000000000000000000D+00

    ytab(1) =   0.00000000000000000000D+00
    ytab(2) =   1.00000000000000000000D+00
    ytab(3) =   0.00000000000000000000D+00

    weight(1) = 0.33333333333333333333D+00
    weight(2) = 0.33333333333333333333D+00
    weight(3) = 0.33333333333333333333D+00
!
!  3 points, precision 2, Strang and Fix formula #1.
!
  else if ( rule == 3 ) then

    xtab(1)   = 0.66666666666666666667D+00
    xtab(2)   = 0.16666666666666666667D+00
    xtab(3)   = 0.16666666666666666667D+00

    ytab(1)   = 0.16666666666666666667D+00
    ytab(2)   = 0.66666666666666666667D+00
    ytab(3)   = 0.16666666666666666667D+00

    weight(1) = 0.33333333333333333333D+00
    weight(2) = 0.33333333333333333333D+00
    weight(3) = 0.33333333333333333333D+00
!
!  3 points, precision 2, Strang and Fix formula #2.
!
  else if ( rule == 4 ) then

    xtab(1)   = 0.50000000000000000000D+00
    xtab(2)   = 0.50000000000000000000D+00
    xtab(3)   = 0.00000000000000000000D+00

    ytab(1)   = 0.00000000000000000000D+00
    ytab(2)   = 0.50000000000000000000D+00
    ytab(3)   = 0.50000000000000000000D+00

    weight(1) = 0.33333333333333333333D+00
    weight(2) = 0.33333333333333333333D+00
    weight(3) = 0.33333333333333333333D+00
!
!  4 points, precision 3, Strang and Fix formula #3.
!
  else if ( rule == 5 ) then

    a =   6.0D+00
    b =  10.0D+00
    c =  18.0D+00
    d =  25.0D+00
    e = -27.0D+00
    f =  30.0D+00
    g =  48.0D+00

    xtab(1:4) =   (/ b, c, a, a /) / f
    ytab(1:4) =   (/ b, a, c, a /) / f
    weight(1:4) = (/ e, d, d, d /) / g
!
!  6 points, precision 3, Strang and Fix formula #4.
!
  else if ( rule == 6 ) then

    a = 0.659027622374092D+00
    b = 0.231933368553031D+00
    c = 0.109039009072877D+00

    xtab(1:6) =   (/ a, a, b, b, c, c /)
    ytab(1:6) =   (/ b, c, a, c, a, b /)

    weight(1) = 0.16666666666666666667D+00
    weight(2) = 0.16666666666666666667D+00
    weight(3) = 0.16666666666666666667D+00
    weight(4) = 0.16666666666666666667D+00
    weight(5) = 0.16666666666666666667D+00
    weight(6) = 0.16666666666666666667D+00
!
!  6 points, precision 3, Stroud T2:3-1.
!
  else if ( rule == 7 ) then

    a = 0.0D+00
    b = 0.5D+00
    c = 2.0D+00 /  3.0D+00
    d = 1.0D+00 /  6.0D+00
    v = 1.0D+00 / 30.0D+00
    w = 3.0D+00 / 10.0D+00

    xtab(1:6) =   (/ a, b, b, c, d, d /)
    ytab(1:6) =   (/ b, a, b, d, c, d /)
    weight(1:6) = (/ v, v, v, w, w, w /)
!
!  6 points, precision 4, Strang and Fix, formula #5.
!
  else if ( rule == 8 ) then

    a = 0.816847572980459D+00
    b = 0.091576213509771D+00
    c = 0.108103018168070D+00
    d = 0.445948490915965D+00
    v = 0.109951743655322D+00
    w = 0.223381589678011D+00

    xtab(1:6) =   (/ a, b, b, c, d, d /)
    ytab(1:6) =   (/ b, a, b, d, c, d /)
    weight(1:6) = (/ v, v, v, w, w, w /)
!
!  7 points, precision 4, Strang and Fix formula #6.
!
  else if ( rule == 9 ) then

    a = 1.0D+00 / 3.0D+00
    c = 0.736712498968435D+00
    d = 0.237932366472434D+00
    e = 0.025355134551932D+00
    v = 0.375000000000000D+00
    w = 0.104166666666667D+00

    xtab(1:7) =   (/ a, c, c, d, d, e, e /)
    ytab(1:7) =   (/ a, d, e, c, e, c, d /)
    weight(1:7) = (/ v, w, w, w, w, w, w /)
!
!  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1
!
  else if ( rule == 10 ) then

    a = 1.0D+00 / 3.0D+00
    b = ( 9.0D+00 + 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    c = ( 6.0D+00 -           sqrt ( 15.0D+00 ) ) / 21.0D+00
    d = ( 9.0D+00 - 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    e = ( 6.0D+00 +           sqrt ( 15.0D+00 ) ) / 21.0D+00
    u = 0.225D+00
    v = ( 155.0D+00 - sqrt ( 15.0D+00 ) ) / 1200.0D+00
    w = ( 155.0D+00 + sqrt ( 15.0D+00 ) ) / 1200.0D+00

    xtab(1:7) =   (/ a, b, c, c, d, e, e /)
    ytab(1:7) =   (/ a, c, b, c, e, d, e /)
    weight(1:7) = (/ u, v, v, v, w, w, w /)
!
!  9 points, precision 6, Strang and Fix formula #8.
!
  else if ( rule == 11 ) then

    a = 0.124949503233232D+00
    b = 0.437525248383384D+00
    c = 0.797112651860071D+00
    d = 0.165409927389841D+00
    e = 0.037477420750088D+00

    u = 0.205950504760887D+00
    v = 0.063691414286223D+00

    xtab(1:9) =   (/ a, b, b, c, c, d, d, e, e /)
    ytab(1:9) =   (/ b, a, b, d, e, c, e, c, d /)
    weight(1:9) = (/ u, u, u, v, v, v, v, v, v /)
!
!  12 points, precision 6, Strang and Fix, formula #9.
!
  else if ( rule == 12 ) then

    a = 0.873821971016996D+00
    b = 0.063089014491502D+00
    c = 0.501426509658179D+00
    d = 0.249286745170910D+00
    e = 0.636502499121399D+00
    f = 0.310352451033785D+00
    g = 0.053145049844816D+00

    u = 0.050844906370207D+00
    v = 0.116786275726379D+00
    w = 0.082851075618374D+00

    xtab(1:12) =   (/ a, b, b, c, d, d, e, e, f, f, g, g /)
    ytab(1:12) =   (/ b, a, b, d, c, d, f, g, e, g, e, f /)
    weight(1:12) = (/ u, u, u, v, v, v, w, w, w, w, w, w /)
!
!  13 points, precision 7, Strang and Fix, formula #10.
!
!  Note that there is a typographical error in Strang and Fix
!  which lists the value of the XSI(3) component of the
!  last generator point as 0.4869... when it should be 0.04869...
!
  else if ( rule == 13 ) then

    h = 1.0D+00 / 3.0D+00
    a = 0.479308067841923D+00
    b = 0.260345966079038D+00
    c = 0.869739794195568D+00
    d = 0.065130102902216D+00
    e = 0.638444188569809D+00
    f = 0.312865496004875D+00
    g = 0.048690315425316D+00

    w = -0.149570044467670D+00
    t =  0.175615257433204D+00
    u =  0.053347235608839D+00
    v =  0.077113760890257D+00

    xtab(1:13) =   (/ h, a, b, b, c, d, d, e, e, f, f, g, g /)
    ytab(1:13) =   (/ h, b, a, b, d, c, d, f, g, e, g, e, f /)
    weight(1:13) = (/ w, t, t, t, u, u, u, v, v, v, v, v, v /)
!
!  7 points, precision 3.
!
  else if ( rule == 14 ) then

    a = 1.0D+00 / 3.0D+00
    b = 1.0D+00
    c = 0.5D+00
    z = 0.0D+00

    u = 27.0D+00 / 60.0D+00
    v =  3.0D+00 / 60.0D+00
    w =  8.0D+00 / 60.0D+00

    xtab(1:7) =   (/ a, b, z, z, z, c, c /)
    ytab(1:7) =   (/ a, z, b, z, c, z, c /)
    weight(1:7) = (/ u, v, v, v, w, w, w /)
!
!  16 points, precision 5, Stroud T2:7-1.
!
  else if ( rule == 15 ) then
!
!  Legendre rule of order 4.
!
    order2 = 4

    xtab1(1:4) = (/ &
      -0.861136311594052575223946488893D+00, &
      -0.339981043584856264802665759103D+00, &
       0.339981043584856264802665759103D+00, &
       0.861136311594052575223946488893D+00 /)

    weight1(1:4) = (/ &
      0.347854845137453857373063949222D+00, &
      0.652145154862546142626936050778D+00, &
      0.652145154862546142626936050778D+00, &
      0.347854845137453857373063949222D+00 /)

    xtab1(1:order2) = 0.5D+00 * ( xtab1(1:order2) + 1.0D+00 )

    weight2(1) = 0.1355069134D+00
    weight2(2) = 0.2034645680D+00
    weight2(3) = 0.1298475476D+00
    weight2(4) = 0.0311809709D+00

    xtab2(1) = 0.0571041961D+00
    xtab2(2) = 0.2768430136D+00
    xtab2(3) = 0.5835904324D+00
    xtab2(4) = 0.8602401357D+00

    k = 0
    do i = 1, order2
      do j = 1, order2
        k = k + 1
        xtab(k) = xtab2(j)
        ytab(k) = xtab1(i) * ( 1.0D+00 - xtab2(j) )
        weight(k) = weight1(i) * weight2(j)
      end do
    end do
!
!  64 points, precision 15.
!
  else if ( rule == 16 ) then
!
!  Legendre rule of order 8.
!
    order2 = 8

    xtab1(1) = -0.960289856497536231683560868569D+00
    xtab1(2) = -0.796666477413626739591553936476D+00
    xtab1(3) = -0.525532409916328985817739049189D+00
    xtab1(4) = -0.183434642495649804939476142360D+00
    xtab1(5) =  0.183434642495649804939476142360D+00
    xtab1(6) =  0.525532409916328985817739049189D+00
    xtab1(7) =  0.796666477413626739591553936476D+00
    xtab1(8) =  0.960289856497536231683560868569D+00

    weight1(1) = 0.101228536290376259152531354310D+00
    weight1(2) = 0.222381034453374470544355994426D+00
    weight1(3) = 0.313706645877887287337962201987D+00
    weight1(4) = 0.362683783378361982965150449277D+00
    weight1(5) = 0.362683783378361982965150449277D+00
    weight1(6) = 0.313706645877887287337962201987D+00
    weight1(7) = 0.222381034453374470544355994426D+00
    weight1(8) = 0.101228536290376259152531354310D+00

    weight2(1) = 0.00329519144D+00
    weight2(2) = 0.01784290266D+00
    weight2(3) = 0.04543931950D+00
    weight2(4) = 0.07919959949D+00
    weight2(5) = 0.10604735944D+00
    weight2(6) = 0.11250579947D+00
    weight2(7) = 0.09111902364D+00
    weight2(8) = 0.04455080436D+00

    xtab2(1) = 0.04463395529D+00
    xtab2(2) = 0.14436625704D+00
    xtab2(3) = 0.28682475714D+00
    xtab2(4) = 0.45481331520D+00
    xtab2(5) = 0.62806783542D+00
    xtab2(6) = 0.78569152060D+00
    xtab2(7) = 0.90867639210D+00
    xtab2(8) = 0.98222008485D+00

    k = 0
    do j = 1, order2
      do i = 1, order2
        k = k + 1
        xtab(k) = 1.0D+00 - xtab2(j)
        ytab(k) = 0.5D+00 * ( 1.0D+00 + xtab1(i) ) * xtab2(j)
        weight(k) = weight1(i) * weight2(j)
      end do
    end do
!
!  19 points, precision 8, from CUBTRI.
!
  else if ( rule == 17 ) then

    a = 1.0D+00 / 3.0D+00
    b = ( 9.0D+00 + 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    c = ( 6.0D+00 -       sqrt ( 15.0D+00 ) ) / 21.0D+00
    d = ( 9.0D+00 - 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    e = ( 6.0D+00 +       sqrt ( 15.0D+00 ) ) / 21.0D+00
    f = ( 40.0D+00 - 10.0D+00 * sqrt ( 15.0D+00 ) &
      + 10.0D+00 * sqrt ( 7.0D+00 ) + 2.0D+00 * sqrt ( 105.0D+00 ) ) / 90.0D+00
    g = ( 25.0D+00 +  5.0D+00 * sqrt ( 15.0D+00 ) &
      -  5.0D+00 * sqrt ( 7.0D+00 ) - sqrt ( 105.0D+00 ) ) / 90.0D+00
    p = ( 40.0D+00 + 10.0D+00 * sqrt ( 15.0D+00 ) &
      + 10.0D+00 * sqrt ( 7.0D+00 ) - 2.0D+00 * sqrt ( 105.0D+00 ) ) / 90.0D+00
    q = ( 25.0D+00 -  5.0D+00 * sqrt ( 15.0D+00 ) &
      -  5.0D+00 * sqrt ( 7.0D+00 ) + sqrt ( 105.0D+00 ) ) / 90.0D+00
    r = ( 40.0D+00 + 10.0D+00 * sqrt ( 7.0D+00 ) ) / 90.0D+00
    s = ( 25.0D+00 +  5.0D+00 * sqrt ( 15.0D+00 ) - 5.0D+00 * sqrt ( 7.0D+00 ) &
      - sqrt ( 105.0D+00 ) ) / 90.0D+00
    t = ( 25.0D+00 -  5.0D+00 * sqrt ( 15.0D+00 ) - 5.0D+00 * sqrt ( 7.0D+00 ) &
      + sqrt ( 105.0D+00 ) ) / 90.0D+00

    w1 = ( 7137.0D+00 - 1800.0D+00 * sqrt ( 7.0D+00 ) ) / 62720.0D+00
    w2 = -9301697.0D+00 / 4695040.0D+00 - 13517313.0D+00 * sqrt ( 15.0D+00 ) &
      / 23475200.0D+00 + 764885.0D+00 * sqrt ( 7.0D+00 ) / 939008.0D+00 &
      + 198763.0D+00 * sqrt ( 105.0D+00 ) / 939008.0D+00
    w2 = w2 / 3.0D+00
    w3 = -9301697.0D+00 / 4695040.0D+00 + 13517313.0D+00 * sqrt ( 15.0D+00 ) &
      / 23475200.0D+00 &
      + 764885.0D+00 * sqrt ( 7.0D+00 ) / 939008.0D+00 &
      - 198763.0D+00 * sqrt ( 105.0D+00 ) / 939008.0D+00
    w3 = w3 / 3.0D+00
    w4 = ( 102791225.0D+00 - 23876225.0D+00 * sqrt ( 15.0D+00 ) &
      - 34500875.0D+00 * sqrt ( 7.0D+00 ) &
      + 9914825.0D+00 * sqrt ( 105.0D+00 ) ) / 59157504.0D+00
    w4 = w4 / 3.0D+00
    w5 = ( 102791225.0D+00 + 23876225.0D+00 * sqrt ( 15.0D+00 ) &
      - 34500875.0D+00 * sqrt ( 7.0D+00 ) &
      - 9914825D+00 * sqrt ( 105.0D+00 ) ) / 59157504.0D+00
    w5 = w5 / 3.0D+00
    w6 = ( 11075.0D+00 - 3500.0D+00 * sqrt ( 7.0D+00 ) ) / 8064.0D+00
    w6 = w6 / 6.0D+00

    xtab(1:19) =   (/  a,  b,  c,  c,  d,  e,  e,  f,  g,  g,  p,  q,  q, &
                       r,  r,  s,  s,  t,  t /)
    ytab(1:19) =   (/  a,  c,  b,  c,  e,  d,  e,  g,  f,  g,  q,  p,  q, &
                       s,  t,  r,  t,  r,  s /)
    weight(1:19) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
      w6, w6, w6, w6, w6, w6 /)
!
!  19 points, precision 9.
!  Lyness and Jesperson.
!
  else if ( rule == 18 ) then

    a = 1.0D+00 / 3.0D+00
    b =  0.02063496160252593D+00
    c =  0.4896825191987370D+00
    d =  0.1258208170141290D+00
    e =  0.4370895914929355D+00
    f =  0.6235929287619356D+00
    g =  0.1882035356190322D+00
    r =  0.9105409732110941D+00
    s =  0.04472951339445297D+00
    t =  0.7411985987844980D+00
    u =  0.03683841205473626D+00
    v =  0.22196298916076574D+00

    w1 = 0.09713579628279610D+00
    w2 = 0.03133470022713983D+00
    w3 = 0.07782754100477543D+00
    w4 = 0.07964773892720910D+00
    w5 = 0.02557767565869810D+00
    w6 = 0.04328353937728940D+00

    xtab(1:19) =   (/  a,  b,  c,  c,  d,  e,  e,  f,  g,  g,  r,  s,  s, &
      t, t, u, u, v, v /)
    ytab(1:19) =   (/  a,  c,  b,  c,  e,  d,  e,  g,  f,  g,  s,  r,  s, &
      u, v, t, v, t, u /)
    weight(1:19) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
      w6, w6, w6, w6, w6, w6 /)
!
!  28 points, precision 11.
!  Lyness and Jesperson.
!
  else if ( rule == 19 ) then

    a = 1.0D+00 / 3.0D+00
    b = 0.9480217181434233D+00
    c = 0.02598914092828833D+00
    d = 0.8114249947041546D+00
    e = 0.09428750264792270D+00
    f = 0.01072644996557060D+00
    g = 0.4946367750172147D+00
    p = 0.5853132347709715D+00
    q = 0.2073433826145142D+00
    r = 0.1221843885990187D+00
    s = 0.4389078057004907D+00
    t = 0.6779376548825902D+00
    u = 0.04484167758913055D+00
    v = 0.27722066752827925D+00
    w = 0.8588702812826364D+00
    x = 0.0D+00
    y = 0.1411297187173636D+00

    w1 = 0.08797730116222190D+00
    w2 = 0.008744311553736190D+00
    w3 = 0.03808157199393533D+00
    w4 = 0.01885544805613125D+00
    w5 = 0.07215969754474100D+00
    w6 = 0.06932913870553720D+00
    w7 = 0.04105631542928860D+00
    w8 = 0.007362383783300573D+00

    xtab(1:28) =   (/  a,  b,  c,  c,  d,  e,  e,  f,  g,  g,  p,  q,  q, &
       r,  s,  s,  t,  t,  u,  u,  v,  v,  w,  w,  x,  x,  y,  y /)
    ytab(1:28) =   (/  a,  c,  b,  c,  e,  d,  e,  g,  f,  g,  q,  p,  q, &
       s,  r,  s,  u,  v,  t,  v,  t,  u,  x,  y,  w,  y,  w,  x /)
    weight(1:28) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
      w6, w6, w6, w7, w7, w7, w7, w7, w7, w8, w8, w8, w8, w8, w8 /)
!
!  37 points, precision 13.
!
  else if ( rule == 20 ) then

    a = 1.0D+00 / 3.0D+00
    b = 0.950275662924105565450352089520D+00
    c = 0.024862168537947217274823955239D+00
    d = 0.171614914923835347556304795551D+00
    e = 0.414192542538082326221847602214D+00
    f = 0.539412243677190440263092985511D+00
    g = 0.230293878161404779868453507244D+00

    w1 = 0.051739766065744133555179145422D+00
    w2 = 0.008007799555564801597804123460D+00
    w3 = 0.046868898981821644823226732071D+00
    w4 = 0.046590940183976487960361770070D+00
    w5 = 0.031016943313796381407646220131D+00
    w6 = 0.010791612736631273623178240136D+00
    w7 = 0.032195534242431618819414482205D+00
    w8 = 0.015445834210701583817692900053D+00
    w9 = 0.017822989923178661888748319485D+00
    wx = 0.037038683681384627918546472190D+00

    xtab(1:10) =   (/ a, b, c, c, d, e, e, f, g, g /)
    ytab(1:10) =   (/ a, c, b, c, e, d, e, g, f, g /)
    weight(1:37) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
                      w6, w6, w6, w7, w7, w7, w8, w8, w8, w8, w8, w8, w9, &
                      w9, w9, w9, w9, w9, wx, wx, wx, wx, wx, wx /)

    a = 0.772160036676532561750285570113D+00
    b = 0.113919981661733719124857214943D+00

    xtab(11) = a
    ytab(11) = b

    xtab(12) = b
    ytab(12) = a

    xtab(13) = b
    ytab(13) = b

    a = 0.009085399949835353883572964740D+00
    b = 0.495457300025082323058213517632D+00

    xtab(14) = a
    ytab(14) = b

    xtab(15) = b
    ytab(15) = a

    xtab(16) = b
    ytab(16) = b

    a = 0.062277290305886993497083640527D+00
    b = 0.468861354847056503251458179727D+00

    xtab(17) = a
    ytab(17) = b

    xtab(18) = b
    ytab(18) = a

    xtab(19) = b
    ytab(19) = b

    a = 0.022076289653624405142446876931D+00
    b = 0.851306504174348550389457672223D+00
    c = 0.126617206172027096933163647918263D+00

    xtab(20) = a
    ytab(20) = b

    xtab(21) = a
    ytab(21) = c

    xtab(22) = b
    ytab(22) = a

    xtab(23) = b
    ytab(23) = c

    xtab(24) = c
    ytab(24) = a

    xtab(25) = c
    ytab(25) = b

    a = 0.018620522802520968955913511549D+00
    b = 0.689441970728591295496647976487D+00
    c = 0.291937506468887771754472382212953D+00

    xtab(26) = a
    ytab(26) = b

    xtab(27) = a
    ytab(27) = c

    xtab(28) = b
    ytab(28) = a

    xtab(29) = b
    ytab(29) = c

    xtab(30) = c
    ytab(30) = a

    xtab(31) = c
    ytab(31) = b

    a = 0.096506481292159228736516560903D+00
    b = 0.635867859433872768286976979827D+00
    c = 0.267625659273967961282458816185681D+00

    xtab(32) = a
    ytab(32) = b

    xtab(33) = a
    ytab(33) = c

    xtab(34) = b
    ytab(34) = a

    xtab(35) = b
    ytab(35) = c

    xtab(36) = c
    ytab(36) = a

    xtab(37) = c
    ytab(37) = b

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_UNIT_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
    stop

  end if

  return
end
function triangle_unit_size ( rule )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_SIZE returns the "size" of a unit triangle quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jarle Berntsen, Terje Espelid,
!    Algorithm 706,
!    DCUTRI: an algorithm for adaptive cubature over a collection of triangles, 
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, September 1992, pages 329-342.
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!    DP Laurie,
!    Algorithm 584,
!    CUBTRI, Automatic Cubature Over a Triangle,
!    ACM Transactions on Mathematical Software,
!    Volume 8, Number 2, 1982, pages 210-218.
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!    Hans Rudolf Schwarz,
!    Methode der Finiten Elemente,
!    Teubner Studienbuecher, 1980,
!    ISBN: 3-519-02349-0.
!
!    Gilbert Strang, George Fix,
!    An Analysis of the Finite Element Method,
!    Prentice Hall, 1973, page 184,
!    ISBN: 096140888X,
!    LC: TA335.S77.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!    Olgierd Zienkiewicz,
!    The Finite Element Method,
!    Sixth Edition,
!    Butterworth-Heinemann, 2005,
!    ISBN: 0750663200,
!    TA640.2.Z54
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!     1, ORDER =  1, precision 1, Zienkiewicz #1.
!     2, ORDER =  2, precision 1, (the "vertex rule").
!     3, ORDER =  3, precision 2, Strang and Fix formula #1.
!     4, ORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
!     5, ORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
!     6, ORDER =  6, precision 3, Strang and Fix formula #4.
!     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
!     8, ORDER =  6, precision 4, Strang and Fix formula #5.
!     9, ORDER =  7, precision 4, Strang and Fix formula #6.
!    10, ORDER =  7, precision 5, Strang and Fix formula #7,
!        Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
!    11, ORDER =  9, precision 6, Strang and Fix formula #8.
!    12, ORDER = 12, precision 6, Strang and Fix formula #9.
!    13, ORDER = 13, precision 7, Strang and Fix formula #10.
!    14, ORDER =  7, precision ?.
!    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
!    16, ORDER = 64, precision 15, triangular product Gauss rule.
!    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
!    18, ORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
!    19, ORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
!    20, ORDER = 37, precision 13, from ACM TOMS #706.
!
!    Output, integer ( kind = 4 ) TRIANGLE_UNIT_SIZE, the order of the rule.
!
  implicit none

  integer ( kind = 4 ) rule
  integer ( kind = 4 ) triangle_unit_size

  if ( rule == 1 ) then
    triangle_unit_size = 1
  else if ( rule == 2 ) then
    triangle_unit_size = 3
  else if ( rule == 3 ) then
    triangle_unit_size = 3
  else if ( rule == 4 ) then
    triangle_unit_size = 3
  else if ( rule == 5 ) then
    triangle_unit_size = 4
  else if ( rule == 6 ) then
    triangle_unit_size = 6
  else if ( rule == 7 ) then
    triangle_unit_size = 6
  else if ( rule == 8 ) then
    triangle_unit_size = 6
  else if ( rule == 9 ) then
    triangle_unit_size = 7
  else if ( rule == 10 ) then
    triangle_unit_size = 7
  else if ( rule == 11 ) then
    triangle_unit_size = 9
  else if ( rule == 12 ) then
    triangle_unit_size = 12
  else if ( rule == 13 ) then
    triangle_unit_size = 13
  else if ( rule == 14 ) then
    triangle_unit_size = 7
  else if ( rule == 15 ) then
    triangle_unit_size = 16
  else if ( rule == 16 ) then
    triangle_unit_size = 64
  else if ( rule == 17 ) then
    triangle_unit_size = 19
  else if ( rule == 18 ) then
    triangle_unit_size = 19
  else if ( rule == 19 ) then
    triangle_unit_size = 28
  else if ( rule == 20 ) then
    triangle_unit_size = 37
  else
    triangle_unit_size = -1
  end if

  return
end
subroutine triangle_unit_sum ( func, order, xtab, ytab, weight, result )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_SUM carries out a quadrature rule in the unit triangle.
!
!  Integration region:
!
!      0 <= X,
!    and
!      0 <= Y, 
!    and
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the user supplied
!    function of two variables which is to be integrated,
!    of the form:
!      function func ( x, y )
!      real ( kind = 8 ) func
!      real ( kind = 8 ) x
!      real ( kind = 8 ) y
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas.
!
!    Input, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral of the function.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) result
  real ( kind = 8 ) triangle_unit_volume
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) ytab(order)

  quad = 0.0D+00

  do i = 1, order
    quad = quad + weight(i) * func ( xtab(i), ytab(i) )
  end do

  volume = triangle_unit_volume ( )
  result = quad * volume

  return
end
function triangle_unit_volume ( )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_VOLUME returns the "volume" of the unit triangle in 2D.
!
!  Integration region:
!
!      0 <= X,
!    and
!      0 <= Y, 
!    and
!      X + Y <= 1.
!
!  Discussion:
!
!    The "volume" of a triangle is usually called its area.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) TRIANGLE_UNIT_VOLUME, the volume of the unit
!    triangle.
!
  implicit none

  real ( kind = 8 ) triangle_unit_volume

  triangle_unit_volume = 1.0D+00 / 2.0D+00

  return
end
function triangle_volume ( x, y )

!*****************************************************************************80
!
!! TRIANGLE_VOLUME returns the "volume" of a triangle in 2D.
!
!  Integration region:
!
!      0 <= X,
!    and
!      0 <= Y, 
!    and
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(3), Y(3), the vertices of the triangle.
!
!    Output, real ( kind = 8 ) TRIANGLE_VOLUME, the volume of the triangle.
!
  implicit none

  real ( kind = 8 ) triangle_volume
  real ( kind = 8 ) x(3)
  real ( kind = 8 ) y(3)

  triangle_volume = 0.5D+00 * abs ( &
    x(1) * ( y(2) - y(3) ) + &
    x(2) * ( y(3) - y(1) ) + &
    x(3) * ( y(1) - y(2) ) )

  return
end
subroutine tvec_even ( nt, t )

!*****************************************************************************80
!
!! TVEC_EVEN computes an evenly spaced set of angles between 0 and 2*PI.
!
!  Discussion:
!
!    The computation realizes that 0 = 2 * PI, and does not include that value.
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
!    14 December 2000
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
    t(i) = real ( 2 * ( i - 1 ), kind = 8 ) * pi / real ( nt, kind = 8 )
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
!    15 December 2000
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
    t(i) = real ( 2 * i - 1, kind = 8 ) * pi / real ( nt, kind = 8 )
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
!    13 January 2001
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
      t(i) = real ( 2 * ( i - 1 ), kind = 8 ) * pi / real ( nt - 1, kind = 8 )
    end do
  end if

  return
end
subroutine tvec_even_bracket ( nt, theta1, theta2, t )

!*****************************************************************************80
!
!! TVEC_EVEN_BRACKET computes evenly spaced angles between THETA1 and THETA2.
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
!    20 May 2004
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
      t(i) = ( real ( nt - i,     kind = 8     ) * theta1   &
             + real (      i - 1, kind = 8     ) * theta2 ) &
             / real ( nt     - 1, kind = 8     )
    end do
  end if

  return
end
subroutine tvec_even_bracket2 ( nt, theta1, theta2, t )

!*****************************************************************************80
!
!! TVEC_EVEN_BRACKET2 computes evenly spaced angles between THETA1 and THETA2.
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
!    13 January 2001
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
           / real ( nt + 1,     kind = 8     )
  end do

  return
end
subroutine tvec_even_bracket3 ( nt, theta1, theta2, t )

!*****************************************************************************80
!
!! TVEC_EVEN_BRACKET3 computes evenly spaced angles between THETA1 and THETA2.
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
!    13 January 2001
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
subroutine vec_lex_next ( dim_num, base, a, more )

!*****************************************************************************80
!
!! VEC_LEX_NEXT generates vectors in lex order.
!
!  Discussion:
!
!    The vectors are produced in lexical order, starting with
!    (0,0,...,0),
!    (0,0,...,1), 
!    ...
!    (BASE-1,BASE-1,...,BASE-1).
!
!  Examples:
!
!    DIM_NUM = 2, 
!    BASE = 3
!
!    0   0
!    0   1
!    0   2
!    1   0
!    1   1
!    1   2
!    2   0
!    2   1
!    2   2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 May 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the size of the vectors to be used.
!
!    Input, integer ( kind = 4 ) BASE, the base to be used.  BASE = 2 will
!    give vectors of 0's and 1's, for instance.
!
!    Input/output, integer ( kind = 4 ) A(DIM_NUM), the next vector.
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  On return, MORE is TRUE if another vector has
!    been computed.  If MORE is returned FALSE, ignore the output 
!    vector and stop calling the routine.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) a(dim_num)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  logical more

  if ( .not. more ) then

    a(1:dim_num) = 0
    more = .true.

  else

    do i = dim_num, 1, -1

      a(i) = a(i) + 1

      if ( a(i) < base ) then
        return
      end if

      a(i) = 0

    end do

    more = .false.

  end if

  return
end
