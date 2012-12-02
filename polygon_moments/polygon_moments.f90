subroutine moment ( n, x, y, p, q, nu_pq )

!*****************************************************************************80
!
!! MOMENT computes an unnormalized moment of a polygon.
!
!  Discussion:
!
!    Nu(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carsten Steger,
!    On the calculation of arbitrary moments of polygons,
!    Technical Report FGBV-96-05,
!    Forschungsgruppe Bildverstehen, Informatik IX,
!    Technische Universitaet Muenchen, October 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the vertex coordinates.
!
!    Input, integer ( kind = 4 ) P, Q, the indices of the moment.
!
!    Output, real ( kind = 8 ) NU_PQ, the unnormalized moment Nu(P,Q).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) nu_pq
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) s_pq
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xi
  real ( kind = 8 ) xj
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yi
  real ( kind = 8 ) yj

  nu_pq = 0.0D+00

  xj = x(n)
  yj = y(n)

  do i = 1, n

    xi = x(i)
    yi = y(i)

    s_pq = 0.0D+00
    do k = 0, p
      do l = 0, q
        s_pq = s_pq &
          + r8_choose ( k + l, l ) * r8_choose ( p + q - k - l, q - l ) &
          * xi ** k * xj ** ( p - k ) &
          * yi ** l * yj ** ( q - l )
      end do
    end do

    nu_pq = nu_pq + ( xj * yi - xi * yj ) * s_pq

    xj = xi
    yj = yi

  end do

  nu_pq = nu_pq / real ( p + q + 2, kind = 8 ) &
    / real ( p + q + 1, kind = 8 ) &
    / r8_choose ( p + q, p )

  return
end
subroutine moment_central ( n, x, y, p, q, mu_pq )

!*****************************************************************************80
!
!! MOMENT_CENTRAL computes central moments of a polygon.
!
!  Discussion:
!
!    The central moment Mu(P,Q) is defined by
!
!      Mu(P,Q) = Integral ( polygon ) (x-Alpha(1,0))^p (y-Alpha(0,1))^q dx dy
!              / Area ( polygon )
!
!    where 
!
!      Alpha(1,0) = Integral ( polygon ) x dx dy / Area ( polygon )
!      Alpha(0,1) = Integral ( polygon ) y dx dy / Area ( polygon )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carsten Steger,
!    On the calculation of arbitrary moments of polygons,
!    Technical Report FGBV-96-05,
!    Forschungsgruppe Bildverstehen, Informatik IX,
!    Technische Universitaet Muenchen, October 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the vertex coordinates.
!
!    Input, integer ( kind = 4 ) P, Q, the indices of the moment.
!
!    Output, real ( kind = 8 ) MU_PQ, the unnormalized moment Mu(P,Q).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha_01
  real ( kind = 8 ) alpha_10
  real ( kind = 8 ) alpha_ij
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mu_pq
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) r8_mop
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  call moment_normalized ( n, x, y, 1, 0, alpha_10 )
  call moment_normalized ( n, x, y, 0, 1, alpha_01 )

  mu_pq = 0.0D+00

  do i = 0, p
    do j = 0, q

      call moment_normalized ( n, x, y, i, j, alpha_ij )

      mu_pq = mu_pq + r8_mop ( p + q - i - j ) &
        * r8_choose ( p, i ) * r8_choose ( q, j ) &
        * alpha_10 ** ( p - i ) * alpha_01 ** ( q - j ) * alpha_ij

    end do
  end do

  return
end
subroutine moment_normalized ( n, x, y, p, q, alpha_pq )

!*****************************************************************************80
!
!! MOMENT_NORMALIZED computes a normalized moment of a polygon.
!
!  Discussion:
!
!    Alpha(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy / Area ( polygon )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carsten Steger,
!    On the calculation of arbitrary moments of polygons,
!    Technical Report FGBV-96-05,
!    Forschungsgruppe Bildverstehen, Informatik IX,
!    Technische Universitaet Muenchen, October 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of vertices of the polygon.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the vertex coordinates.
!
!    Input, integer ( kind = 4 ) P, Q, the indices of the moment.
!
!    Output, real ( kind = 8 ) ALPHA_PQ, the normalized moment Alpha(P,Q).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha_pq
  real ( kind = 8 ) nu_00
  real ( kind = 8 ) nu_pq
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  call moment ( n, x, y, p, q, nu_pq )
  call moment ( n, x, y, 0, 0, nu_00 )

  alpha_pq = nu_pq / nu_00

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
function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8.
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
