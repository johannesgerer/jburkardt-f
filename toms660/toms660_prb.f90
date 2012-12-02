program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS660_PRB.
!
!  Discussion:
!
!    TOMS660_PRB is a test for TOMS660.
!
!    A quadratic function is sampled on a 6x6 regular grid of [0,1]x[0,1].
!
!    The interpolant is formed from this data.
!
!    The interpolant is evaluated on a 10x10 regular grid of [0,1]x[0,1]
!    and compared to the original function.
!
!    Since the function and interpolant are quadratic, they should be
!    numerically equal for this test.
!
!  Modified:
!
!    26 January 2012
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 36
  integer ( kind = 4 ), parameter :: nq = 13
  integer ( kind = 4 ), parameter :: nr = 3
  integer ( kind = 4 ), parameter :: nw = 19

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) eps
  real ( kind = 8 ) eq
  real ( kind = 8 ) eqx
  real ( kind = 8 ) eqy
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fq
  real ( kind = 8 ) fx
  real ( kind = 8 ) fy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lcell(3,3)
  integer ( kind = 4 ) lnext(n)
  real ( kind = 8 ) p(10)
  real ( kind = 8 ) px
  real ( kind = 8 ) py
  real ( kind = 8 ) q
  real ( kind = 8 ) q1
  real ( kind = 8 ) qs2val
  real ( kind = 8 ) qx
  real ( kind = 8 ) qy
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rq
  real ( kind = 8 ) rsq(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xx
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yy
!
!  Quadratic test function and partial derivatives.
!
  fq(xx,yy) = (         ( xx + 2.0D+00 * yy ) / 3.0E+00 )**2
  fx(xx,yy) = 2.0D+00 * ( xx + 2.0E+00 * yy ) / 9.0E+00
  fy(xx,yy) = 4.0D+00 * ( xx + 2.0E+00 * yy ) / 9.0E+00

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS660_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS660 library.'
!
!  Generate a 6 by 6 grid of nodes in the unit square with the natural ordering.
!
  k = 0
  do j = 5, 0, -1
    do i = 0, 5
      k = k + 1
      x(k) = real ( i, kind = 8 ) / 5.0D+00
      y(k) = real ( j, kind = 8 ) / 5.0D+00
    end do
  end do
!
!  Compute the data values.
!
  do k = 1, n
    f(k) = fq ( x(k), y(k) )
  end do
!
!  Call QSHEP2 to define the interpolant Q to this data.
!
  call qshep2 ( n, x, y, f, nq, nw, nr, lcell, lnext, xmin, ymin, &
    dx, dy, rmax, rsq, a, ier )

  if ( ier /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TOMS660_PRB - Error!'
    write ( *, '(a,i8)' ) '  Error in TOMS660, IER = ', ier
    stop
  end if
!
!  Generate a 10 by 10 uniform grid of interpolation points
!  (p(i),p(j)) in the unit square.
!
  do i = 1, 10
    p(i) = real ( i - 1, kind = 8 ) / 9.0D+00
  end do
!
!  Compute the machine precision EPS.
!
  eps = epsilon ( eps )
!
!  Compute the interpolation errors.
!
  eq = 0.0D+00
  eqx = 0.0D+00
  eqy = 0.0D+00

  do j = 1, 10

    py = p(j)

    do i = 1, 10

      px = p(i)

      q1 = qs2val ( px, py, n, x, y, f, nr, lcell, lnext, xmin, &
        ymin, dx, dy, rmax, rsq, a )

      call qs2grd ( px, py, n, x, y, f, nr, lcell, lnext, xmin, &
        ymin, dx, dy, rmax, rsq, a, q, qx, qy, ier )

      if ( ier /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TOMS660_PRB - Error!'
        write ( *, '(a,i6)' ) '  Error in QS2GRD, IER = ', ier
        stop
      end if

      if ( 3.0D+00 * abs ( q ) * eps < abs ( q1 - q ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TOMS660_PRB - Error!'
        write ( *, '(a)' ) '  The interpolated values Q1 (QS2VAL)'
        write ( *, '(a)' ) '  and Q (QS2GRD) differ.'
        write ( *, '(a,g14.6)' ) '  Q1 = ', q1
        write ( *, '(a,g14.6)' ) '  Q  = ', q
        stop
      end if

      eq = max ( eq, abs ( fq ( px, py ) - q ) )
      eqx = max ( eqx, abs ( fx ( px, py ) - qx ) )
      eqy = max ( eqy, abs ( fy ( px, py ) - qy ) )

    end do

  end do
!
!  Print the maximum errors and the ratio EQ / EPS.
!
  rq = eq / eps

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Maximum absolute errors in the interpolant Q and'
  write ( *, '(a)' ) '  partial derivatives QX and QY relative to machine'
  write ( *, '(a)' ) '  precision EPS.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Function   Max error   Max error/EPS'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a4,7x,e9.3,5x,f4.2)' ) 'Q   ', eq, rq
  write ( *, '(2x,a4,7x,e9.3)'      ) 'dQdX', eqx
  write ( *, '(2x,a4,7x,e9.3)'      ) 'dQdY', eqy
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS660_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
