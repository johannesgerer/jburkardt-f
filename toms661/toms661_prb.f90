program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS661_PRB.
!
!  Discussion:
!
!    TOMS661_PRB tests TOMS661.
!
!    TOMS661 is algorithm 661, collected algorithms from acm.
!    This work is published in the
!    ACM Transactions on Mathematical Software,
!    Volume 14, Number 2, page 151.
!
!    This program tests the scattered data interpolation
!    package TOMS661 by printing the maximum errors associated
!    with interpolated values and gradients on a 5 by 5 by 5
!    uniform grid in the unit cube.  The data set consists
!    of 64 nodes with data values taken from a quadratic
!    function for which the method is exact.  The ratio of maximum
!    interpolation error relative to the machine precision is
!    also printed.  This should be O(1).  The interpolated
!    values from QS3VAL and QS3GRD are compared for agreement.
!
!  Modified:
!
!    11 February 2007
!
  implicit none

  integer, parameter :: n = 64
  integer, parameter :: nq = 17
  integer, parameter :: nr = 3
  integer, parameter :: nw = 32

  real ( kind = 8 ) a(9,n)
  real ( kind = 8 ) eps
  real ( kind = 8 ) eq
  real ( kind = 8 ) eqx
  real ( kind = 8 ) eqy
  real ( kind = 8 ) eqz
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fq
  real ( kind = 8 ) fx
  real ( kind = 8 ) fy
  real ( kind = 8 ) fz
  integer i
  integer ier
  integer j
  integer k
  integer l
  integer lcell(3,3,3)
  integer lnext(n)
  real ( kind = 8 ) p(5)
  real ( kind = 8 ) px
  real ( kind = 8 ) py
  real ( kind = 8 ) pz
  real ( kind = 8 ) q
  real ( kind = 8 ) q1
  real ( kind = 8 ) qs3val
  real ( kind = 8 ) qx
  real ( kind = 8 ) qy
  real ( kind = 8 ) qz
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rq
  real ( kind = 8 ) rsq(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xx
  real ( kind = 8 ) xyzdel(3)
  real ( kind = 8 ) xyzmin(3)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yl
  real ( kind = 8 ) yy
  real ( kind = 8 ) z(n)
  real ( kind = 8 ) zl
  real ( kind = 8 ) zz
!
!  The quadratic test function and its partial derivatives.
!
  fq(xx,yy,zz) = ( ( xx + 2.0D+00 * yy + 3.0D+00 * zz ) / 6.0D+00 )**2
  fx(xx,yy,zz) = ( xx + 2.0D+00 * yy + 3.0D+00 * zz ) / 18.0D+00
  fy(xx,yy,zz) = ( xx + 2.0D+00 * yy + 3.0D+00 * zz ) / 9.0D+00
  fz(xx,yy,zz) = ( xx + 2.0D+00 * yy + 3.0D+00 * zz ) / 6.0D+00

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS661_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS661 library.'
!
!  Generate a 4 by 4 by 4 grid of nodes in the unit cube.
!
  l = 0
  do k = 1, 4
    zl = real ( k - 1, kind = 8 ) / 3.0D+00
    do j = 1, 4
      yl = real ( j - 1, kind = 8 ) / 3.0D+00
      do i = 1, 4
        l = l + 1
        x(l) = real ( i - 1, kind = 8 ) / 3.0D+00
        y(l) = yl
        z(l) = zl
      end do
    end do
  end do
!
!  Compute the data values.
!
  do l = 1, n
    f(l) = fq ( x(l), y(l), z(l) )
  end do
!
!  Compute parameters defining the interpolant Q.
!
  call qshep3 ( n, x, y, z, f, nq, nw, nr, lcell, lnext, xyzmin, &
    xyzdel, rmax, rsq, a, ier )

  if ( ier /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TOMS661_PRB'
    write ( *, '(a,i6)' ) '  Error return from QSHEP3, IER = ', ier
    stop
  end if
!
!  Generate a 5 by 5 by 5 uniform grid of interpolation
!  points (p(i),p(j),p(k)) in the unit cube.  The eight
!  corners coincide with nodes.
!
  do i = 1, 5
    p(i) = real ( i - 1, kind = 8 ) / 4.0D+00
  end do
!
!  Compute the machine precision.
!
  eps = epsilon ( eps )
!
!  Compute interpolation errors and test for agreement in the
!  Q values returned by qs3val and qs3grd.
!
  eq = 0.0D+00
  eqx = 0.0D+00
  eqy = 0.0D+00
  eqz = 0.0D+00

  do k = 1, 5
    pz = p(k)
    do j = 1, 5
      py = p(j)
      do i = 1, 5
        px = p(i)

        q1 = qs3val ( px, py, pz, n, x, y, z, f, nr, lcell, lnext, &
          xyzmin, xyzdel, rmax, rsq, a )

        call qs3grd ( px, py, pz, n, x, y, z, f, nr, lcell, lnext, &
          xyzmin, xyzdel, rmax, rsq, a, q, qx, qy, qz, ier )

        if ( ier /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TOMS661_PRB'
          write ( *, '(a,i6)' ) '  Error return from QS3GRD, IER = ', ier
          stop
        end if

        if ( abs ( q1 - q ) > 3.0D+00 * abs ( q ) * eps ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TOMS661_PRB - Error.'
          write ( *, '(a)' ) '  The interpolated values Q1 (QS3VAL) and'
          write ( *, '(a)' ) '  Q (QS3GRD) differ.'
          write ( *, '(a,g14.6)' ) '  Q1 = ', q1
          write ( *, '(a,g14.6)' ) '  Q  = ', q
          stop
        end if

        eq  = max ( eq,  abs ( fq(px,py,pz) - q  ) )
        eqx = max ( eqx, abs ( fx(px,py,pz) - qx ) )
        eqy = max ( eqy, abs ( fy(px,py,pz) - qy ) )
        eqz = max ( eqz, abs ( fz(px,py,pz) - qz ) )

      end do
    end do
  end do
!
!  Print errors and the ratio eq/eps.
!
  rq = eq / eps

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Maximum absolute errors in the interpolant Q '
  write ( *, '(a)' ) '  and partial derivatives (Qx,Qy,Qz) relative '
  write ( *, '(a)' ) '  to machine precision EPS.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Function   Max error  Max error/EPS'
  write ( *, '(a)' ) ' '
  write ( *, '(a,e9.3,a,f4.2)' ) '      Q       ', eq, '       ', rq
  write ( *, '(a,e9.3)'        ) '      Qx      ', eqx
  write ( *, '(a,e9.3)'        ) '      Qy      ', eqy
  write ( *, '(a,e9.3)'        ) '      Qz      ', eqz
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS661_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
