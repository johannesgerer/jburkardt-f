program main

!*****************************************************************************80
!
!! MAIN is the main program for HERMITE_PRB.
!
!  Discussion:
!
!    HERMITE_PRB tests HERMITE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the HERMITE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses f(x) = 1 + 2x + 3x^2 at x = 0, 1, 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), dimension ( n ) :: x =  (/ 0.0D+00, 1.0D+00,  2.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: y =  (/ 1.0D+00, 6.0D+00, 17.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: yp = (/ 2.0D+00, 8.0D+00, 14.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  HERMITE computes the Hermite interpolant to data.'
  write ( *, '(a)' ) '  Here, f(x) = 1 + 2x + 3x^2.'

  call hermite_demo ( n, x, y, yp )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 uses f(x) = 6 + 5x + 4x^2 + 3x^3 + 2x^4 + x^5 at x = 0, 1, 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yp(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  HERMITE computes the Hermite interpolant to data.'
  write ( *, '(a)' ) '  Here, f(x) = 6 + 5x + 4x^2 + 3x^3 + 2x^4 + x^5.'

  do i = 1, 3

    x(i) = real ( i - 1, kind = 8 )

    y(i) = 6.0D+00 + x(i) * ( &
           5.0D+00 + x(i) * ( &
           4.0D+00 + x(i) * ( &
           3.0D+00 + x(i) * ( &
           2.0D+00 + x(i) ) ) ) )

    yp(i) = 5.0D+00 + x(i) * ( &
            8.0D+00 + x(i) * ( &
            9.0D+00 + x(i) * ( &
            8.0D+00 + x(i) *   &
            5.0D+00 ) ) )
  end do

  call hermite_demo ( n, x, y, yp )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 uses f(x) = r1 + r2x + r3x^2 + r4x^3 + r5x^4 + r6x^5 at x = r7 r8 r9
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) c(0:2*n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yp(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  HERMITE computes the Hermite interpolant to data.'
  write ( *, '(a)' ) '  Here, f(x) is a fifth order polynomial with random'
  write ( *, '(a)' ) '  coefficients, and the abscissas are random.'

  seed = 123456789

  call r8vec_uniform_01 ( n, seed, x )
  call r8vec_print ( n, x, '  Random abscissas' )

  call r8vec_uniform_01 ( 2 * n, seed, c )
  call r8vec_print ( 2 * n, c, '  Random polynomial coefficients.' )

  do i = 1, 3

    y(i) = c(0) + x(i) * ( &
           c(1) + x(i) * ( &
           c(2) + x(i) * ( &
           c(3) + x(i) * ( &
           c(4) + x(i) * ( &
           c(5) ) ) ) ) )

    yp(i) = c(1)           + x(i) * ( &
            c(2) * 2.0D+00 + x(i) * ( &
            c(3) * 3.0D+00 + x(i) * ( &
            c(4) * 4.0D+00 + x(i) *   &
            c(5) * 5.0D+00 ) ) )
  end do

  call hermite_demo ( n, x, y, yp )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 interpolates the Runge function at equally spaced points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) max_dif
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ns
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xdp(:)
  real ( kind = 8 ), allocatable :: xs(:)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xt
  real ( kind = 8 ), allocatable :: y(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: ydp(:)
  real ( kind = 8 ), allocatable :: yp(:)
  real ( kind = 8 ), allocatable :: ys(:)
  real ( kind = 8 ) yt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  HERMITE computes the Hermite interpolant to data.'
  write ( *, '(a)' ) '  Here, f(x) is the Runge function'
  write ( *, '(a)' ) '  and the data is evaluated at equally spaced points.'
  write ( *, '(a)' ) '  As N increases, the maximum error grows.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     Max | F(X) - H(F(X)) |'
  write ( *, '(a)' ) ' '

  do n = 3, 15, 2

    allocate ( x(1:n) )
    allocate ( y(1:n) )
    allocate ( yp(1:n) )

    nd = 2 * n

    allocate ( xd(1:nd) )
    allocate ( yd(1:nd) )
    allocate ( xdp(1:nd-1) )
    allocate ( ydp(1:nd-1) )

    ns = 10 * ( n - 1 ) + 1
    allocate ( xs(1:ns) )
    allocate ( ys(1:ns) )

    xlo = -5.0D+00
    xhi = +5.0D+00
    call r8vec_linspace ( n, xlo, xhi, x )

    y(1:n) = 1.0D+00 / ( 1.0D+00 + x(1:n)**2 )
    yp(1:n) = - 2.0D+00 * x(1:n) / ( 1.0D+00 + x(1:n)**2 )**2

    call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
!
!  Compare exact and interpolant at sample points.
!
    call r8vec_linspace ( ns, xlo, xhi, xs )

    call dif_vals ( nd, xd, yd, ns, xs, ys )

    max_dif = 0.0D+00
    do i = 1, ns
      xt = xs(i)
      yt = 1.0D+00 / ( 1.0D+00 + xt * xt )
      max_dif = max ( max_dif, abs ( ys(i) - yt ) )
    end do

    write ( *, '(2x,i4,2x,g14.6)' ) n, max_dif

    deallocate ( x )
    deallocate ( xd )
    deallocate ( xdp )
    deallocate ( xs )
    deallocate ( y )
    deallocate ( yd )
    deallocate ( ydp )
    deallocate ( yp )
    deallocate ( ys )

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 interpolates the Runge function at Chebyshev points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) max_dif
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ns
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xdp(:)
  real ( kind = 8 ), allocatable :: xs(:)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xt
  real ( kind = 8 ), allocatable :: y(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: ydp(:)
  real ( kind = 8 ), allocatable :: yp(:)
  real ( kind = 8 ), allocatable :: ys(:)
  real ( kind = 8 ) yt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  HERMITE computes the Hermite interpolant to data.'
  write ( *, '(a)' ) '  Here, f(x) is the Runge function'
  write ( *, '(a)' ) '  and the data is evaluated at Chebyshev points.'
  write ( *, '(a)' ) '  As N increases, the maximum error goes down.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     Max | F(X) - H(F(X)) |'
  write ( *, '(a)' ) ' '

  do n = 3, 15, 2

    allocate ( x(1:n) )
    allocate ( y(1:n) )
    allocate ( yp(1:n) )

    nd = 2 * n

    allocate ( xd(1:nd) )
    allocate ( yd(1:nd) )

    allocate ( xdp(1:nd-1) )
    allocate ( ydp(1:nd-1) )

    ns = 10 * ( n - 1 ) + 1
    allocate ( xs(1:ns) )
    allocate ( ys(1:ns) )

    xlo = -5.0D+00
    xhi = +5.0D+00
    call r8vec_chebyshev ( n, xlo, xhi, x )

    y(1:n) = 1.0D+00 / ( 1.0D+00 + x(1:n)**2 )
    yp(1:n) = - 2.0D+00 * x(1:n) / ( 1.0D+00 + x(1:n)**2 )**2

    call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
!
!  Compare exact and interpolant at sample points.
!
    call r8vec_linspace ( ns, xlo, xhi, xs )

    call dif_vals ( nd, xd, yd, ns, xs, ys )

    max_dif = 0.0D+00
    do i = 1, ns
      xt = xs(i)
      yt = 1.0D+00 / ( 1.0D+00 + xt * xt )
      max_dif = max ( max_dif, abs ( ys(i) - yt ) )
    end do

    write ( *, '(2x,i4,2x,g14.6)' ) n, max_dif

    deallocate ( x )
    deallocate ( xd )
    deallocate ( xdp )
    deallocate ( xs )
    deallocate ( y )
    deallocate ( yd )
    deallocate ( ydp )
    deallocate ( yp )
    deallocate ( ys )

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests HERMITE_BASIS_0 and HERMITE_BASIS_1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nd = 2

  real ( kind = 8 ) f01
  real ( kind = 8 ) f02
  real ( kind = 8 ) f11
  real ( kind = 8 ) f12
  integer ( kind = 4 ) i
  real ( kind = 8 ), dimension ( nd ) :: xd = (/ 0.0D+00, 10.0D+00 /)
  real ( kind = 8 ) xv
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yh
  real ( kind = 8 ) ypd(nd)
  real ( kind = 8 ) yv

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  HERMITE_BASIS_0 and HERMITE_BASIS_1 evaluate the'
  write ( *, '(a)' ) '  Hermite global polynomial basis functions'
  write ( *, '(a)' ) '  of type 0: associated with function values, and'
  write ( *, '(a)' ) '  of type 1: associated with derivative values.'
!
!  Let y = x^3 + x^2 + x + 1,
!  and compute the Hermite global polynomial interpolant based on two 
!  abscissas:
!
  yd(1:nd) = xd(1:nd)**3 + xd(1:nd)**2 + xd(1:nd) + 1.0D+00
  ypd(1:nd) = 3.0D+00 * xd(1:nd)**2 + 2.0D+00 * xd(1:nd) + 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Interpolate y = x^3 + x^2 + x + 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     XD         Y(XD)      Y''(XD)'
  write ( *, '(a)' ) ' '
  do i = 1, nd
    write ( *, '(2x,g10.4,2x,g10.4,2x,g10.4)' ) xd(i), yd(i), ypd(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     XV         Y(XV)      H(XV)'
  write ( *, '(a)' ) ' '

  do i = 1, 11

    xv = real ( i - 1, kind = 8 )

    yv = xv**3 + xv**2 + xv + 1.0D+00

    call hermite_basis_0 ( 2, xd, 1, xv, f01 )
    call hermite_basis_1 ( 2, xd, 1, xv, f11 )
    call hermite_basis_0 ( 2, xd, 2, xv, f02 )
    call hermite_basis_1 ( 2, xd, 2, xv, f12 )

    yh = yd(1) * f01 + ypd(1) * f11 + yd(2) * f02 + ypd(2) * f12

    write ( *, '(2x,g10.4,2x,g10.4,2x,g10.4)' ), xv, yv, yh

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests HERMITE_INTERPOLANT_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) q
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  HERMITE_INTERPOLANT_RULE'
  write ( *, '(a)' ) '  is given a set of N abscissas for a Hermite interpolant'
  write ( *, '(a)' ) '  and returns N pairs of quadrature weights'
  write ( *, '(a)' ) '  for function and derivative values at the abscissas.'

  n = 3
  allocate ( x(1:n) )
  allocate ( w(1:2*n) )
  a = 0.0D+00
  b = 10.0D+00
  call r8vec_linspace ( n, a, b, x )
  call hermite_interpolant_rule ( n, a, b, x, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I       X               W(F(X))        W(F''(X))'
  write ( *, '(a)' ) ' '
  k = 1
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, x(i), w(k), w(k+1)
    k = k + 2
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Use the quadrature rule over interval ', a, ' to ' , b
  write ( *, '(a)' ) ' '

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * 1 + w(k+1) * 0.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of 1 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i) + w(k+1) * 1.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X^2 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * sin ( x(i) ) + w(k+1) * cos ( x(i) )
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of SIN(X) = ', q

  deallocate ( w )
  deallocate ( x )

  n = 3
  allocate ( x(1:n) )
  allocate ( w(1:2*n) )
  a = 0.0D+00
  b = 1.0D+00
  call r8vec_linspace ( n, a, b, x )
  call hermite_interpolant_rule ( n, a, b, x, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I       X               W(F(X))        W(F''(X))'
  write ( *, '(a)' ) ' '
  k = 1
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, x(i), w(k), w(k+1)
    k = k + 2
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Use the quadrature rule over interval ', a, ' to ' , b
  write ( *, '(a)' ) ' '

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * 1 + w(k+1) * 0.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of 1 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i) + w(k+1) * 1.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X^2 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * sin ( x(i) ) + w(k+1) * cos ( x(i) )
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of SIN(X) = ', q

  deallocate ( w )
  deallocate ( x )

  n = 11
  allocate ( x(1:n) )
  allocate ( w(1:2*n) )
  a = 0.0D+00
  b = 10.0D+00
  call r8vec_linspace ( n, a, b, x )
  call hermite_interpolant_rule ( n, a, b, x, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I       X               W(F(X))        W(F''(X))'
  write ( *, '(a)' ) ' '
  k = 1
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, x(i), w(k), w(k+1)
    k = k + 2
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Use the quadrature rule over interval ', a, ' to ' , b
  write ( *, '(a)' ) ' '

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * 1 + w(k+1) * 0.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of 1 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i) + w(k+1) * 1.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X^2 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * sin ( x(i) ) + w(k+1) * cos ( x(i) )
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of SIN(X) = ', q

  deallocate ( w )
  deallocate ( x )

  n = 11
  allocate ( x(1:n) )
  allocate ( w(1:2*n) )
  a = 0.0D+00
  b = 1.0D+00
  call r8vec_linspace ( n, a, b, x )
  call hermite_interpolant_rule ( n, a, b, x, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I       X               W(F(X))        W(F''(X))'
  write ( *, '(a)' ) ' '
  k = 1
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, x(i), w(k), w(k+1)
    k = k + 2
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Use the quadrature rule over interval ', a, ' to ' , b
  write ( *, '(a)' ) ' '

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * 1 + w(k+1) * 0.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of 1 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i) + w(k+1) * 1.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X^2 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * sin ( x(i) ) + w(k+1) * cos ( x(i) )
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of SIN(X) = ', q

  deallocate ( w )
  deallocate ( x )

  n = 11
  allocate ( x(1:n) )
  allocate ( w(1:2*n) )
  a = 0.0D+00
  b = 1.0D+00
  call r8vec_chebyshev ( n, a, b, x )
  call hermite_interpolant_rule ( n, a, b, x, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I       X               W(F(X))        W(F''(X))'
  write ( *, '(a)' ) ' '
  k = 1
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) i, x(i), w(k), w(k+1)
    k = k + 2
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Use the quadrature rule over interval ', a, ' to ' , b
  write ( *, '(a)' ) ' '

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * 1 + w(k+1) * 0.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of 1 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i) + w(k+1) * 1.0D+00
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of X^2 = ', q

  q = 0.0D+00
  k = 1
  do i = 1, n
    q = q + w(k) * sin ( x(i) ) + w(k+1) * cos ( x(i) )
    k = k + 2
  end do
  write ( *, '(a,g14.6)' ) '  Estimate integral of SIN(X) = ', q

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tabulates the interpolant and its derivative. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ndp
  integer ( kind = 4 ) ns
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xd(:)
  real ( kind = 8 ), allocatable :: xdp(:)
  real ( kind = 8 ), allocatable :: xs(:)
  real ( kind = 8 ), allocatable :: y(:)
  real ( kind = 8 ), allocatable :: yd(:)
  real ( kind = 8 ), allocatable :: ydp(:)
  real ( kind = 8 ), allocatable :: yp(:)
  real ( kind = 8 ), allocatable :: ys(:)
  real ( kind = 8 ), allocatable :: ysp(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  HERMITE_INTERPOLANT sets up the Hermite interpolant.'
  write ( *, '(a)' ) '  HERMITE_INTERPOLANT_VALUE evaluates it.'
  write ( *, '(a)' ) '  Consider data for y=sin(x) at x=0,1,2,3,4.'

  n = 5
  allocate ( x(1:n) )
  allocate ( y(1:n) )
  allocate ( yp(1:n) )

  nd = 2 * n
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )

  ndp = 2 * n - 1
  allocate ( xdp(1:ndp) )
  allocate ( ydp(1:ndp) )

  call r8vec_linspace ( n, 0.0D+00, 4.0D+00, x )
  y(1:n) = sin ( x(1:n) )
  yp(1:n) = cos ( x(1:n) ) 

  call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
!
!  Now sample the interpolant at NS points, which include data values.
!
  ns = 4 * ( n - 1 ) + 1
  allocate ( xs(1:ns) )
  allocate ( ys(1:ns) )
  allocate ( ysp(1:ns) )

  call r8vec_linspace ( ns, 0.0D+00, 4.0D+00, xs )

  call hermite_interpolant_value ( nd, xd, yd, xdp, ydp, ns, xs, ys, ysp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the following table, there should be perfect'
  write ( *, '(a)' ) '  agreement between F and H, and F'' and H'''
  write ( *, '(a)' ) '  at the data points X = 0, 1, 2, 3, and 4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In between, H and H'' approximate F and F''.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I       X(I)          F(X(I))         H(X(I))' // &
    '        F''(X(I))        H''(X(I))'
  write ( *, '(a)' ) ' '
  do i = 1, ns
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, xs(i), sin ( xs(i) ), ys(i), &
      cos ( xs(i) ), ysp(i)
  end do

  deallocate ( x )
  deallocate ( xd )
  deallocate ( xdp )
  deallocate ( xs )
  deallocate ( y )
  deallocate ( yd )
  deallocate ( ydp )
  deallocate ( yp )
  deallocate ( ys )
  deallocate ( ysp )

  return
end
