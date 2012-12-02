program main

!*****************************************************************************80
!
!! QW_GOLUB_WELSCH_PRB tests QW_GOLUB_WELSCH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QW_GOLUB_WELSCH_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the QW_GOLUB_WELSCH library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QW_GOLUB_WELSCH_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests QW_GOLUB_WELSCH for the Chebyshev Type 1 weight.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), allocatable :: aj(:)
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: bj(:)
  integer ( kind = 4 ) j
  real ( kind = 8 ) mu0
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  a = -1.0D+00
  b = +1.0D+00
  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Compute points and weights for Gauss quadrature'
  write ( *, '(a)' ) '  with the Chebyshev Type 1 weight w(x) = 1/sqrt(1-x^2).'
  write ( *, '(a,i4)' ) '  Order N = ', n
  write ( *, '(a,f10.4,a,f10.4,a)' ) '  Interval = [', a, ',', b, ']'

  allocate ( aj(1:n) )
  allocate ( bj(1:n) )

  aj(1:n) = 0.0D+00

  bj(1) = 1.0D+00 / 2.0D+00
  do j = 2, n - 1
    bj(j) = 1.0D+00 / 4.0D+00
  end do
  bj(n) = 0.0D+00

  bj(1:n) = sqrt ( bj(1:n) )
!
!  Note that Kautsky and Elhay seem to indicate mu0 = 1.0, but
!  PI seems correct to me.
!
  mu0 = pi

  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call sgqf ( n, aj, bj, mu0, x, w )

  call r8vec_print ( n, x, '  Abscissas:' )
  call r8vec_print ( n, w, '  Weights:' )

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests QW_GOLUB_WELSCH for the Chebyshev Type 2 weight.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), allocatable :: aj(:)
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: bj(:)
  integer ( kind = 4 ) j
  real ( kind = 8 ) mu0
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  a = -1.0D+00
  b = +1.0D+00
  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Compute points and weights for Gauss quadrature'
  write ( *, '(a)' ) '  with the Chebyshev Type 2 weight w(x) = sqrt(1-x^2).'
  write ( *, '(a,i4)' ) '  Order N = ', n
  write ( *, '(a,f10.4,a,f10.4,a)' ) '  Interval = [', a, ',', b, ']'

  allocate ( aj(1:n) )
  allocate ( bj(1:n) )

  aj(1:n) = 0.0D+00

  do j = 1, n - 1
    bj(j) = 1.0D+00 / 4.0D+00
  end do
  bj(n) = 0.0D+00

  bj(1:n) = sqrt ( bj(1:n) )
!
!  This value of MU0 seems correct to me.
!
  mu0 = pi / 2.0D+00

  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call sgqf ( n, aj, bj, mu0, x, w )

  call r8vec_print ( n, x, '  Abscissas:' )
  call r8vec_print ( n, w, '  Weights:' )

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests QW_GOLUB_WELSCH for the Gegenbauer weight.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), allocatable :: aj(:)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: bj(:)
  integer ( kind = 4 ) j
  real ( kind = 8 ) jr
  real ( kind = 8 ) mu0
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  a = -1.0D+00
  b = +1.0D+00
  alpha = 0.25D+00
  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  Compute points and weights for Gauss quadrature'
  write ( *, '(a)' ) '  with the Gegenbauer weight w(x) = (1-x^2)^alpha.'
  write ( *, '(a,i4)' ) '  Order N = ', n
  write ( *, '(a,f10.4)' ) '  ALPHA = ', alpha
  write ( *, '(a,f10.4,a,f10.4,a)' ) '  Interval = [', a, ',', b, ']'

  allocate ( aj(1:n) )
  allocate ( bj(1:n) )

  aj(1:n) = 0.0D+00

  do j = 1, n - 1
    jr = real ( j, kind = 8 )
    bj(j) = ( jr * ( 2.0D+00 * alpha + jr ) ) &
      / ( 4.0D+00 * ( alpha + jr )**2 - 1.0D+00 )
  end do
  bj(n) = 0.0D+00

  bj(1:n) = sqrt ( bj(1:n) )

  mu0 = gamma ( alpha + 1.0D+00 ) * gamma ( 0.5D+00 ) / gamma ( alpha + 1.5D+00 )

  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call sgqf ( n, aj, bj, mu0, x, w )

  call r8vec_print ( n, x, '  Abscissas:' )
  call r8vec_print ( n, w, '  Weights:' )

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests QW_GOLUB_WELSCH for the generalized Hermite weight.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: aj(:)
  real ( kind = 8 ) alpha
  real ( kind = 8 ), allocatable :: bj(:)
  integer ( kind = 4 ) j
  real ( kind = 8 ) mu0
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  n = 5
  alpha = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  Compute points and weights for Gauss quadrature'
  write ( *, '(a)' ) '  with the generalized Hermite weight w(x) = |x|^alpha * exp(-x^2).'
  write ( *, '(a,f10.4)' ) '  ALPHA = ', alpha
  write ( *, '(a,i4)' ) '  Order N = ', n
  write ( *, '(a,f10.4,a,f10.4,a)' ) '  Interval = (-oo,+oo)'

  allocate ( aj(1:n) )
  allocate ( bj(1:n) )

  aj(1:n) = 0.0D+00

  do j = 1, n - 1
    if ( mod ( j, 2 ) == 1 ) then
      bj(j) = ( real ( j, kind = 8 ) + alpha ) / 2.0D+00
    else
      bj(j) = real ( j, kind = 8 ) / 2.0D+00
    end if
  end do
  bj(n) = 0.0D+00

  bj(1:n) = sqrt ( bj(1:n) )

  mu0 = gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )

  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call sgqf ( n, aj, bj, mu0, x, w )

  call r8vec_print ( n, x, '  Abscissas:' )
  call r8vec_print ( n, w, '  Weights:' )

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests QW_GOLUB_WELSCH for the generalized Laguerre weight.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), allocatable :: aj(:)
  real ( kind = 8 ) alpha
  real ( kind = 8 ), allocatable :: bj(:)
  integer ( kind = 4 ) j
  real ( kind = 8 ) mu0
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  a = 0.0D+00
  alpha = 2.0D+00
  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  Compute points and weights for Gauss quadrature'
  write ( *, '(a)' ) '  with the generalized Laguerre weight w(x) = x^alpha * exp(-x).'
  write ( *, '(a,i4)' ) '  Order N = ', n
  write ( *, '(a,f10.4)' ) '  ALPHA = ', alpha
  write ( *, '(a,f10.4,a,f10.4,a)' ) '  Interval = [0,+oo)'

  allocate ( aj(1:n) )
  allocate ( bj(1:n) )

  do j = 1, n
    aj(j) = alpha + real ( 2 * j - 1, kind = 8 )
  end do

  do j = 1, n - 1
    bj(j) = real ( j, kind = 8 ) * ( alpha + real ( j, kind = 8 ) )
  end do
  bj(n) = 0.0D+00

  bj(1:n) = sqrt ( bj(1:n) )

  mu0 = gamma ( alpha + 1 )

  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call sgqf ( n, aj, bj, mu0, x, w )

  call r8vec_print ( n, x, '  Abscissas:' )
  call r8vec_print ( n, w, '  Weights:' )

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests QW_GOLUB_WELSCH for the Hermite weight.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: aj(:)
  real ( kind = 8 ), allocatable :: bj(:)
  integer ( kind = 4 ) j
  real ( kind = 8 ) mu0
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  Compute points and weights for Gauss quadrature'
  write ( *, '(a)' ) '  with the Hermite weight w(x) = exp(-x^2).'
  write ( *, '(a,i4)' ) '  Order N = ', n
  write ( *, '(a,f10.4,a,f10.4,a)' ) '  Interval = (-oo,+oo)'

  allocate ( aj(1:n) )
  allocate ( bj(1:n) )

  aj(1:n) = 0.0D+00

  do j = 1, n - 1
    bj(j) = real ( j, kind = 8 ) / 2.0D+00
  end do
  bj(n) = 0.0D+00

  bj(1:n) = sqrt ( bj(1:n) )

  mu0 = sqrt ( pi )

  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call sgqf ( n, aj, bj, mu0, x, w )

  call r8vec_print ( n, x, '  Abscissas:' )
  call r8vec_print ( n, w, '  Weights:' )

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests QW_GOLUB_WELSCH for the Jacobi weight.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), allocatable :: aj(:)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ), allocatable :: bj(:)
  integer ( kind = 4 ) j
  real ( kind = 8 ) jr
  real ( kind = 8 ) mu0
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  a = -1.0D+00
  b = +1.0D+00
  alpha = 0.25D+00
  beta = 0.75D+00
  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  Compute points and weights for Gauss quadrature'
  write ( *, '(a)' ) '  with the Jacobi weight w(x) = (1-x^2)^alpha*(1+x)^beta'
  write ( *, '(a,i4)' ) '  Order N = ', n
  write ( *, '(a,f10.4)' ) '  ALPHA = ', alpha
  write ( *, '(a,f10.4)' ) '  BETA =  ', beta
  write ( *, '(a,f10.4,a,f10.4,a)' ) '  Interval = [', a, ',', b, ']'

  allocate ( aj(1:n) )
  allocate ( bj(1:n) )

  do j = 1, n
    jr = real ( j, kind = 8 )
    aj(j) = ( beta - alpha ) * ( beta + alpha ) &
      / ( alpha + beta + 2.0D+00 * jr - 2.0D+00 ) &
      / ( alpha + beta + 2.0D+00 * jr )
  end do

  do j = 1, n - 1
    jr = real ( j, kind = 8 )
    bj(j) = 4.0D+00 * jr * ( alpha + jr ) * ( beta + jr ) &
      * ( alpha + beta + jr ) &
      / ( ( alpha + beta + 2.0D+00 * jr )**2 - 1.0D+00 ) &
      / ( alpha + beta + 2.0D+00 * jr )**2
  end do
  bj(n) = 0.0D+00

  bj(1:n) = sqrt ( bj(1:n) )

  mu0 = 2.0D+00**( alpha + beta + 1.0 ) &
    * gamma ( alpha + 1.0D+00 ) * gamma ( beta + 1.0D+00 ) &
    / gamma ( alpha + beta + 2.0D+00 )

  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call sgqf ( n, aj, bj, mu0, x, w )

  call r8vec_print ( n, x, '  Abscissas:' )
  call r8vec_print ( n, w, '  Weights:' )

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests QW_GOLUB_WELSCH for the Laguerre weight.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), allocatable :: aj(:)
  real ( kind = 8 ), allocatable :: bj(:)
  integer ( kind = 4 ) j
  real ( kind = 8 ) mu0
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  a = 0.0D+00
  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  Compute points and weights for Gauss quadrature'
  write ( *, '(a)' ) '  with the Laguerre weight w(x) = exp(-x).'
  write ( *, '(a,i4)' ) '  Order N = ', n
  write ( *, '(a,f10.4,a,f10.4,a)' ) '  Interval = [0,+oo)'

  allocate ( aj(1:n) )
  allocate ( bj(1:n) )

  do j = 1, n
    aj(j) = real ( 2 * j - 1, kind = 8 )
  end do

  do j = 1, n - 1
    bj(j) = real ( j**2, kind = 8 )
  end do
  bj(n) = 0.0D+00

  bj(1:n) = sqrt ( bj(1:n) )

  mu0 = 1.0D+00

  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call sgqf ( n, aj, bj, mu0, x, w )

  call r8vec_print ( n, x, '  Abscissas:' )
  call r8vec_print ( n, w, '  Weights:' )

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests QW_GOLUB_WELSCH for the Legendre weight.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ), allocatable :: aj(:)
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: bj(:)
  integer ( kind = 4 ) j
  real ( kind = 8 ) mu0
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  a = -1.0D+00
  b = +1.0D+00
  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09:'
  write ( *, '(a)' ) '  Compute points and weights for Gauss quadrature'
  write ( *, '(a)' ) '  with the Legendre weight w(x) = 1.'
  write ( *, '(a,i4)' ) '  Order N = ', n
  write ( *, '(a,f10.4,a,f10.4,a)' ) '  Interval = [', a, ',', b, ']'

  allocate ( aj(1:n) )
  allocate ( bj(1:n) )

  aj(1:n) = 0.0D+00

  do j = 1, n - 1
    bj(j) = real ( j**2, kind = 8 ) / real ( 4 * j**2 - 1, kind = 8 )
  end do
  bj(n) = 0.0D+00

  bj(1:n) = sqrt ( bj(1:n) )

  mu0 = 2.0D+00

  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call sgqf ( n, aj, bj, mu0, x, w )

  call r8vec_print ( n, x, '  Abscissas:' )
  call r8vec_print ( n, w, '  Weights:' )

  deallocate ( w )
  deallocate ( x )

  return
end
