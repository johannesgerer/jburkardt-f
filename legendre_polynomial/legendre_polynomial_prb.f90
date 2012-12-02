program main

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_PRB tests the LEGENDRE_POLYNOMIAL library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) p

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the LEGENDRE_POLYNOMIAL library.'

  call legendre_polynomial_test01 ( )
  call legendre_polynomial_test02 ( )
  call legendre_polynomial_test03 ( )
  call legendre_polynomial_test04 ( )

  p = 5
  b = 0.0D+00
  call legendre_polynomial_test05 ( p, b )

  p = 5
  b = 1.0D+00
  call legendre_polynomial_test05 ( p, b )

  p = 5
  e = 0
  call legendre_polynomial_test06 ( p, e )

  p = 5
  e = 1
  call legendre_polynomial_test06 ( p, e )

  call legendre_polynomial_test07 ( )
  call legendre_polynomial_test08 ( )
  call legendre_polynomial_test09 ( )

  p = 5
  call legendre_polynomial_test10 ( p )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine legendre_polynomial_test01 ( )

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_TEST01 tests P_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) e
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ), parameter :: m = 1
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: v(:,:)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_TEST01:'
  write ( *, '(a)' ) '  P_POLYNOMIAL_VALUES stores values of'
  write ( *, '(a)' ) '  the Legendre polynomial P(n,x).'
  write ( *, '(a)' ) '  P_POLYNOMIAL evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                        Tabulated                 Computed'
  write ( *, '(a)' ) '     N        X           P(N,X)                    P(N,X)                     Error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call p_polynomial_values ( n_data, n, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( v(m,n+1) )
    call p_polynomial ( m, n, x, v )
    fx2 = v(1,n+1)
    deallocate ( v )

    e = fx1 - fx2

    write ( *, '(2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, x, fx1, fx2, e

  end do

  return
end
subroutine legendre_polynomial_test02 ( )

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_TEST02 tests P_POLYNOMIAL_COEFFICIENTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_TEST02'
  write ( *, '(a)' ) &
    '  P_POLYNOMIAL_COEFFICIENTS determines polynomial coefficients of P(n,x).'

  call p_polynomial_coefficients ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  P(', i, ',x) = '
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( c(i,j) == 0.0D+00 ) then

      else if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x^', j
      end if
    end do
  end do
 
  return
end
subroutine legendre_polynomial_test03 ( )

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_TEST03 tests P_POLYNOMIAL_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) degree
  real ( kind = 8 ), allocatable :: lz(:,:)
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable :: z(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_TEST03:'
  write ( *, '(a)' ) '  P_POLYNOMIAL_ZEROS computes the zeros of P(n,x)'
  write ( *, '(a)' ) '  Check by calling L_POLYNOMIAL there.'

  do degree = 1, 5

    allocate ( z(1:degree) )
    call p_polynomial_zeros ( degree, z )
    write ( title, '(a,i1,a)' ) '  Computed zeros for P(', degree, ',z):'
    call r8vec_print ( degree, z, title )

    allocate ( lz(degree,0:degree) )
    call p_polynomial ( degree, degree, z, lz )
    write ( title, '(a,i1,a)' ) '  Evaluate P(', degree, ',z):'
    call r8vec_print ( degree, lz(1:degree,degree), title )

    deallocate ( lz )
    deallocate ( z )

  end do

  return
end
subroutine legendre_polynomial_test04 ( )

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_TEST04 tests P_QUADRATURE_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e
  real ( kind = 8 ), allocatable :: f(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) q
  real ( kind = 8 ) q_exact
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_TEST04:'
  write ( *, '(a)' ) '  P_QUADRATURE_RULE computes the quadrature rule'
  write ( *, '(a)' ) '  associated with P(n,x)'

  n = 5
  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call p_quadrature_rule ( n, x, w )

  call r8vec2_print ( n, x, w,  '      X            W' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use the quadrature rule to estimate:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Q = Integral ( -1 <= X < +1 ) X^E dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   E       Q_Estimate      Q_Exact'
  write ( *, '(a)' ) ' '

  allocate ( f(1:n) )

  do e = 0, 2 * n - 1
    if ( e == 0 ) then
      f(1:n) = 1.0D+00
    else
      f(1:n) = x(1:n)**e
    end if
    q = dot_product ( w(1:n), f(1:n) )
    call p_integral ( e, q_exact )
    write ( *, '(2x,i2,2x,g14.6,2x,g14.6)' ) e, q, q_exact
  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine legendre_polynomial_test05 ( p, b )

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_TEST05 tests P_EXPONENTIAL_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polynomial 
!    factors.
!
!    Input, real ( kind = 8 ) B, the coefficient of X in the exponential factor.
!
  implicit none

  real ( kind = 8 ) b
  integer ( kind = 4 ) p
  real ( kind = 8 ), allocatable :: table(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_TEST05'
  write ( *, '(a)' ) '  Compute an exponential product table for P(n,x):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -1 <= x <= +1 ) exp(b*x) P(i,x) P(j,x) dx'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponential argument coefficient B = ', b

  allocate ( table(0:p,0:p) )

  call p_exponential_product ( p, b, table )

  call r8mat_print ( p + 1, p + 1, table, '  Exponential product table:' )

  deallocate ( table )

  return
end
subroutine legendre_polynomial_test06 ( p, e )

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_TEST06 tests P_POWER_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polynomial 
!    factors.
!
!    Input, integer ( kind = 4 ) E, the exponent of X.
!
  implicit none

  integer ( kind = 4 ) e
  integer ( kind = 4 ) p
  real ( kind = 8 ), allocatable :: table(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_TEST06'
  write ( *, '(a)' ) '  Compute a power product table for P(n,x):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -1 <= x <= +1 ) x^e P(i,x) P(j,x) dx'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponent of X, E = ', e

  allocate ( table(0:p,0:p) )

  call p_power_product ( p, e, table )

  call r8mat_print ( p + 1, p + 1, table, '  Power product table:' )

  deallocate ( table )

  return
end
subroutine legendre_polynomial_test07 ( )

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_TEST07 tests PM_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) e
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: mm = 1
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: v(:,:)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_TEST07:'
  write ( *, '(a)' ) '  PM_POLYNOMIAL_VALUES stores values of'
  write ( *, '(a)' ) '  the Legendre polynomial Pm(n,m,x).'
  write ( *, '(a)' ) '  PM_POLYNOMIAL evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                                Tabulated                 Computed'
  write ( *, '(a)' ) '     N     M        X           Pm(N,M,X)                 Pm(N,M,X)             Error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call pm_polynomial_values ( n_data, n, m, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( v(mm,n+1) )
    call pm_polynomial ( mm, n, m, x, v )
    fx2 = v(1,n+1)
    deallocate ( v )

    e = fx1 - fx2

    write ( *, '(2x,i4,2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, m, x, fx1, fx2, e

  end do

  return
end
subroutine legendre_polynomial_test08 ( )

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_TEST08 tests PMN_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) e
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: mm = 1
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: v(:,:)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_TEST08:'
  write ( *, '(a)' ) '  PMN_POLYNOMIAL_VALUES stores values of'
  write ( *, '(a)' ) '  the Legendre polynomial Pmn(n,m,x).'
  write ( *, '(a)' ) '  PMN_POLYNOMIAL evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                                Tabulated                 Computed'
  write ( *, '(a)' ) '     N     M        X           Pmn(N,M,X)                Pmn(N,M,X)             Error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call pmn_polynomial_values ( n_data, n, m, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( v(mm,n+1) )
    call pmn_polynomial ( mm, n, m, x, v )
    fx2 = v(1,n+1)
    deallocate ( v )

    e = fx1 - fx2

    write ( *, '(2x,i4,2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, m, x, fx1, fx2, e

  end do

  return
end
subroutine legendre_polynomial_test09 ( )

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_TEST09 tests PMNS_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) e
  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: mm = 1
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: v(:,:)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_TEST09:'
  write ( *, '(a)' ) '  PMNS_POLYNOMIAL_VALUES stores values of'
  write ( *, '(a)' ) '  the Legendre polynomial Pmns(n,m,x).'
  write ( *, '(a)' ) '  PMNS_POLYNOMIAL evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                                Tabulated                 Computed'
  write ( *, '(a)' ) '     N     M        X           Pmns(N,M,X)                Pmns(N,M,X)             Error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call pmns_polynomial_values ( n_data, n, m, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( v(mm,n+1) )
    call pmns_polynomial ( mm, n, m, x, v )
    fx2 = v(1,n+1)
    deallocate ( v )

    e = fx1 - fx2

    write ( *, '(2x,i4,2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, m, x, fx1, fx2, e

  end do

  return
end
subroutine legendre_polynomial_test10 ( p )

!*****************************************************************************80
!
!! LEGENDRE_POLYNOMIAL_TEST10 tests PN_PAIR_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polynomial 
!    factors.
!
  implicit none

  integer ( kind = 4 ) p
  real ( kind = 8 ), allocatable :: table(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_TEST10'
  write ( *, '(a)' ) '  Compute a pair product table for Pn(n,x):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -1 <= x <= +1 ) Pn(i,x) Pn(j,x) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Pn(n,x) polynomials are orthonormal,'
  write ( *, '(a)' ) '  so T should be the identity matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p

  allocate ( table(0:p,0:p) )

  call pn_pair_product ( p, table )

  call r8mat_print ( p + 1, p + 1, table, '  Pair product table:' )

  deallocate ( table )

  return
end
