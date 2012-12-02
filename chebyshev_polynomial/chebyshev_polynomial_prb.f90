program main

!*****************************************************************************80
!
!! MAIN is the main program for CHEBYSHEV_POLYNOMIAL_PRB.
!
!  Discussion:
!
!    CHEBYSHEV_POLYNOMIAL_PRB tests the CHEBYSHEV_POLYNOMIAL library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBYSHEV_POLYNOMIAL_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CHEBYSHEV_POLYNOMIAL library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBYSHEV_POLYNOMIAL_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests T_PROJECT_COEFFICIENTS_DATA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: c(:)
  real ( kind = 8 ), allocatable :: d(:)
  real ( kind = 8 ), allocatable :: d2(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBYSHEV_POLYNOMIAL_TEST01:'
  write ( *, '(a)' ) '  T_PROJECT_COEFFICIENTS_DATA estimates the Chebyshev polynomial'
  write ( *, '(a)' ) '  coefficients for a function given as data (x,fx).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we use fx = f(x) = x^2 for the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Since T(0,x) = 1 and T(2,x) = 2*x^2 - 1, the correct expansion is'
  write ( *, '(a)' ) '  f(x) = 1/2 T(0,x) + 0 T(1,x) + 1/2 T(2,x) + 0 * all other polys.'
!
!  Data in [0,1];
!
  a = 0.0D+00
  b = 1.0D+00
  m = 20
  seed = 123456789
  allocate ( x(1:m) )
  call r8vec_uniform_01 ( m, seed, x )
  allocate ( d(1:m) )
  d(1:m) = x(1:m)**2

  call r8vec2_print ( m, x, d, '  Data ( X, D ):' )

  n = 4
  allocate ( c(0:n) )
  call t_project_coefficients_data ( a, b, m, n, x, d, c )
  
  call r8vec_print ( n, c, '  Coefficients of Chebyshev expansion of degree 4.' )
!
!  Compare Chebyshev expansion and original function.
!
  allocate ( d2(1:m) )

  call t_project_value ( m, n, x, c, d2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I      X(I)     Data(I)      Chebyshev(X(I))'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(2x,i2,2x,g12.4,2x,g12.4,2x,g12.4)' ) i, x(i), d(i), d2(i)
  end do

  deallocate ( c )
  deallocate ( d )
  deallocate ( d2 )
  deallocate ( x )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests T_POLYNOMIAL_COEFFICIENTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  T_POLYNOMIAL_COEFFICIENTS determines coefficients for'
  write ( *, '(a)' ) '  Chebyshev polynomials T(n,x).'

  call t_polynomial_coefficients ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  T(', i, ',x) = '
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( c(i,j) == 0.0D+00 ) then

      else if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x ^ ', j
      end if
    end do
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests T_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ), allocatable :: fx2(:)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  T_POLYNOMIAL evaluates Chebyshev polynomials T(n,x).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                    Tabulated     Computed'
  write ( *, '(a)' ) '     N      X        T(n,x)        T(n,x)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call t_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( fx2(0:n) )

    call t_polynomial ( 1, n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

    deallocate ( fx2 )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests T_POLYNOMIAL_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: fx(:,:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: n_max = 5
  real ( kind = 8 ), allocatable :: z(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  T_POLYNOMIAL_ZEROS returns zeroes of T(n,x).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      X        T(n,x)'

  do n = 1, n_max

    allocate ( z(1:n) )

    call t_polynomial_zeros ( n, z )

    allocate ( fx(1:n,0:n) )

    call t_polynomial ( n, n, z, fx )

    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i8,2x,f8.4,2x,g14.6)' ) n, z(i), fx(i,n)
    end do

    deallocate ( fx )
    deallocate ( z )

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests T_QUADRATURE_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e
  real ( kind = 8 ), allocatable :: f(:)
  integer ( kind = 4 ) n
  real ( kind = 8 ) q
  real ( kind = 8 ) q_exact
  real ( kind = 8 ) t_integral
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  T_QUADRATURE_RULE computes the quadrature rule'
  write ( *, '(a)' ) '  associated with T(n,x);'

  n = 7
  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call t_quadrature_rule ( n, x, w )

  call r8vec2_print ( n, x, w, '      X            W' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use the quadrature rule to estimate:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Q = Integral ( -1 <= X <= +1 ) X^E / sqrt ( 1-x^2) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   E       Q_Estimate      Q_Exact'
  write ( *, '(a)' ) ' '

  allocate ( f(1:n) )

  do e = 0, 2 * n - 1
    if ( e == 0 ) then
      f(1:n) = 1.0D+00
    else
      f(1:n) = x(1:n) ** e
    end if
    q = dot_product ( w, f )
    q_exact = t_integral ( e )
    write ( *, '(2x,i2,g14.6,2x,g14.6)' ) e, q, q_exact
  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! CHEBYSHEV_POLYNOMIAL_TEST06 tests the projection of T(i,x) and T(j,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: c(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: phi(:,:)
  real ( kind = 8 ), allocatable :: phiw(:,:)
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  As a sanity check, make sure that the projection of:'
  write ( *, '(a)' ) '  T(i,x) onto T(j,x) is:'
  write ( *, '(a)' ) '  0 if i is not equal to j;'
  write ( *, '(a)' ) '  pi if i = j = 0;'
  write ( *, '(a)' ) '  pi/2 if i = j =/= 0.'

  n = 3

  allocate ( x(1:n+1) )
  allocate ( w(1:n+1) )

  call t_quadrature_rule ( n + 1, x, w )

  allocate ( c(0:n) )
  allocate ( phi(1:n+1,0:n) )
  allocate ( phiw(0:n,1:n+1) )

  call t_polynomial ( n + 1, n, x, phi )

  do j = 0, n

    do i = 1, n + 1
      phiw(0:n,i) = phi(i,0:n) * w(i)
    end do

    c(0:n) = matmul ( phiw(0:n,1:n+1), phi(1:n+1,j) )

    write ( title, '(a,i2,a)' ) &
      '  Chebyshev expansion coefficients for T(', j, ',x)'

    call r8vec_print ( n + 1, c, title )

  end do

  deallocate ( c )
  deallocate ( phi )
  deallocate ( phiw )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests T_PROJECT_COEFFICIENTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: c(:)
  intrinsic dexp
  intrinsic dsin
  intrinsic dsqrt
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  T_PROJECT_COEFFICIENTS computes the Chebyshev coefficients'
  write ( *, '(a)' ) '  of a function defined over [-1,+1].'
  write ( *, '(a)' ) '  T_PROJECT_COEFFICIENTS_AB works in [A,B].'

  n = 3
  allocate ( c(0:n) )
  call t_project_coefficients ( n, dexp, c )
  call r8vec_print ( n + 1, c, '  Chebyshev coefficients for exp(x) in [-1,+1]' )
  deallocate ( c )

  n = 5
  allocate ( c(0:n) )
  call t_project_coefficients ( n, dexp, c )
  call r8vec_print ( n + 1, c, '  Chebyshev coefficients for exp(x) in [-1,+1]' )
  deallocate ( c )

  n = 5
  allocate ( c(0:n) )
  call t_project_coefficients ( n, dsin, c )
  call r8vec_print ( n + 1, c, '  Chebyshev coefficients for sin(x) in [-1,+1]' )
  deallocate ( c )
!
!  Repeat calculation with T_PROJECT_COEFFICIENTS_AB.
!
  n = 5
  allocate ( c(0:n) )
  a = -1.0D+00
  b = +1.0D+00
  call t_project_coefficients_ab ( n, dsin, a, b, c )
  call r8vec_print ( n + 1, c, '  Chebyshev coefficients for sin(x) in [-1,+1]' )
  deallocate ( c )
!
!  Now try a different interval.
!
  n = 5
  allocate ( c(0:n) )
  a = 0.0D+00
  b = 1.0D+00
  call t_project_coefficients_ab ( n, dsqrt, a, b, c )
  call r8vec_print ( n + 1, c, '  Chebyshev coefficients for sqrt(x) in [0,+1]' )
  deallocate ( c )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests T_PROJECT_COEFFICIENTS_DATA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: c(:)
  real ( kind = 8 ), allocatable :: d(:)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  T_PROJECT_COEFFICIENTS_DATA computes the Chebyshev'
  write ( *, '(a)' ) '  coefficients of a function defined by data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We are looking for an approximation that is good in [-1,+1].'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Begin by using equally spaced points in [-1,+1].'

  a = -1.0D+00
  b = +1.0D+00
  m = 10
  allocate ( x(1:m) )
  allocate ( d(1:m) )
  call r8vec_linspace ( m, a, b, x )
  d(1:m) = exp ( x(1:m) )
  n = 3
  allocate ( c(0:n) )
  call t_project_coefficients_data ( a, b, m, n, x, d, c )
  call r8vec_print ( n + 1, c, '  Chebyshev coefficients for exp(x) on [-1,+1]' )
  deallocate ( c )
  deallocate ( d )
  deallocate ( x )

  a = -1.0D+00
  b = +1.0D+00
  m = 10
  allocate ( x(1:m) )
  allocate ( d(1:m) )
  call r8vec_linspace ( m, a, b, x )
  d(1:m) = exp ( x(1:m) )
  n = 5
  allocate ( c(0:n) )
  call t_project_coefficients_data ( a, b, m, n, x, d, c )
  call r8vec_print ( n + 1, c, '  Chebyshev coefficients for exp(x) on [-1,+1]' )
  deallocate ( c )
  deallocate ( d )
  deallocate ( x )

  a = -1.0D+00
  b = +1.0D+00
  m = 10
  allocate ( x(1:m) )
  allocate ( d(1:m) )
  call r8vec_linspace ( m, a, b, x )
  d(1:m) = sin ( x(1:m) )
  n = 5
  allocate ( c(0:n) )
  call t_project_coefficients_data ( a, b, m, n, x, d, c )
  call r8vec_print ( n + 1, c, '  Chebyshev coefficients for sin(x) on [-1,+1]' )
  deallocate ( c )
  deallocate ( d )
  deallocate ( x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now sample equally spaced points in [0,+1].'
  write ( *, '(a)' ) '  The approximation still applies to the interval [-1,+1].'

  a = 0.0D+00
  b = +1.0D+00
  m = 10
  allocate ( x(1:m) )
  allocate ( d(1:m) )
  call r8vec_linspace ( m, a, b, x )
  d(1:m) = sin ( x(1:m) )
  n = 5
  allocate ( c(0:n) )
  call t_project_coefficients_data ( a, b, m, n, x, d, c )
  call r8vec_print ( n + 1, c, '  Chebyshev coefficients for sin(x) on [0,+1]' )
  deallocate ( c )
  deallocate ( d )
  deallocate ( x )

  a = 0.0D+00
  b = +1.0D+00
  m = 10
  allocate ( x(1:m) )
  allocate ( d(1:m) )
  call r8vec_linspace ( m, a, b, x )
  d(1:m) = sqrt ( x(1:m) )
  n = 5
  allocate ( c(0:n) )
  call t_project_coefficients_data ( a, b, m, n, x, d, c )
  call r8vec_print ( n + 1, c, '  Chebyshev coefficients for sqrt(x) on [0,+1]' )
  deallocate ( c )
  deallocate ( d )
  deallocate ( x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now random points in [-1,+1].'

  a = -1.0D+00
  b = +1.0D+00
  m = 10
  allocate ( x(1:m) )
  allocate ( d(1:m) )
  seed = 123456789
  call r8vec_uniform_01 ( m, seed, x )
  x = x * b + ( 1.0D+00 - x ) * a
  d(1:m) = sin ( x(1:m) )
  n = 5
  allocate ( c(0:n) )
  call t_project_coefficients_data ( a, b, m, n, x, d, c )
  call r8vec_print ( n + 1, c, '  Chebyshev coefficients for sin(x) on [-1,+1]' )
  deallocate ( c )
  deallocate ( d )
  deallocate ( x )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 compares a function and projection over [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: c(:)
  intrinsic dexp
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: r(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09:'
  write ( *, '(a)' ) '  T_PROJECT_COEFFICIENTS computes the Chebyshev interpolant C(F)(N,X)'
  write ( *, '(a)' ) '  of a function F(X) defined over [-1,+1].'
  write ( *, '(a)' ) '  T_PROJECT_VALUE evaluates that projection.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute projections of order N to exp(x) over [-1,+1],'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N   Max||F(X)-C(F)(N,X)||'
  write ( *, '(a)' ) ' '

  a = -1.0D+00
  b = +1.0D+00

  do n = 0, 10
    allocate ( c(0:n) )
    call t_project_coefficients ( n, dexp, c )
    m = 101
    allocate ( x(1:m) )
    call r8vec_linspace ( m, a, b, x )
    allocate ( v(1:m) )
    call t_project_value ( m, n, x, c, v )
    allocate ( r(1:m) )
    r(1:m) = v(1:m) - exp ( x(1:m) )
    write ( *, '(2x,i2,2x,g12.4)' ) n, maxval ( abs ( r ) )
    deallocate ( c )
    deallocate ( r )
    deallocate ( v )
    deallocate ( x )
  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 compares a function and projection over [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable :: c(:)
  intrinsic dexp
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: r(:)
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10:'
  write ( *, '(a)' ) '  T_PROJECT_COEFFICIENTS_AB computes the Chebyshev interpolant C(F)(N,X)'
  write ( *, '(a)' ) '  of a function F(X) defined over [A,B].'
  write ( *, '(a)' ) '  T_PROJECT_VALUE_AB evaluates that projection.'

  a = 0.0D+00
  b = 1.5D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,f4.1,a,f4.1,a)' ) &
    '  Compute projections of order N to exp(x) over [', a, ',', b, ']'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N   Max||F(X)-C(F)(N,X)||'
  write ( *, '(a)' ) ' '

  do n = 0, 10
    allocate ( c(0:n) )
    call t_project_coefficients_ab ( n, dexp, a, b, c )
    m = 101
    allocate ( x(1:m) )
    call r8vec_linspace ( m, a, b, x )
    allocate ( v(1:m) )
    call t_project_value_ab ( m, n, x, c, a, b, v )
    allocate ( r(1:m) )
    r(1:m) = v(1:m) - exp ( x(1:m) )
    write ( *, '(2x,i2,2x,g12.4)' ) n, maxval ( abs ( r ) )
    deallocate ( c )
    deallocate ( r )
    deallocate ( v )
    deallocate ( x )
  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests U_POLYNOMIAL_COEFFICIENTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  U_POLYNOMIAL_COEFFICIENTS determines coefficients'
  write ( *, '(a)' ) '  for Chebyshev polynomials U(n,x).'

  call u_polynomial_coefficients ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  T(', i, ',x) = '
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( c(i,j) == 0.0D+00 ) then

      else if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x ^ ', j
      end if
    end do
  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests U_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ), allocatable :: fx2(:)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12:'
  write ( *, '(a)' ) '  U_POLYNOMIAL evaluates the Chebyshev polynomials U(n,x).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                    Tabulated     Computed'
  write ( *, '(a)' ) '     N      X        U(n,x)        U(n,x)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call u_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( fx2(0:n) )

    call u_polynomial ( 1, n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

    deallocate ( fx2 )

  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests U_POLYNOMIAL_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: fx(:,:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: n_max = 5
  real ( kind = 8 ), allocatable :: z(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13:'
  write ( *, '(a)' ) '  U_POLYNOMIAL_ZEROS returns zeroes of U(n,x).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      X        U(n,x)'

  do n = 1, n_max

    allocate ( z(1:n) )

    call u_polynomial_zeros ( n, z )

    allocate ( fx(1:n,0:n) )

    call u_polynomial ( n, n, z, fx )

    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i8,2x,f8.4,2x,g14.6)' ) n, z(i), fx(i,n)
    end do

    deallocate ( fx )
    deallocate ( z )

  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests U_QUADRATURE_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e
  real ( kind = 8 ), allocatable :: f(:)
  integer ( kind = 4 ) n
  real ( kind = 8 ) q
  real ( kind = 8 ) q_exact
  real ( kind = 8 ) u_integral
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14:'
  write ( *, '(a)' ) '  U_QUADRATURE_RULE computes the quadrature rule'
  write ( *, '(a)' ) '  associated with U(n,x);'

  n = 7
  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call u_quadrature_rule ( n, x, w )

  call r8vec2_print ( n, x, w, '      X            W' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use the quadrature rule to estimate:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Q = Integral ( -1 <= X <= +1 ) X^E * sqrt ( 1-x^2) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   E       Q_Estimate      Q_Exact'
  write ( *, '(a)' ) ' '

  allocate ( f(1:n) )

  do e = 0, 2 * n - 1
    if ( e == 0 ) then
      f(1:n) = 1.0D+00
    else
      f(1:n) = x(1:n) ** e
    end if
    q = dot_product ( w, f )
    q_exact = u_integral ( e )
    write ( *, '(2x,i2,g14.6,2x,g14.6)' ) e, q, q_exact
  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests V_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ), allocatable :: fx2(:)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15:'
  write ( *, '(a)' ) '  V_POLYNOMIAL evaluates Chebyshev polynomials V(n,x).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                    Tabulated     Computed'
  write ( *, '(a)' ) '     N      X        V(n,x)        V(n,x)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call v_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( fx2(0:n) )

    call v_polynomial ( 1, n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

    deallocate ( fx2 )

  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests W_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ), allocatable :: fx2(:)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16:'
  write ( *, '(a)' ) '  W_POLYNOMIAL evaluates Chebyshev polynomials W(n,x).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                    Tabulated     Computed'
  write ( *, '(a)' ) '     N      X        W(n,x)        W(n,x)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call w_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( fx2(0:n) )

    call w_polynomial ( 1, n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

    deallocate ( fx2 )

  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests T_TRIPLE_PRODUCT_INTEGRAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t_polynomial_value
  real ( kind = 8 ) t_triple_product_integral
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 20
  real ( kind = 8 ) ti
  real ( kind = 8 ) tj
  real ( kind = 8 ) tk
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17:'
  write ( *, '(a)' ) '  T_TRIPLE_PRODUCT_INTEGRAL computes the triple integral'
  write ( *, '(a)' ) '  Tijk = integral ( -1 <= x <= 1 ) T(i,x) T(j,x) T(k,x) / sqrt ( 1-x^2) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J   K     Tijk           Tijk'
  write ( *, '(a)' ) '                 computed       exact'
  write ( *, '(a)' ) ' '

  n = 15
  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call t_quadrature_rule ( n, x, w )

  seed = 123456789

  do test = 1, test_num
    i = i4_uniform ( 2, 6, seed )
    j = i4_uniform ( 1, 3, seed )
    k = i4_uniform ( 0, 4, seed )
    fx1 = t_triple_product_integral ( i, j, k )
    fx2 = 0.0D+00
    do l = 1, n
      ti = t_polynomial_value ( i, x(l) )
      tj = t_polynomial_value ( j, x(l) )
      tk = t_polynomial_value ( k, x(l) )
      fx2 = fx2 + w(l) * ti * tj * tk
    end do
    write ( *, '(2x,i2,2x,i2,2x,i2,2x,g14.6,2x,g14.6)' ) i, j, k, fx1, fx2
  end do

  deallocate ( x )
  deallocate ( w )

  return
end
