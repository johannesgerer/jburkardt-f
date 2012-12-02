program main

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_PRB tests the LAGUERRE_POLYNOMIAL library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2012
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
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the LAGUERRE_POLYNOMIAL library.'

  call laguerre_polynomial_test01 ( )
  call laguerre_polynomial_test02 ( )
  call laguerre_polynomial_test03 ( )
  call laguerre_polynomial_test04 ( )
  call laguerre_polynomial_test05 ( )
  call laguerre_polynomial_test06 ( )

  p = 5
  b = 0.0D+00
  call laguerre_polynomial_test07 ( p, b )

  p = 5
  b = 1.0D+00
  call laguerre_polynomial_test07 ( p, b )

  p = 5
  e = 0
  call laguerre_polynomial_test08 ( p, e )

  p = 5
  e = 1
  call laguerre_polynomial_test08 ( p, e )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine laguerre_polynomial_test01 ( )

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_TEST01 tests L_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_TEST01:'
  write ( *, '(a)' ) '  L_POLYNOMIAL_VALUES stores values of'
  write ( *, '(a)' ) '  the Laguerre polynomials.'
  write ( *, '(a)' ) '  L_POLYNOMIAL evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                        Tabulated                 Computed'
  write ( *, '(a)' ) '     N        X           L(N,X)                    L(N,X)                     Error'

  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call l_polynomial_values ( n_data, n, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( v(m,n+1) )
    call l_polynomial ( m, n, x, v )
    fx2 = v(1,n+1)
    deallocate ( v )

    e = fx1 - fx2

    write ( *, '(2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, x, fx1, fx2, e

  end do

  return
end
subroutine laguerre_polynomial_test02 ( )

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_TEST02 tests L_POLYNOMIAL_COEFFICIENTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2012
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
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_TEST02'
  write ( *, '(a)' ) &
    '  L_POLYNOMIAL_COEFFICIENTS determines polynomial coefficients of L(n,x).'

  call l_polynomial_coefficients ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  L(', i, ') = '
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
subroutine laguerre_polynomial_test03 ( )

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_TEST03 tests L_POLYNOMIAL_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_TEST03:'
  write ( *, '(a)' ) '  L_POLYNOMIAL_ZEROS computes the zeros of L(n,x)'
  write ( *, '(a)' ) '  Check by calling L_POLYNOMIAL there.'

  do degree = 1, 5

    allocate ( z(1:degree) )
    call l_polynomial_zeros ( degree, z )
    write ( title, '(a,i1,a)' ) '  Computed zeros for L(', degree, ',z):'
    call r8vec_print ( degree, z, title )

    allocate ( lz(degree,0:degree) )
    call l_polynomial ( degree, degree, z, lz )
    write ( title, '(a,i1,a)' ) '  Evaluate L(', degree, ',z):'
    call r8vec_print ( degree, lz(1:degree,degree), title )

    deallocate ( lz )
    deallocate ( z )

  end do

  return
end
subroutine laguerre_polynomial_test04 ( )

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_TEST04 tests L_QUADRATURE_RULE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_TEST04:'
  write ( *, '(a)' ) '  L_QUADRATURE_RULE computes the quadrature rule'
  write ( *, '(a)' ) '  associated with L(n,x)'

  n = 7
  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call l_quadrature_rule ( n, x, w )

  call r8vec2_print ( n, x, w,  '      X            W' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use the quadrature rule to estimate:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Q = Integral ( 0 <= X < +00 ) X^E exp(-X) dx'
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
    call l_integral ( e, q_exact )
    write ( *, '(2x,i2,2x,g14.6,2x,g14.6)' ) e, q, q_exact
  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine laguerre_polynomial_test05 ( )

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_TEST05 tests LM_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2012
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
  integer ( kind = 4 ), parameter :: mm = 1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: v(:,:)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_TEST05:'
  write ( *, '(a)' ) '  LM_POLYNOMIAL_VALUES stores values of'
  write ( *, '(a)' ) '  the Laguerre polynomial Lm(n,m,x)'
  write ( *, '(a)' ) '  LM_POLYNOMIAL evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                                 Tabulated                 Computed'
  write ( *, '(a)' ) '     N     M        X            Lm(N,M,X)                 Lm(N,M,X)               Error'

  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lm_polynomial_values ( n_data, n, m, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( v(mm,n+1) )
    call lm_polynomial ( mm, n, m, x, v )
    fx2 = v(1,n+1)
    deallocate ( v )

    e = fx1 - fx2

    write ( *, '(2x,i4,2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, m, x, fx1, fx2, e

  end do

  return
end
subroutine laguerre_polynomial_test06 ( )

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_TEST06 tests LM_POLYNOMIAL_COEFFICIENTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2012
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
  integer ( kind = 4 ) m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_TEST06'
  write ( *, '(a)' ) &
    '  LM_POLYNOMIAL_COEFFICIENTS determines polynomial coefficients of Lm(n,m,x).'

  do m = 0, 4

    call lm_polynomial_coefficients ( n, m, c )
 
    do i = 0, n
      write ( *, '(a)' ) ' '
      write ( *, '(a,i2,a,i2,a)' ) '  Lm(', i, ',', m, ') = '
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

  end do
 
  return
end
subroutine laguerre_polynomial_test07 ( p, b )

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_TEST07 tests L_EXPONENTIAL_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2012
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
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_TEST07'
  write ( *, '(a)' ) '  Compute an exponential product table for L(n,x):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( 0 <= x < +oo ) exp(b*x) Ln(i,x) Ln(j,x) exp(-x) dx'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponential argument coefficient B = ', b

  allocate ( table(0:p,0:p) )

  call l_exponential_product ( p, b, table )

  call r8mat_print ( p + 1, p + 1, table, '  Exponential product table:' )

  deallocate ( table )

  return
end
subroutine laguerre_polynomial_test08 ( p, e )

!*****************************************************************************80
!
!! LAGUERRE_POLYNOMIAL_TEST08 tests L_POWER_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2012
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
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_TEST08'
  write ( *, '(a)' ) '  Compute a power product table for L(n,x):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( 0 <= x < +oo ) x^e L(i,x) L(j,x) exp(-x) dx'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponent of X, E = ', e

  allocate ( table(0:p,0:p) )

  call l_power_product ( p, e, table )

  call r8mat_print ( p + 1, p + 1, table, '  Power product table:' )

  deallocate ( table )

  return
end
