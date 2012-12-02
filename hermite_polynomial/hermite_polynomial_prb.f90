program main

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_PRB tests the HERMITE_POLYNOMIAL library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
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
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the HERMITE_POLYNOMIAL library.'

  call hermite_polynomial_test01 ( )
  call hermite_polynomial_test02 ( )
  call hermite_polynomial_test03 ( )
  call hermite_polynomial_test04 ( )
  call hermite_polynomial_test05 ( )
  call hermite_polynomial_test06 ( )
  call hermite_polynomial_test07 ( )

  p = 5
  b = 0.0D+00
  call hermite_polynomial_test08 ( p, b )

  p = 5
  b = 1.0D+00
  call hermite_polynomial_test08 ( p, b )

  p = 5
  e = 0
  call hermite_polynomial_test09 ( p, e )

  p = 5
  e = 1
  call hermite_polynomial_test09 ( p, e )

  p = 5
  b = 0.0D+00
  call hermite_polynomial_test10 ( p, b )

  p = 5
  b = 1.0D+00
  call hermite_polynomial_test10 ( p, b )

  p = 5
  e = 0
  call hermite_polynomial_test11 ( p, e )

  p = 5
  e = 1
  call hermite_polynomial_test11 ( p, e )

  p = 5
  b = 0.0D+00
  call hermite_polynomial_test12 ( p, b )

  p = 5
  b = 1.0D+00
  call hermite_polynomial_test12 ( p, b )

  p = 5
  e = 0
  call hermite_polynomial_test13 ( p, e )

  p = 5
  e = 1
  call hermite_polynomial_test13 ( p, e )

  call hermite_polynomial_test14 ( )

  call hermite_polynomial_test15 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine hermite_polynomial_test01 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST01 tests H_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2012
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
  real ( kind = 8 ), allocatable :: fx2_vec(:)
  integer ( kind = 4 ) n
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST01:'
  write ( *, '(a)' ) '  H_POLYNOMIAL_VALUES stores values of'
  write ( *, '(a)' ) '  the physicist''s Hermite polynomials.'
  write ( *, '(a)' ) '  H_POLYNOMIAL evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                        Tabulated                 Computed'
  write ( *, '(a)' ) '     N        X           H(N,X)                    H(N,X)                     Error'

  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call h_polynomial_values ( n_data, n, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( fx2_vec(n+1) )
    call h_polynomial ( 1, n, x, fx2_vec )
    fx2 = fx2_vec(n+1)
    deallocate ( fx2_vec )

    e = fx1 - fx2

    write ( *, '(2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, x, fx1, fx2, e

  end do

  return
end
subroutine hermite_polynomial_test02 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST02 tests HE_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2012
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
  real ( kind = 8 ), allocatable :: fx2_vec(:)
  integer ( kind = 4 ) n
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST02:'
  write ( *, '(a)' ) '  HE_POLYNOMIAL_VALUES stores values of'
  write ( *, '(a)' ) '  the probabilist''s Hermite polynomials.'
  write ( *, '(a)' ) '  HE_POLYNOMIAL evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                        Tabulated                 Computed'
  write ( *, '(a)' ) '     N        X          He(N,X)' // &
    '                   He(N,X)                     Error'

  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call he_polynomial_values ( n_data, n, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( fx2_vec(1:n+1) )
    call he_polynomial ( 1, n, x, fx2_vec )
    fx2 = fx2_vec(n+1)
    deallocate ( fx2_vec )
    e = fx1 - fx2

    write ( *, '(2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, x, fx1, fx2, e

  end do

  return
end
subroutine hermite_polynomial_test03 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST03 tests HF_FUNCTION.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
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
  real ( kind = 8 ), allocatable :: fx2_vec(:)
  integer ( kind = 4 ) n
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST03:'
  write ( *, '(a)' ) '  HF_FUNCTION_VALUES stores values of'
  write ( *, '(a)' ) '  the Hermite function Hf(n,x).'
  write ( *, '(a)' ) '  HF_FUNCTION evaluates the function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                        Tabulated                 Computed'
  write ( *, '(a)' ) '     N        X          Hf(N,X)' // &
    '                   Hf(N,X)                     Error'

  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hf_function_values ( n_data, n, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( fx2_vec(n+1) )
    call hf_function ( 1, n, x, fx2_vec )
    fx2 = fx2_vec(n+1)
    deallocate ( fx2_vec )

    e = fx1 - fx2

    write ( *, '(2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, x, fx1, fx2, e

  end do

  return
end
subroutine hermite_polynomial_test04 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST04 tests H_POLYNOMIAL_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) degree
  real ( kind = 8 ), allocatable :: hz(:,:)
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable :: z(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST04:'
  write ( *, '(a)' ) '  H_POLYNOMIAL_ZEROS computes the zeros of H(n,x)'
  write ( *, '(a)' ) '  Check by calling H_POLYNOMIAL there.'

  do degree = 1, 5

    allocate ( z(1:degree) )
    call h_polynomial_zeros ( degree, z )
    write ( title, '(a,i1,a)' ) '  Computed zeros for H(', degree, ',z):'
    call r8vec_print ( degree, z, title )

    allocate ( hz(degree,0:degree) )
    call h_polynomial ( degree, degree, z, hz )
    write ( title, '(a,i1,a)' ) '  Evaluate H(', degree, ',z):'
    call r8vec_print ( degree, hz(1:degree,degree), title )

    deallocate ( hz )
    deallocate ( z )

  end do

  return
end
subroutine hermite_polynomial_test05 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST05 tests HE_POLYNOMIAL_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) degree
  real ( kind = 8 ), allocatable :: hz(:,:)
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable :: z(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST05:'
  write ( *, '(a)' ) '  HE_POLYNOMIAL_ZEROS computes the zeros of He(n,x)'
  write ( *, '(a)' ) '  Check by calling HE_POLYNOMIAL there.'

  do degree = 1, 5

    allocate ( z(1:degree) )
    call he_polynomial_zeros ( degree, z )
    write ( title, '(a,i1,a)' ) '  Computed zeros for He(', degree, ',z):'
    call r8vec_print ( degree, z, title )

    allocate ( hz(degree,0:degree) )
    call he_polynomial ( degree, degree, z, hz )
    write ( title, '(a,i1,a)' ) '  Evaluate He(', degree, ',z):'
    call r8vec_print ( degree, hz(1:degree,degree), title )

    deallocate ( hz )
    deallocate ( z )

  end do

  return
end
subroutine hermite_polynomial_test06 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST06 tests H_QUADRATURE_RULE.
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

  integer ( kind = 4 ) e
  real ( kind = 8 ), allocatable :: f(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) q
  real ( kind = 8 ) q_exact
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST06:'
  write ( *, '(a)' ) '  H_QUADRATURE_RULE computes the quadrature rule'
  write ( *, '(a)' ) '  associated with H(n,x)'

  n = 7
  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call h_quadrature_rule ( n, x, w )

  call r8vec2_print ( n, x, w, '      X            W' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use the quadrature rule to estimate:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Q = Integral ( -oo < X < +00 ) X^E exp(-X^2) dx'
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
    call h_integral ( e, q_exact )
    write ( *, '(2x,i2,2x,g14.6,2x,g14.6)' ) e, q, q_exact
  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine hermite_polynomial_test07 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST07 tests HE_QUADRATURE_RULE.
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
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST07:'
  write ( *, '(a)' ) '  HE_QUADRATURE_RULE computes the quadrature rule'
  write ( *, '(a)' ) '  associated with He(n,x)'

  n = 7
  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call he_quadrature_rule ( n, x, w )

  call r8vec2_print ( n, x, w, '      X            W' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use the quadrature rule to estimate:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Q = Integral ( -oo < X < +00 ) X^E exp(-0.5*X^2) dx'
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
    call he_integral ( e, q_exact )
    write ( *, '(2x,i2,2x,g14.6,2x,g14.6)' ) e, q, q_exact
  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine hermite_polynomial_test08 ( p, b )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST08 tests HN_EXPONENTIAL_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST08'
  write ( *, '(a)' ) '  Compute a normalized physicist''s Hermite exponential product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) exp(B*X) Hn(I,X) Hn(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hn(I,X) = normalized physicist''s Hermite polynomial of degree I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponential argument coefficient B = ', b

  allocate ( table(0:p,0:p) )

  call hn_exponential_product ( p, b, table )

  call r8mat_print ( p + 1, p + 1, table, '  Exponential product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test09 ( p, e )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST09 tests HN_POWER_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST09'
  write ( *, '(a)' ) '  Compute a normalized physicist''s Hermite power product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) X^E Hn(I,X) Hn(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hn(I,X) = normalized physicist''s Hermite polynomial of degree I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponent of X, E = ', e

  allocate ( table(0:p,0:p) )

  call hn_power_product ( p, e, table )

  call r8mat_print ( p + 1, p + 1, table, '  Power product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test10 ( p, b )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST10 tests HEN_EXPONENTIAL_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST10'
  write ( *, '(a)' ) '  Compute a normalized probabilist''s Hermite exponential product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) exp(B*X) Hen(I,X) Hen(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hen(I,X) = normalized probabilist''s Hermite polynomial of degree I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponential argument coefficient B = ', b

  allocate ( table(0:p,0:p) )

  call hen_exponential_product ( p, b, table )

  call r8mat_print ( p + 1, p + 1, table, '  Exponential product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test11 ( p, e )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST11 tests HEN_POWER_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST11'
  write ( *, '(a)' ) '  Compute a normalized probabilist''s Hermite power product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) X^E Hen(I,X) Hen(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hn(I,X) = normalized probabilist''s Hermite polynomial of degree I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponent of X, E = ', e

  allocate ( table(0:p,0:p) )

  call hen_power_product ( p, e, table )

  call r8mat_print ( p + 1, p + 1, table, '  Power product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test12 ( p, b )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST12 tests HF_EXPONENTIAL_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST12'
  write ( *, '(a)' ) '  Compute a Hermite function exponential product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) exp(B*X) Hf(I,X) Hf(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hf(I,X) = Hermite function of "degree" I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponential argument coefficient B = ', b

  allocate ( table(0:p,0:p) )

  call hf_exponential_product ( p, b, table )

  call r8mat_print ( p + 1, p + 1, table, '  Exponential product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test13 ( p, e )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST13 tests HF_POWER_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
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
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST13'
  write ( *, '(a)' ) '  Compute a Hermite function power product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) X^E Hf(I,X) Hf(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hf(I,X) = Hermite function of "degree" I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponent of X, E = ', e

  allocate ( table(0:p,0:p) )

  call hf_power_product ( p, e, table )

  call r8mat_print ( p + 1, p + 1, table, '  Power product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test14 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST14 tests H_POLYNOMIAL_COEFFICIENTS.
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
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST14'
  write ( *, '(a)' ) &
    '  H_POLYNOMIAL_COEFFICIENTS determines the physicist''s Hermite polynomial coefficients.'

  call h_polynomial_coefficients ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  H(', i, ',x) ='
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
subroutine hermite_polynomial_test15 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST15 tests HE_POLYNOMIAL_COEFFICIENTS.
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
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST15'
  write ( *, '(a)' ) &
    '  HE_POLYNOMIAL_COEFFICIENTS determines the probabilist''s Hermite polynomial coefficients.'

  call he_polynomial_coefficients ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  He(', i, ',x) ='
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
