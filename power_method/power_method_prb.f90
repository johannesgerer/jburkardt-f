program main

!*****************************************************************************80
!
!! MAIN is the main program for POWER_METHOD_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POWER_METHOD_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the POWER_METHOD library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POWER_METHOD_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses POWER_METHOD on the Fibonacci2 matrix.
!
!  Discussion:
!
!    This matrix, despite having a single dominant eigenvalue, will generally
!    converge only very slowly under the power method.  This has to do with
!    the fact that the matrix has only 3 eigenvectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) cos_x1x2
  real    ( kind = 8 ) ctime
  real    ( kind = 8 ) ctime1
  real    ( kind = 8 ) ctime2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real    ( kind = 8 ) lambda
  real    ( kind = 8 ) phi
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) sin_x1x2
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) x2(n)

  call fibonacci2 ( n, a )

  seed = 123456789

  call r8vec_uniform_01 ( n, seed, x )

  it_max = 300
  tol = 0.000001D+0

  phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use POWER_METHOD on the Fibonacci2 matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Matrix order N         = ', n
  write ( *, '(a,i8)'      ) '  Maximum iterations     = ', it_max
  write ( *, '(a,g14.6)'   ) '  Error tolerance        = ', tol

  call cpu_time ( ctime1 )

  call power_method ( n, a, x, it_max, tol, lambda, it_num )

  call cpu_time ( ctime2 )
  ctime = ctime2 - ctime1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Number of iterations   = ', it_num
  write ( *, '(a,g14.6)'   ) '  CPU time               = ', ctime
  write ( *, '(a,f14.10)'  ) '  Estimated eigenvalue   = ', lambda
  write ( *, '(a,f14.10)'  ) '  Correct value          = ', phi
  write ( *, '(a,g14.6)'   ) '  Error                  = ', abs ( lambda - phi )
!
!  X2 is the exact eigenvector.
!
  x2(1) = 1.0
  do i = 2, n
    x2(i) = phi * x2(i-1)
  end do

  x2(1:n) = x2(1:n) / sqrt ( sum ( x2(1:n)**2 ) )
!
!  The sine of the angle between X and X2 is a measure of error.
!
  cos_x1x2 = dot_product ( x(1:n), x2(1:n) )
  sin_x1x2 = sqrt ( ( 1.0 - cos_x1x2 ) * ( 1.0 + cos_x1x2 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.10)' ) &
    '  Sine of angle between true and estimated vectors = ', sin_x1x2

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 uses POWER_METHOD2 on the Fibonacci2 matrix.
!
!  Discussion:
!
!    This matrix, despite having a single dominant eigenvalue, will generally
!    converge only very slowly under the power method.  This has to do with
!    the fact that the matrix has only 3 eigenvectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) ctime
  real    ( kind = 8 ) ctime1
  real    ( kind = 8 ) ctime2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  complex ( kind = 8 ) lambda
  real    ( kind = 8 ) phi
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) tol
  complex ( kind = 8 ) v(n)
  real    ( kind = 8 ) x(n)

  call fibonacci2 ( n, a )

  seed = 123456789
  call r8vec_uniform_01 ( n, seed, x )

  it_max = 300
  tol = 0.000001D+0

  phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Use POWER_METHOD2 on the Fibonacci2 matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Matrix order N         = ', n
  write ( *, '(a,i8)'      ) '  Maximum iterations     = ', it_max
  write ( *, '(a,g14.6)'   ) '  Error tolerance        = ', tol

  call cpu_time ( ctime1 )

  call power_method2 ( n, a, x, it_max, tol, lambda, v, it_num )

  call cpu_time ( ctime2 )
  ctime = ctime2 - ctime1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Number of iterations   = ', it_num
  write ( *, '(a,g14.6)'   ) '  CPU time               = ', ctime
  write ( *, '(a,2f14.10)' ) '  Estimated eigenvalue   = ', lambda
  write ( *, '(a,f14.10)'  ) '  Correct value          = ', phi
  write ( *, '(a,g14.6)'   ) '  Error                  = ', abs ( lambda - phi )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 uses POWER_METHOD2 on the TRIS matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  real    ( kind = 8 ) cos_y1y2
  real    ( kind = 8 ) ctime
  real    ( kind = 8 ) ctime1
  real    ( kind = 8 ) ctime2
  real    ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  complex ( kind = 8 ) lambda
  complex ( kind = 8 ) lambda_max
  complex ( kind = 8 ) lambda_vec(n)
  real    ( kind = 8 ) phi
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) sin_y1y2
  real    ( kind = 8 ) tol
  complex ( kind = 8 ) v(n)
  real    ( kind = 8 ) x(n)
!
!  If ALPHA * GAMMA is negative, this matrix will have complex eigenvalues.
!
  alpha = -1.0D+00
  beta = 10.0D+00
  gamma = +8.0D+0

  call tris ( n, n, alpha, beta, gamma, a )

  seed = 123456789

  call r8vec_uniform_01 ( n, seed, x )

  it_max = 4000
  tol = 0.000001D+0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Use POWER_METHOD2 on the TRIS (tridiagonal scalar) matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Matrix order N         = ', n
  write ( *, '(a,i8)'      ) '  Maximum iterations     = ', it_max
  write ( *, '(a,g14.6)'   ) '  Error tolerance        = ', tol

  call cpu_time ( ctime1 )

  call power_method2 ( n, a, x, it_max, tol, lambda, v, it_num )

  call cpu_time ( ctime2 )
  ctime = ctime2 - ctime1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Number of iterations   = ', it_num
  write ( *, '(a,g14.6)'   ) '  CPU time               = ', ctime
  write ( *, '(a,2f14.10)' ) '  Estimated eigenvalue   = ', lambda

  call tris_eigenvalues ( n, alpha, beta, gamma, lambda_vec )

  lambda_max = lambda_vec(1)
  do i = 2, n
    if ( abs ( lambda_max ) < abs ( lambda_vec(i) ) ) then
      lambda_max = lambda_vec(i)
    end if
  end do

  write ( *, '(a,2f14.10)' ) '  Correct max eigenvalue = ', lambda_max
  write ( *, '(a,g14.6)'   ) '  Error                  = ', abs ( lambda - lambda_max )

  return
end
