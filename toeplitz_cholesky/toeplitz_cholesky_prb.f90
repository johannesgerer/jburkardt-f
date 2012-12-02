program main

!*****************************************************************************80
!
!! MAIN tests TOEPLITZ_CHOLESKY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOEPLITZ_CHOLESKY library.'

  call toeplitz_cholesky_prb01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine toeplitz_cholesky_prb01 ( )

!*****************************************************************************80
!
!! TOEPLITZ_CHOLESKY_PRB01 tests TOEPLITZ_CHOLESKY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) g(2,n)
  real ( kind = 8 ) l(n,n)
  real ( kind = 8 ) r(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_PRB01'
  write ( *, '(a)' ) '  Test the factorization of a simple example.'
!
!  TOEPLITZ_CHOLESKY_UPPER.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_UPPER:'

  a = reshape ( (/ &
        1.0,   0.5, -0.375, &
        0.5,   1.0,  0.5, &
       -0.375, 0.5,  1.0 /), (/ 3, 3 /) )

  call r8mat_print ( n, n, a, '  Toeplitz matrix A:' )

  call toeplitz_cholesky_upper ( n, a, r )
  call r8mat_print ( n, n, r, '  Computed upper Cholesky factor R:' )

  b = matmul ( transpose ( r ), r )
  call r8mat_print ( n, n, b, '  Product R''R:' )
!
!  TOEP_CHOLESKY_UPPER.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEP_CHOLESKY_UPPER:'
  
  g = reshape ( (/ &
    1.0, 0.0,  &
    0.5, 0.5,  &
   -0.375, -0.375 /), (/ 2, 3 /) )

  call r8mat_print ( 2, n, g, '  Compressed Toeplitz matrix G:' )

  call toep_cholesky_upper ( n, g, r )
  call r8mat_print ( n, n, r, '  Computed upper Cholesky factor R:' )

  b = matmul ( transpose ( r ), r )
  call r8mat_print ( n, n, b, '  Product R''R:' )
!
!  TOEPLITZ_CHOLESKY_LOWER.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_LOWER:'

  a = reshape ( (/ &
        1.0,   0.5, -0.375, &
        0.5,   1.0,  0.5, &
       -0.375, 0.5,  1.0 /), (/ 3, 3 /) )

  call r8mat_print ( n, n, a, '  Toeplitz matrix A:' )

  call toeplitz_cholesky_lower ( n, a, l )
  call r8mat_print ( n, n, l, '  Computed lower Cholesky factor L:' )

  b = matmul ( l, transpose ( l ) )
  call r8mat_print ( n, n, b, '  Product LL'':' )
!
!  TOEP_CHOLESKY_LOWER.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEP_CHOLESKY_LOWER:'

  g = reshape ( (/ &
    1.0, 0.0,  &
    0.5, 0.5,  &
   -0.375, -0.375 /), (/ 2, 3 /) )

  call r8mat_print ( 2, n, g, '  Compressed Toeplitz matrix G:' )

  call toep_cholesky_lower ( n, g, l )
  call r8mat_print ( n, n, l, '  Computed lower Cholesky factor L:' )

  b = matmul ( l, transpose ( l ) )
  call r8mat_print ( n, n, b, '  Product LL'':' )

  return
end



