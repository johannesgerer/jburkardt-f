subroutine expm1 ( n, a, e )

!*****************************************************************************80
!
!! EXPM1 is essentially MATLAB's built-in matrix exponential algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    Cleve Moler, Charles Van Loan
!
!  Reference:
!
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) E(N,N), the estimate for exp(A).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a2(n,n)
  real ( kind = 8 ) a_norm
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n,n)
  real ( kind = 8 ) e(n,n)
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) k
  logical p
  integer ( kind = 4 ) , parameter :: q = 6
  real ( kind = 8 ) r8_log_2
  real ( kind = 8 ) r8mat_norm_li
  integer ( kind = 4 ) s
  real ( kind = 8 ) x(n,n)

  a2(1:n,1:n) = a(1:n,1:n)

  a_norm = r8mat_norm_li ( n, n, a2 )

  ee = int ( r8_log_2 ( a_norm ) ) + 1

  s = max ( 0, ee + 1 )

  a2(1:n,1:n) = a2(1:n,1:n) / 2.0D+00**s

  x(1:n,1:n) = a2(1:n,1:n)

  c = 0.5D+00

  call r8mat_identity ( n, e )
  e(1:n,1:n) = e(1:n,1:n) + c * a2(1:n,1:n)

  call r8mat_identity ( n, d )
  d(1:n,1:n) = d(1:n,1:n) - c * a2(1:n,1:n)

  p = .true.

  do k = 2, q

    c = c * real ( q - k + 1, kind = 8 ) &
      / real ( k * ( 2 * q - k + 1 ), kind = 8 )

    x(1:n,1:n) = matmul ( a2(1:n,1:n), x(1:n,1:n) )

    e(1:n,1:n) = e(1:n,1:n) + c * x(1:n,1:n)

    if ( p ) then
      d(1:n,1:n) = d(1:n,1:n) + c * x(1:n,1:n)
    else
      d(1:n,1:n) = d(1:n,1:n) - c * x(1:n,1:n)
    end if

    p = .not. p

  end do
!
!  E -> inverse(D) * E
!
  call r8mat_minvm ( n, n, d, e, e )
!
!  E -> E^(2*S)
!
  do k = 1, s
    e(1:n,1:n) = matmul ( e(1:n,1:n), e(1:n,1:n) )
  end do

  return
end
subroutine expm2 ( n, a, e )

!*****************************************************************************80
!
!! EXPM2 uses the Taylor series for the matrix exponential.
!
!  Discussion:
!
!    Formally,
!
!      exp ( A ) = I + A + 1/2 A^2 + 1/3! A^3 + ...
!
!    This function sums the series until a tolerance is satisfied.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    Cleve Moler, Charles Van Loan
!
!  Reference:
!
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) E(N,N), the estimate for exp(A).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) e(n,n)
  real ( kind = 8 ) f(n,n)
  real ( kind = 8 ) g(n,n)
  integer ( kind = 4 ) k
  logical r8mat_insignificant

  e(1:n,1:n) = 0.0D+00

  call r8mat_identity ( n, f )

  k = 1

  do

    if ( r8mat_insignificant ( n, n, e, f ) ) then
      exit
    end if

    e(1:n,1:n) = e(1:n,1:n) + f(1:n,1:n)

    f(1:n,1:n) = matmul ( a(1:n,1:n), f(1:n,1:n) ) / real ( k, kind = 8 )
    k = k + 1

  end do

  return
end
subroutine expm3 ( n, a, e )

!*****************************************************************************80
!
!! EXPM3 approximates the matrix exponential using an eigenvalue approach.
!
!  Discussion:
!
!    exp(A) = V * D * V
!
!    where V is the matrix of eigenvectors of A, and D is the diagonal matrix
!    whose i-th diagonal entry is exp(lambda(i)), for lambda(i) an eigenvalue
!    of A.
!
!    This function is accurate for matrices which are symmetric, orthogonal,
!    or normal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    Cleve Moler, Charles Van Loan
!
!  Reference:
!
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) E(N,N), the estimate for exp(A).
!
! [ V, D ] = eig ( A );
! E = V * diag ( exp ( diag ( D ) ) ) / V;
  return
end
