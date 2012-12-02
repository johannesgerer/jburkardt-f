subroutine fibonacci2 ( n, a )

!*****************************************************************************80
!
!! FIBONACCI2 returns the FIBONACCI2 matrix.
!
!  Example:
!
!    N = 5
!
!    0 1 0 0 0
!    1 1 0 0 0
!    0 1 1 0 0
!    0 0 1 1 0
!    0 0 0 1 1
!
!  Properties:
!
!    A is generally not symmetric: A' /= A.
!
!    A is tridiagonal.
!
!    Because A is tridiagonal, it has property A (bipartite).
!
!    A is banded, with bandwidth 3.
!
!    A is integral, therefore det ( A ) is integral, and
!    det ( A ) * inverse ( A ) is integral.
!
!    A is a zero/one matrix.
!
!    If N = 1 then
!      det ( A ) = 0
!    else
!      det ( A ) = -1
!
!    If 1 < N, then A is unimodular.
!
!    When applied to a Fibonacci1 matrix B, the Fibonacci2 matrix
!    A produces the "next" Fibonacci1 matrix C = A*B.
!
!    Let PHI be the golden ratio (1+sqrt(5))/2.
!
!    For 2 <= N, the eigenvalues and eigenvectors are:
!
!    LAMBDA(1)     = PHI,     vector = (1,PHI,PHI^2,...PHI^(N-1));
!    LAMBDA(2:N-1) = 1        vector = (0,0,0,...,0,1);
!    LAMBDA(N)     = 1 - PHI. vector = ((-PHI)^(N-1),(-PHI)^(N-2),...,1)
!
!    Note that there is only one eigenvector corresponding to 1.
!    Hence, for 3 < N, the matrix is defective.  This fact means,
!    for instance, that the convergence of the eigenvector in the power
!    method will be very slow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n

      if ( i == 1 ) then

        if ( j == 2 ) then
          a(i,j) = 1.0D+00
        else
          a(i,j) = 0.0D+00
        end if

      else

        if ( j == i - 1 .or. j == i ) then
          a(i,j) = 1.0D+00
        else
          a(i,j) = 0.0D+00
        end if

      end if

    end do
  end do

  return
end
subroutine power_method ( n, a, y, it_max, tol, lambda, it_num )

!*****************************************************************************80
!
!! POWER_METHOD applies the power method for a real eigenvalue.
!
!  Discussion:
!
!    For a given NxN matrix A and an N vector Y, the power method produces
!    a series of estimates for LAMBDA, the largest eigenvalue, and Y,
!    the eigenvector corresponding to LAMBDA.
!
!    The iteration repeats the following steps
!
!      AY     = A * Y
!      LAMBDA = || AY ||
!      Y      = AY / LAMBDA
!
!    If the matrix A has a single real eigenvalue of maximum modulus,
!    then this iteration will generally produce a good estimate for that
!    eigenvalue and its corresponding eigenvector.
!
!    If there are multiple distinct eigenvalues of the same modulus,
!    perhaps two values of opposite sign, or complex eigenvalues, then
!    the situation is more complicated.
!
!    Separate issues:
!
!    * when estimating the value of LAMBDA, we use the Rayleigh quotient,
!    LAMBDA = ( y' * A * y ) / ( y' * y ).  Since we normalize Y, the
!    bottom of the fraction is 1.  Using this estimate allows us to
!    easily capture the sign of LAMDBA.  Using the eucldean norm
!    instead, for instance, would always give a positive value.
!
!    * If the dominant eigenvalue is negative, then the iteration
!    as given will produce eigenvector iterates that alternate in sign.
!
!    * It is worth knowing whether the successive eigenvector estimates
!    are tending to some value.  Since an eigenvector is really a direction,
!    we need to normalize the vectors, and we need to somehow treat both
!    a vector and its negative as holding the same information.  This
!    means that the proper way of measuring the difference between two
!    eigenvector estimates is to normalize them both, and then compute
!    the cosine between them as y1'y2, followed by the sine, which is
!    sqrt ( 1 - ( y1'y2)^2 ).  If this sine is small, the vectors y1 and y2
!    are "close" in the sense of direction.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Input/output, real ( kind = 8 ) Y(N), the estimate for the eigenvector.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations to take.
!    1 <= IT_MAX.
!
!    Input, real ( kind = 8 ) TOL, an error tolerance.
!
!    Output, real ( kind = 8 ) LAMBDA, the estimate for the eigenvalue.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) ay(n)
  real ( kind = 8 ) cos_y1y2
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  real ( kind = 8 ) lambda
  real ( kind = 8 ) lambda_old
  real ( kind = 8 ) sin_y1y2
  real ( kind = 8 ) tol
  real ( kind = 8 ) val_dif
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y_old(n)

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     IT      Lambda          Delta-Lambda    Delta-Y'
    write ( *, '(a)' ) ' '
  end if
!
!  Force Y to be a vector of unit norm.
!
  y(1:n) = y(1:n) / sqrt ( sum ( y(1:n)**2 ) )

  it_num = 0

  y_old(1:n) = y(1:n)
!
!  Compute AY = A*Y.
!
  ay(1:n) = matmul ( a(1:n,1:n), y(1:n) )
!
!  Estimate LAMBDA = (AY,Y)/(Y,Y).
!
  lambda = dot_product ( y(1:n), ay(1:n) )
!
!  Force AY to have unit norm.
!  Replace Y by AY.
!
  y(1:n) = ay(1:n) / sqrt ( sum ( ay(1:n)**2 ) )
!
!  The sign of Y is optional.  If LAMBDA is probably negative,
!  switch sign of new Y to match old one.
!
  if ( lambda < 0.0D+00 ) then
    y(1:n) = - y(1:n)
  end if

  val_dif = 0.0D+00
  cos_y1y2 = dot_product ( y(1:n), y_old(1:n) )
  sin_y1y2 = sqrt ( ( 1.0D+00 - cos_y1y2 ) * ( 1.0D+00 + cos_y1y2 ) )

  if ( debug ) then
    write ( *, '(2x,i5,2x,g14.6,2x,g14.6,2x,g14.6)' ) it_num, lambda, val_dif, sin_y1y2
  end if
!
!  Now repeat these steps in an iteration.
!
  do it_num = 1, it_max

    lambda_old = lambda
    y_old(1:n) = y(1:n)

    ay(1:n) = matmul ( a(1:n,1:n), y(1:n) )
    lambda = dot_product ( y(1:n), ay(1:n) )
    y(1:n) = ay(1:n) / sqrt ( sum ( ay(1:n)**2 ) )
    if ( lambda < 0.0D+00 ) then
      y(1:n) = - y(1:n)
    end if

    val_dif = abs ( lambda - lambda_old )
    cos_y1y2 = dot_product ( y(1:n), y_old(1:n) )
    sin_y1y2 = sqrt ( ( 1.0D+00 - cos_y1y2 ) * ( 1.0D+00 + cos_y1y2 ) )

    if ( debug ) then
      write ( *, '(2x,i5,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        it_num, lambda, val_dif, sin_y1y2
    end if

    if ( val_dif <= tol ) then
      exit
    end if

  end do

  y(1:n) = ay(1:n) / lambda

  return
end
subroutine power_method2 ( n, a, x_init, it_max, tol, lambda, v, it_num )

!*****************************************************************************80
!
!! POWER_METHOD2 applies the power method for possibly complex eigenvalues.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric VanDeVelde,
!    Concurrent Scientific Programming,
!    Springer, 1994,
!    ISBN: 0-387-94195-9,
!    LC: QA76.58.V35.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Input, real ( kind = 8 ) X_INIT(N), the initial estimate for the eigenvector.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations to take.
!    1 <= IT_MAX.
!
!    Input, real ( kind = 8 ) TOL, an error tolerance.
!
!    Output, complex ( kind = 8 ) LAMBDA, the estimate for the eigenvalue.
!
!    Output, complex ( kind = 8 ) V(N), the estimate for the eigenvector.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  complex ( kind = 8 ) lambda
  real ( kind = 8 ) lambda_imag
  real ( kind = 8 ) lambda_real
  real ( kind = 8 ) pi_xx
  real ( kind = 8 ) pi_xy
  real ( kind = 8 ) pi_xz
  real ( kind = 8 ) pi_yy
  real ( kind = 8 ) pi_yz
  real ( kind = 8 ) pi_zz
  real ( kind = 8 ) tol
  complex ( kind = 8 ) v(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_init(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  it_num = 0
!
!  Compute data necessary to start the iteration.
!
  x(1:n) = x_init(1:n)

  pi_xx = dot_product ( x(1:n), x(1:n) )
  x(1:n) = x(1:n) / pi_xx
  y(1:n) = matmul ( a(1:n,1:n), x(1:n) )
  pi_xy = dot_product ( x(1:n), y(1:n) )
  pi_yy = dot_product ( y(1:n), y(1:n) )

  do it = 1, it_max

    if ( pi_yy - pi_xy * pi_xy < tol * tol * pi_yy ) then
      lambda = pi_xy
      v(1:n) = y(1:n) / sqrt ( pi_yy )
      return
    end if

    z(1:n) = matmul ( a(1:n,1:n), y(1:n) )

    pi_xz = dot_product ( x(1:n), z(1:n) )
    pi_yz = dot_product ( y(1:n), z(1:n) )
    pi_zz = dot_product ( z(1:n), z(1:n) )

    alpha = - ( pi_yz - pi_xy * pi_xz ) / ( pi_yy - pi_xy * pi_xy )
    beta = ( pi_xy * pi_yz - pi_yy * pi_xz ) / ( pi_yy - pi_xy * pi_xy )
    gamma = pi_zz + alpha * alpha * pi_yy + beta * beta &
      + 2.0D+00 * ( alpha * pi_yz + beta * pi_xz + alpha * beta * pi_xy )

    if ( gamma < tol * tol * pi_zz .and. alpha * alpha < 4.0D+00 * beta ) then

      lambda_real = - alpha / 2.0D+00
      lambda_imag = sqrt ( 4.0D+00 * beta - alpha * alpha ) / 2.0D+00
      lambda = complex ( lambda_real, lambda_imag )

      v(1:n) = ( lambda * y(1:n) - z(1:n) ) &
       / sqrt ( beta * pi_yy + alpha * pi_yz + pi_zz )

      return
    end if

    x(1:n) = y(1:n) / sqrt ( pi_yy )
    y(1:n) = z(1:n) / sqrt ( pi_yy )

    pi_xy = pi_yz / pi_yy
    pi_yy = pi_zz / pi_yy

    it_num = it

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POWER_METHOD2 - Fatal error!'
  write ( *, '(a)' ) '  Convergence was not reached.'

  stop
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tris ( m, n, x, y, z, a )

!*****************************************************************************80
!
!! TRIS returns the TRIS matrix.
!
!  Discussion:
!
!    The matrix is a tridiagonal matrix defined by three scalars.
!
!    See page 155 of the Todd reference.
!
!  Formula:
!
!    if ( J = I-1 )
!      A(I,J) = X
!    else if ( J = I )
!      A(I,J) = Y
!    else if ( J = I + 1 )
!      A(I,J) = Z
!    else
!      A(I,J) = 0
!
!  Example:
!
!    M = 5, N = 5, X = 1, Y = 2, Z = 3
!
!    2 3 0 0 0
!    1 2 3 0 0
!    0 1 2 3 0
!    0 0 1 2 3
!    0 0 0 1 2
!
!  Properties:
!
!    A is generally not symmetric: A' /= A.
!
!    A is tridiagonal.
!
!    Because A is tridiagonal, it has property A (bipartite).
!
!    A is banded, with bandwidth 3.
!
!    A is Toeplitz: constant along diagonals.
!
!    If Y is not zero, then for A to be singular, it must be the case that
!
!      0.5 * Y / sqrt ( X * Z ) < 1
!
!    and
!
!      cos (K*PI/(N+1)) = - 0.5 * Y / sqrt ( X * Z ) for some 1 <= K <= N.
!
!    If Y is zero, then A is singular when N is odd, or if X or Z is zero.
!
!    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
!
!    A has eigenvalues
!
!      LAMBDA(I) = Y + 2 * sqrt(X*Z) * COS(I*PI/(N+1))
!
!    The eigenvalues will be complex if X * Z < 0.
!
!    If X = Z, the matrix is symmetric.
!
!    As long as X and Z are nonzero, the matrix is irreducible.
!
!    If X = Z = -1, and Y = 2, the matrix is a symmetric, positive
!    definite M matrix, the negative of the second difference matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Todd,
!    Basic Numerical Mathematics,
!    Volume 2: Numerical Algebra,
!    Birkhauser, 1980,
!    ISBN: 0817608117,
!    LC: QA297.T58.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) X, Y, Z, the scalars that define A.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j = 1, n
    do i = 1, m

      if ( j == i - 1 ) then
        a(i,j) = x
      else if ( j == i ) then
        a(i,j) = y
      else if ( j == i + 1 ) then
        a(i,j) = z
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine tris_eigenvalues ( n, x, y, z, lambda )

!*****************************************************************************80
!
!! TRIS_EIGENVALUES returns the eigenvalues of the TRIS matrix.
!
!  Discussion:
!
!    The eigenvalues will be complex if X * Z < 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) X, Y, Z, the scalars that define A.
!
!    Output, complex ( kind = 8 ) LAMBDA(N), the eigenvalues.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  complex ( kind = 8 ) arg
  integer ( kind = 4 ) i
  complex ( kind = 8 ) lambda(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do i = 1, n
    angle = real ( i, kind = 8 ) * pi / real ( n + 1, kind = 8 )
    arg = cmplx ( x * z, 0.0D+00, kind = 8 )
    lambda(i) = y + 2.0D+00 * sqrt ( arg ) * cos ( angle )
  end do

  return
end
