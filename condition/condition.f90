subroutine combin ( alpha, beta, n, a )

!*****************************************************************************80
!
!! COMBIN returns the COMBIN matrix.
!
!  Discussion:
!
!    This matrix is known as the combinatorial matrix.
!
!  Formula:
!
!    If ( I = J ) then
!      A(I,J) = ALPHA + BETA
!    else
!      A(I,J) = BETA
!
!  Example:
!
!    N = 5, ALPHA = 2, BETA = 3
!
!    5 3 3 3 3
!    3 5 3 3 3
!    3 3 5 3 3
!    3 3 3 5 3
!    3 3 3 3 5
!
!  Properties:
!
!    A is symmetric: A' = A.
!
!    Because A is symmetric, it is normal.
!
!    Because A is normal, it is diagonalizable.
!
!    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
!
!    A is a circulant matrix: each row is shifted once to get the next row.
!
!    det ( A ) = ALPHA^(N-1) * ( ALPHA + N * BETA ).
!
!    A has constant row sums.
!
!    Because A has constant row sums,
!    it has an eigenvalue with this value,
!    and a (right) eigenvector of ( 1, 1, 1, ..., 1 ).
!
!    A has constant column sums.
!
!    Because A has constant column sums,
!    it has an eigenvalue with this value,
!    and a (left) eigenvector of ( 1, 1, 1, ..., 1 ).
!
!    LAMBDA(1:N-1) = ALPHA,
!    LAMBDA(N) = ALPHA + N * BETA.
!
!    The eigenvector associated with LAMBDA(N) is (1,1,1,...,1)/sqrt(N).
!
!    The other N-1 eigenvectors are simply any (orthonormal) basis
!    for the space perpendicular to (1,1,1,...,1).
!
!    A is nonsingular if ALPHA /= 0 and ALPHA + N * BETA /= 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Gregory, David Karney,
!    A Collection of Matrices for Testing Computational Algorithms,
!    Wiley, 1969,
!    ISBN: 0882756494,
!    LC: QA263.68
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms, Second Edition,
!    Addison-Wesley, Reading, Massachusetts, 1973, page 36.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, BETA, scalars that define A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i

  a(1:n,1:n) = beta

  do i = 1, n
    a(i,i) = a(i,i) + alpha
  end do

  return
end
subroutine combin_inverse ( alpha, beta, n, a )

!*****************************************************************************80
!
!! COMBIN_INVERSE returns the inverse of the COMBIN matrix.
!
!  Formula:
!
!    if ( I = J )
!      A(I,J) = (ALPHA+(N-1)*BETA) / (ALPHA*(ALPHA+N*BETA))
!    else
!      A(I,J) =             - BETA / (ALPHA*(ALPHA+N*BETA))
!
!  Example:
!
!    N = 5, ALPHA = 2, BETA = 3
!
!           14 -3 -3 -3 -3
!           -3 14 -3 -3 -3
!   1/34 *  -3 -3 14 -3 -3
!           -3 -3 -3 14 -3
!           -3 -3 -3 -3 14
!
!  Properties:
!
!    A is symmetric: A' = A.
!
!    Because A is symmetric, it is normal.
!
!    Because A is normal, it is diagonalizable.
!
!    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
!
!    A is a circulant matrix: each row is shifted once to get the next row.
!
!    A is Toeplitz: constant along diagonals.
!
!    det ( A ) = 1 / (ALPHA^(N-1) * (ALPHA+N*BETA)).
!
!    A is well defined if ALPHA /= 0D+00 and ALPHA+N*BETA /= 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms, Second Edition,
!    Addison-Wesley, Reading, Massachusetts, 1973, page 36.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, BETA, scalars that define the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) bot
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( alpha == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMBIN_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  The entries of the matrix are undefined '
    write ( *, '(a)' ) '  because ALPHA = 0.'
    stop
  else if ( alpha + n * beta == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMBIN_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  The entries of the matrix are undefined '
    write ( *, '(a)' ) '  because ALPHA+N*BETA is zero.'
    stop
  end if

  bot = alpha * ( alpha + real ( n, kind = 8 ) * beta )

  do j = 1, n
    do i = 1, n

      if ( i == j ) then
        a(i,j) = ( alpha + real ( n - 1, kind = 8 ) * beta ) / bot
      else
        a(i,j) = - beta / bot
      end if

    end do
  end do

  return
end
subroutine condition_hager ( n, a, cond )

!*****************************************************************************80
!
!! CONDITION_HAGER estimates the L1 condition number of a matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Hager,
!    Condition Estimates,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, June 1984, pages 311-316.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) COND, the estimated L1 condition.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) cond
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) info
  integer ( kind = 4) job
  integer ( kind = 4 ) pivot(n)
  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) r8mat_norm_l1

  i1 = -1
  c1 = 0.0D+00
!
!  Factor the matrix.
!
  a_lu(1:n,1:n) = a(1:n,1:n)
  call r8ge_fa ( n, a_lu, pivot, info )

  b(1:n) = 1.0D+00 / real ( n, kind = 8 )

  do

    job = 0
    call r8ge_sl ( n, a_lu, pivot, b, job )

    c2 = sum ( abs ( b(1:n) ) )

    do i = 1, n
      b(i) = r8_sign ( b(i) )
    end do

    job = 1
    call r8ge_sl ( n, a_lu, pivot, b, job )

    call r8vec_max_abs_index ( n, b, i2 )

    if ( 1 <= i1 ) then
      if ( i1 == i2 .or. c2 <= c1 ) then
        exit
      end if
    end if

    i1 = i2
    c1 = c2

    b(1:n) = 0.0D+00
    b(i1) = 1.0D+00

  end do
  
  cond = c2 * r8mat_norm_l1 ( n, n, a )

  return
end
subroutine condition_linpack ( n, a, pivot, cond, z )

!*****************************************************************************80
!
!! CONDITION_LINPACK estimates the L1 condition number of a matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!    For the system A * X = B, relative perturbations in A and B
!    of size EPSILON may cause relative perturbations in X of size
!    EPSILON * COND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 March 2004
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input/output, real ( kind = 8 ) A(N,N).  On input, a matrix to be factored.
!    On output, the LU factorization of the matrix.
!
!    Output, integer ( kind = 4 ) PIVOT(N), the pivot indices.
!
!    Output, real ( kind = 8 ) COND, the estimated L1 condition.
!
!    Output, real ( kind = 8 ) Z(N), a work vector whose contents are 
!    usually unimportant.  If A is close to a singular matrix, then Z is
!    an approximate null vector in the sense that
!      norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) anorm
  real ( kind = 8 ) cond
  real ( kind = 8 ) ek
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) pivot(n)
  real ( kind = 8 ) s
  real ( kind = 8 ) sm
  real ( kind = 8 ) t
  real ( kind = 8 ) wk
  real ( kind = 8 ) wkm
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z(n)
!
!  Compute the L1 norm of A.
!
  anorm = 0.0D+00
  do j = 1, n
    anorm = max ( anorm, sum ( abs ( a(1:n,j) ) ) )
  end do
!
!  Compute the LU factorization.
!
  call r8ge_fa ( n, a, pivot, info )
!
!  COND = norm(A) * (estimate of norm(inverse(A)))
!
!  estimate of norm(inverse(A)) = norm(Z) / norm(Y)
!
!  where
!    A * Z = Y
!  and
!    A' * Y = E
!
!  The components of E are chosen to cause maximum local growth in the
!  elements of W, where U'*W = E.  The vectors are frequently rescaled
!  to avoid overflow.
!
!  Solve U' * W = E.
!
  ek = 1.0D+00
  z(1:n) = 0.0D+00

  do k = 1, n

    if ( z(k) /= 0.0D+00 ) then
      ek = sign ( ek, -z(k) )
    end if

    if ( abs ( a(k,k) ) < abs ( ek - z(k) ) ) then
      s = abs ( a(k,k) ) / abs ( ek - z(k) )
      z(1:n) = s * z(1:n)
      ek = s * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = abs ( wk )
    sm = abs ( wkm )

    if ( a(k,k) /= 0.0D+00 ) then
      wk = wk / a(k,k)
      wkm = wkm / a(k,k)
    else
      wk = 1.0D+00
      wkm = 1.0D+00
    end if

    if ( k + 1 <= n ) then

      do j = k+1, n
        sm = sm + abs ( z(j) + wkm * a(k,j) )
        z(j) = z(j) + wk * a(k,j)
        s = s + abs ( z(j) )
      end do

      if ( s < sm ) then
        t = wkm - wk
        wk = wkm
        z(k+1:n) = z(k+1:n) + t * a(k,k+1:n)
      end if

    end if

    z(k) = wk

  end do

  t = sum ( abs ( z(1:n) ) )
  z(1:n) = z(1:n) / t
!
!  Solve L' * Y = W
!
  do k = n, 1, -1

    z(k) = z(k) + sum ( a(k+1:n,k) * z(k+1:n) )

    t = abs ( z(k) )

    if ( 1.0D+00 < t ) then
      z(1:n) = z(1:n) / t
    end if

    l = pivot(k)

    t    = z(l)
    z(l) = z(k)
    z(k) = t

  end do

  z(1:n) = z(1:n) / sum ( abs ( z(1:n) ) )

  ynorm = 1.0D+00
!
!  Solve L * V = Y.
!
  do k = 1, n

    l = pivot(k)

    t    = z(l)
    z(l) = z(k)
    z(k) = t

    z(k+1:n) = z(k+1:n) + z(k) * a(k+1:n,k)

    if ( 1.0D+00 < abs ( z(k) ) ) then
      ynorm = ynorm / abs ( z(k) )
      z(1:n) = z(1:n) / abs ( z(k) )
    end if

  end do

  s = sum ( abs ( z(1:n) ) )
  z(1:n) = z(1:n) / s
  ynorm = ynorm / s
!
!  Solve U * Z = V.
!
  do k = n, 1, -1

    if ( abs ( a(k,k) ) < abs ( z(k) ) ) then
      s = abs ( a(k,k) ) / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

    if ( a(k,k) /= 0.0D+00 ) then
      z(k) = z(k) / a(k,k)
    else
      z(k) = 1.0D+00
    end if

    z(1:k-1) = z(1:k-1) - z(k) * a(1:k-1,k)

  end do
!
!  Normalize Z in the L1 norm.
!
  s = 1.0D+00 / sum ( abs ( z(1:n) ) )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm

  cond = anorm / ynorm

  return
end
subroutine condition_sample1 ( n, a, m, cond )

!*****************************************************************************80
!
!! CONDITION_SAMPLE1 estimates the L1 condition number of a matrix.
!
!  Discussion:
!
!    A naive sampling method is used.
!
!    Only "forward" sampling is used, that is, we only look at results
!    of the form y=A*x.
!
!    Presumably, solving systems A*y=x would give us a better idea of 
!    the inverse matrix.
!
!    Moreover, a power sequence y1 = A*x, y2 = A*y1, ... and the same for
!    the inverse might work better too.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Input, integer ( kind = 4 ) M, the number of samples to use.
!
!    Output, real ( kind = 8 ) COND, the estimated L1 condition.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a_norm
  real ( kind = 8 ) ainv_norm
  real ( kind = 8 ) ax(n)
  real ( kind = 8 ) ax_norm
  real ( kind = 8 ) cond
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) r8vec_norm_l1
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_norm

  a_norm = 0.0D+00
  ainv_norm = 0.0D+00
  seed = 123456789

  do i = 1, m

    call r8vec_uniform_unit ( n, seed, x )
    x_norm = r8vec_norm_l1 ( n, x )
    ax = matmul ( a, x )
    ax_norm = r8vec_norm_l1 ( n, ax )

    if ( ax_norm == 0.0D+00 ) then
      cond = 0.0D+00
      return
    end if

    a_norm    = max ( a_norm,    ax_norm / x_norm  )
    ainv_norm = max ( ainv_norm, x_norm  / ax_norm )

  end do

  cond = a_norm * ainv_norm

  return
end
subroutine conex1 ( alpha, a )

!*****************************************************************************80
!
!! CONEX1 returns the CONEX1 matrix.
!
!  Discussion:
!
!    The CONEX1 matrix is a counterexample to the LINPACK condition
!    number estimator RCOND available in the LINPACK routine DGECO.
!
!  Formula:
!
!    1  -1 -2*ALPHA   0
!    0   1    ALPHA    -ALPHA
!    0   1  1+ALPHA  -1-ALPHA
!    0   0  0           ALPHA
!
!  Example:
!
!    ALPHA = 100
!
!    1  -1  -200     0
!    0   1   100  -100
!    0   1   101  -101
!    0   0     0   100
!
!  Properties:
!
!    A is generally not symmetric: A' /= A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Cline, Russell Rew,
!    A set of counterexamples to three condition number estimators,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 4, Number 4, December 1983, pages 602-611.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the scalar defining A.  
!    A common value is 100.0.
!
!    Output, real ( kind = 8 ) A(4,4), the matrix.
!
  implicit none

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ) alpha

  a(1,1) = 1.0D+00
  a(2,1) = 0.0D+00
  a(3,1) = 0.0D+00
  a(4,1) = 0.0D+00

  a(1,2) = -1.0D+00
  a(2,2) = 1.0D+00
  a(3,2) = 1.0D+00
  a(4,2) = 0.0D+00

  a(1,3) = -2.0D+00 * alpha
  a(2,3) = alpha
  a(3,3) = 1.0D+00 + alpha
  a(4,3) = 0.0D+00

  a(1,4) = 0.0D+00
  a(2,4) = -alpha
  a(3,4) = -1.0D+00 - alpha
  a(4,4) = alpha

  return
end
subroutine conex1_inverse ( alpha, a )

!*****************************************************************************80
!
!! CONEX1_INVERSE returns the inverse of the CONEX1 matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the scalar defining A.  
!
!    Output, real ( kind = 8 ) A(4,4), the matrix.
!
  implicit none

  real ( kind = 8 )    a(4,4)
  real ( kind = 8 )    alpha

  a(1,1) =  1.0D+00
  a(1,2) =  1.0D+00 - alpha
  a(1,3) =            alpha
  a(1,4) =  2.0D+00

  a(2,1) =  0.0D+00
  a(2,2) =  1.0D+00 + alpha
  a(2,3) =          - alpha
  a(2,4) =  0.0D+00

  a(3,1) =  0.0D+00
  a(3,2) = -1.0D+00
  a(3,3) =  1.0D+00
  a(3,4) =  1.0D+00 / alpha

  a(4,1) = 0.0D+00
  a(4,2) = 0.0D+00
  a(4,3) = 0.0D+00
  a(4,4) = 1.0D+00 / alpha

  return
end
subroutine conex2 ( alpha, a )

!*****************************************************************************80
!
!! CONEX2 returns the CONEX2 matrix.
!
!  Formula:
!
!    1   1-1/ALPHA^2  -2
!    0   1/ALPHA      -1/ALPHA
!    0   0             1
!
!  Example:
!
!    ALPHA = 100
!
!    1  0.9999  -2
!    0  0.01    -0.01
!    0  0        1
!
!  Properties:
!
!    A is generally not symmetric: A' /= A.
!
!    A is upper triangular.
!
!    det ( A ) = 1 / ALPHA.
!
!    LAMBDA = ( 1, 1/ALPHA, 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Cline, Russell Rew,
!    A set of counterexamples to three condition number estimators,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 4, Number 4, December 1983, pages 602-611.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the scalar defining A.  
!    A common value is 100.0.  ALPHA must not be zero.
!
!    Output, real ( kind = 8 ) A(3,3), the matrix.
!
  implicit none

  real ( kind = 8 ) a(3,3)
  real ( kind = 8 ) alpha

  if ( alpha == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CONEX2 - Fatal error!'
    write ( *, '(a)' ) '  The input value of ALPHA was zero!'
    stop
  end if

  a(1,1) = 1.0D+00
  a(1,2) = ( alpha**2 - 1.0D+00 ) / alpha**2
  a(1,3) = -2.0D+00

  a(2,1) = 0.0D+00
  a(2,2) = 1.0D+00 / alpha
  a(2,3) = -1.0D+00 / alpha

  a(3,1) = 0.0D+00
  a(3,2) = 0.0D+00
  a(3,3) = 1.0D+00

  return
end
subroutine conex2_inverse ( alpha, a )

!*****************************************************************************80
!
!! CONEX2_INVERSE returns the inverse of the CONEX2 matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the scalar defining A.  
!    A common value is 100.0.  ALPHA must not be zero.
!
!    Output, real ( kind = 8 ) A(3,3), the matrix.
!
  implicit none

  real ( kind = 8 ) a(3,3)
  real ( kind = 8 ) alpha

  if ( alpha == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CONEX2_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  The input value of ALPHA was zero!'
    stop
  end if

  a(1,1) = 1.0D+00
  a(1,2) = ( 1.0D+00 - alpha**2 ) / alpha
  a(1,3) = ( 1.0D+00 + alpha**2 ) / alpha**2

  a(2,1) = 0.0D+00
  a(2,2) = alpha
  a(2,3) = 1.0D+00

  a(3,1) = 0.0D+00
  a(3,2) = 0.0D+00
  a(3,3) = 1.0D+00

  return
end
subroutine conex3 ( n, a )

!*****************************************************************************80
!
!! CONEX3 returns the CONEX3 matrix.
!
!  Formula:
!
!    if ( I = J and I < N )
!      A(I,J) =  1.0D+00 for 1<=I<N
!    else if ( I = J = N )
!      A(I,J) = -1.0D+00
!    else if ( J < I )
!      A(I,J) = -1.0D+00
!    else
!      A(I,J) =  0.0D+00
!
!  Example:
!
!    N = 5
!
!     1  0  0  0  0
!    -1  1  0  0  0
!    -1 -1  1  0  0
!    -1 -1 -1  1  0
!    -1 -1 -1 -1 -1
!
!  Properties:
!
!    A is generally not symmetric: A' /= A.
!
!    A is integral, therefore det ( A ) is integral, and 
!    det ( A ) * inverse ( A ) is integral.
!
!    A is lower triangular.
!
!    det ( A ) = -1.
!
!    A is unimodular.
!
!    LAMBDA = ( 1, 1, 1, 1, -1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Cline, Russell Rew,
!    A set of counterexamples to three condition number estimators,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 4, Number 4, December 1983, pages 602-611.
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

  do j = 1, n
    do i = 1, n

      if ( j < i ) then
        a(i,j) = -1.0D+00
      else if ( j == i .and. i /= n ) then
        a(i,j) = 1.0D+00
      else if ( j == i .and. i == n ) then
        a(i,j) = - 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine conex3_inverse ( n, a )

!*****************************************************************************80
!
!! CONEX3_INVERSE returns the inverse of the CONEX3 matrix.
!
!  Example:
!
!    N = 5
!
!     1  0  0  0  0
!     1  1  0  0  0
!     2  1  1  0  0
!     4  2  1  1  0
!    -8 -4 -2 -1 -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Cline, Russell Rew,
!    A set of counterexamples to three condition number estimators,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 4, Number 4, December 1983, pages 602-611.
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

  do j = 1, n
    do i = 1, n

      if ( i < n ) then
      
        if ( j < i ) then
          a(i,j) = 2.0D+00**( i - j - 1 )
        else if ( i == j ) then
          a(i,j) = 1.0D+00
        else
          a(i,j) = 0.0D+00
        end if
	
      else if ( i == n ) then
      
        if ( j < i ) then
          a(i,j) = - 2.0D+00**( i - j - 1 )
        else
          a(i,j) = -1.0D+00
        end if
	
      end if
      
    end do
  end do

  return
end
subroutine conex4 ( a )

!*****************************************************************************80
!
!! CONEX4 returns the CONEX4 matrix.
!
!  Discussion:
!
!    7  10   8   7
!    6   8  10   9
!    5   7   9  10
!    5   7   6   5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A(4,4), the matrix.
!
  implicit none

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ), dimension ( 4, 4 ) :: a_save = reshape ( (/ &
     7.0D+00,  6.0D+00,  5.0D+00,  5.0D+00, &
    10.0D+00,  8.0D+00,  7.0D+00,  7.0D+00, &
     8.0D+00, 10.0D+00,  9.0D+00,  6.0D+00, &
     7.0D+00,  9.0D+00, 10.0D+00,  5.0D+00 /), (/ 4, 4 /) )

  a(1:4,1:4) = a_save(1:4,1:4)

  return
end
subroutine conex4_inverse ( a )

!*****************************************************************************80
!
!! CONEX4_INVERSE returns the inverse CONEX4 matrix.
!
!  Discussion:
!
!   -41  -17   10   68
!    25   10   -6  -41
!    10    5   -3  -17
!    -6   -3    2   10
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A(4,4), the matrix.
!
  implicit none

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ), dimension ( 4, 4 ) :: a_save = reshape ( (/ &
   -41.0D+00,  25.0D+00,  10.0D+00, -6.0D+00, &
   -17.0D+00,  10.0D+00,   5.0D+00, -3.0D+00, &
    10.0D+00,  -6.0D+00,  -3.0D+00,  2.0D+00, &
    68.0D+00, -41.0D+00, -17.0D+00, 10.0D+00 /), (/ 4, 4 /) )

  a(1:4,1:4) = a_save(1:4,1:4)

  return
end
subroutine kahan ( alpha, m, n, a )

!*****************************************************************************80
!
!! KAHAN returns the KAHAN matrix.
!
!  Formula:
!
!    if ( I = J )
!      A(I,I) =  sin(ALPHA)^I
!    else if ( I < J )
!      A(I,J) = - sin(ALPHA)^I * cos(ALPHA)
!    else
!      A(I,J) = 0
!
!  Example:
!
!    ALPHA = 0.25, N = 4
!
!    S  -C*S    -C*S      -C*S
!    0     S^2  -C*S^2    -C*S^2
!    0     0       S^3    -C*S^3
!    0     0       0         S^4
!
!    where
!
!      S = sin(ALPHA), C=COS(ALPHA)
!
!  Properties:
!
!    A is upper triangular.
!
!    A = B * C, where B is a diagonal matrix and C is unit upper triangular.
!    For instance, for the case M = 3, N = 4:
!
!    A = | S 0    0    |  * | 1 -C -C  -C |
!        | 0 S^2  0    |    | 0  1 -C  -C |
!        | 0 0    S^3  |    | 0  0  1  -C |
!
!    A is generally not symmetric: A' /= A.
!
!    A has some interesting properties regarding estimation of
!    condition and rank.
!
!    det ( A ) = sin(ALPHA)^(N*(N+1)/2).
!
!    LAMBDA(I) = sin ( ALPHA )^I
!
!    A is nonsingular if and only if sin ( ALPHA ) =/= 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Nicholas Higham,
!    A survey of condition number estimation for triangular matrices,
!    SIAM Review,
!    Volume 9, Number 4, December 1987, pages 575-596.
!
!    W Kahan,
!    Numerical Linear Algebra,
!    Canadian Mathematical Bulletin,
!    Volume 9, 1966, pages 757-801.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the scalar that defines A.  A typical
!    value is 1.2.  The "interesting" range of ALPHA is 0 < ALPHA < PI.
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) csi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) si

  do i = 1, m

    si = sin ( alpha )**i
    csi = - cos ( alpha ) * si

    do j = 1, n

      if ( j < i ) then
        a(i,j) = 0.0D+00
      else if ( j == i ) then
        a(i,j) = si
      else
        a(i,j) = csi
      end if

    end do
  end do

  return
end
subroutine kahan_inverse ( alpha, n, a )

!*****************************************************************************80
!
!! KAHAN_INVERSE returns the inverse of the KAHAN matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the scalar that defines A.  A typical 
!    value is 1.2.  The "interesting" range of ALPHA is 0 < ALPHA < PI.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) ci
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) si

  ci = cos ( alpha )

  do i = 1, n
    do j = 1, n

      if ( i == j ) then
        a(i,j) = 1.0D+00
      else if ( i == j - 1 ) then
        a(i,j) = ci
      else if ( i < j ) then
        a(i,j) = ci * ( 1.0D+00 + ci )**( j - i - 1 )
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do
!
!  Scale the columns.
!
  do j = 1, n
    si = sin ( alpha )**j
    a(1:n,j) = a(1:n,j) / si
  end do

  return
end
subroutine r8ge_fa ( n, a, pivot, info )

!*****************************************************************************80
!
!! R8GE_FA performs a LINPACK style PLU factorization of an R8GE matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N,N), the matrix to be factored.
!    On output, A contains an upper triangular matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = L * U, where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Output, integer ( kind = 4 ) PIVOT(N), a vector of pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t

  info = 0

  do k = 1, n - 1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k + 1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      t      = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = -a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k + 1, n

      if ( l /= k ) then
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  pivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop
  end if

  return
end
subroutine r8ge_inverse ( n, a, pivot )

!*****************************************************************************80
!
!! R8GE_INVERSE computes the inverse of a matrix factored by R8GE_FA.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!    R8GE_INVERSE is a simplified standalone version of the LINPACK routine
!    SGEDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, the factor information computed by R8GE_FA.
!    On output, the inverse matrix.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from R8GE_FA.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) work(n)
!
!  Compute Inverse(U).
!
  do k = 1, n

    a(k,k) = 1.0D+00 / a(k,k)
    a(1:k-1,k) = -a(1:k-1,k) * a(k,k)

    do j = k + 1, n

      t = a(k,j)
      a(k,j) = 0.0D+00
      a(1:k,j) = a(1:k,j) + a(1:k,k) * t

    end do

  end do
!
!  Form Inverse(U) * Inverse(L).
!
  do k = n - 1, 1, -1

    work(k+1:n) = a(k+1:n,k)
    a(k+1:n,k) = 0.0D+00

    do j = k + 1, n
      a(1:n,k) = a(1:n,k) + a(1:n,j) * work(j)
    end do

    if ( pivot(k) /= k ) then

      do i = 1, n
        t             = a(i,k)
        a(i,k)        = a(i,pivot(k))
        a(i,pivot(k)) = t
      end do

    end if

  end do

  return
end
subroutine r8ge_sl ( n, a_lu, pivot, b, job )

!*****************************************************************************80
!
!! R8GE_SL solves a system factored by R8GE_FA.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!    R8GE_SL is a simplified version of the LINPACK routine SGESL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A_LU(N,N), the LU factors from R8GE_FA.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from R8GE_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB, specifies the operation.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve PL * Y = B.
!
    do k = 1, n - 1

      l = pivot(k)

      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if

      b(k+1:n) = b(k+1:n) + a_lu(k+1:n,k) * b(k)

    end do
!
!  Solve U * X = Y.
!
    do k = n, 1, -1
      b(k) = b(k) / a_lu(k,k)
      b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
    end do
!
!  Solve A' * X = B.
!
  else
!
!  Solve U' * Y = B.
!
    do k = 1, n
      b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
    end do
!
!  Solve ( PL )' * X = Y.
!
    do k = n - 1, 1, -1

      b(k) = b(k) + sum ( b(k+1:n) * a_lu(k+1:n,k) )

      l = pivot(k)

      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

  return
end
