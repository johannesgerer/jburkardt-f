subroutine assign ( n, pnodes, recur, t, ierr )

!*****************************************************************************80
!
!! ASSIGN generates the polynomial whose roots are the preassigned nodes.
!
!  Discussion:
!
!    This routine generates the initial polynomial T whose roots are the
!    required preassigned nodes.
!
!    It requires a user-supplied routine RECUR(K,C,D,E) which defines the
!    recursion coefficients of the orthogonal polynomials.  The routine takes
!    as input the polynomial index K, and returns coefficients C, D and E
!    such that:
!
!      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
!
!  Modified:
!
!    30 April 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of preassigned nodes.
!
!    Input, real ( kind = 8 ) PNODES(N), the preassigned nodes.
!
!    Input, external RECUR, the user-supplied function which defines the
!    orthogonal polynomials.
!
!    Output, real ( kind = 8 ) T(0:N), the coefficients of the polynomial whose
!    roots define the preassigned nodes of the quadrature rule and expressed
!    as:
!      H0 * Sum ( 0 <= I <= N ) T(I)/HI * P(I,X)
!
!    Output, integer ( kind = 4 ) IERR.
!    * 0, No error detected
!    * 1, The linear system of equations used to generate the polynomial T
!      became singular or very ill-conditioned.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) c0
  real    ( kind = 8 ) d0
  real    ( kind = 8 ) e0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real    ( kind = 8 ) p
  real    ( kind = 8 ) p0
  real    ( kind = 8 ) p1
  real    ( kind = 8 ) pnodes(n)
  external             recur
  real    ( kind = 8 ) t(0:n)
  real    ( kind = 8 ) x

  ierr = 0
!
!  Set up the linear system of equations.
!
  do i = 1, n

    x = pnodes(i)
    p0 = 0.0D+00
    p1 = 1.0D+00
    p = p1

    do j = 1, n
      a(i,j) = p
      call recur ( j, c0, d0, e0 )
      p = ( c0 * x + d0 ) * p1 + e0 * p0
      p0 = p1
      p1 = p
    end do

    b(i) = p

  end do
!
!  Solve the linear system.
!
  call dgefa ( a, n, n, ipvt, info )

  if ( info /= 0 ) then
    iflag = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ASSIGN - Fatal error!'
    write ( *, '(a)' ) '  The linear system is too ill-conditioned.'
    return
  end if

  call dgesl ( a, n, n, ipvt, b, 0 )
!
!  Set T.
!
  t(0:n-1) = - b(1:n)
  t(n) = 1.0D+00
!
!  Weight with moments.
!
  m = 0
  iflag = 1
  call transf ( t, m, n, recur, iflag )

  return
end
subroutine bair ( n, polin, polout, a0, a1, a2, recur, idigit, errval, ifail )

!*****************************************************************************80
!
!! BAIR seeks roots of a polynomial.
!
!  Discussion:
!
!    This function carries out a generalized Bairstow root extraction
!    for the polynomial:
!
!      SUM(I = 0 to N)  POLIN(I) * P(I,X).
!
!    It calculates the root as a quadratic factor:
!
!      A2 * P(2,X) - A1 * P(1,X) - A0 * P(0,X)
!
!    where P(I,X) is a general orthogonal polynomial of degree I.
!
!  Modified:
!
!    27 February 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Gene Golub, Thomas Robertson,
!    A generalized Bairstow Algorithm,
!    Communications of the ACM,
!    Volume 10, Number 6, June 1967, pages 371-373.
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of input polynomial POLIN.
!
!    Input, real ( kind = 8 ) POLIN(0:N), coefficients of the polynomial whose
!    quadratic factor is to be found, i.e.
!      POLIN = SUM(I = 0 to N) POLIN(I) * P(I,X)
!
!    Output, real ( kind = 8 ) POLOUT(0:N-2), coefficients of the deflated
!    polynomial of degree N-2 with the quadratic factor removed, i.e.
!      POLOUT = SUM(I = 0 to N-2) POLOUT(I) * P(I,X)
!
!    Input/output, real ( kind = 8 ) A0, A1, A2, on input, the estimated
!    quadratic factors.  On output, the estimate has been improved.
!
!    Input, external RECUR ( ), the function which defines the orthogonal
!    polynomials.  See EXTEND for full description.
!
!    Input, integer ( kind = 4 ) IDIGIT, the node convergence parameter, an
!    integer greater than 0.  An attempt is made to calculate the nodes to
!    the maximum accuracy possible by the machine precision available.
!    IDIGIT controls the assessment procedure to take account of
!    round-off errors and specifies the number of least signific
!    decimal digits that can be ignored (i.e. attributed
!    to round-off) in the computed relative error.  A typical value is 5.
!
!    Output, real ( kind = 8 ) ERRVAL, the mean value of the correction to
!    the coefficients of the quadratic factor. May be used as a measure of the
!    root accuracy when convergence is not achieved.
!
!    Output, integer ( kind = 4 ) IFAIL, error flag.
!    * 0, Quadratic factor found.
!    * 1, Convergence not achieved after 50 iterations.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a
  real    ( kind = 8 ) a0
  real    ( kind = 8 ) a1
  real    ( kind = 8 ) a2
  real    ( kind = 8 ) aa
  real    ( kind = 8 ) ab
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) b
  real    ( kind = 8 ) ba
  real    ( kind = 8 ) bb
  real    ( kind = 8 ) beta
  real    ( kind = 8 ) delta
  real    ( kind = 8 ) eps
  real    ( kind = 8 ) errval
  real    ( kind = 8 ) f1
  real    ( kind = 8 ) f2
  integer ( kind = 4 ) idigit
  integer ( kind = 4 ) ifail
  integer ( kind = 4 ) iter
  real    ( kind = 8 ) polin(0:n)
  real    ( kind = 8 ) polout(0:n-2)
  external             recur
  real    ( kind = 8 ) scale
  real    ( kind = 8 ) tol

  ifail = 0
  iter = 50
!
!  Temporary
!
  iter = 100
  errval = 0.0D+00
!
!  Special cases.
!
  if ( n == 1 ) then
    a0 = - polin(0)
    a1 = - polin(1)
    a2 = 0.0D+00
    return
  end if

  if ( n == 2 ) then
    a0 = - polin(0)
    a1 = - polin(1)
    a2 =  polin(2)
    return
  end if
!
!  Estimated ALPHA and BETA.
!
  tol = 10.0D+00**( - max ( 1, idigit ) )
  alpha = a1 / a2
  beta = a0 / a2

  do

    iter = iter - 1

    if ( iter < 0 ) then
      ifail = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BAIR - Warning!'
      write ( *, '(a)' ) '  Iteration did not meet convergence criterion.'
      write ( *, '(a,g14.6)' ) 'TOL = ', tol
      write ( *, '(a,g14.6)' ) 'TOL * ( abs ( ALPHA ) + 1 ) = ', &
        tol * ( abs ( alpha ) + 1 )
      write ( *, '(a,g14.6)' ) '  abs ( DELTA) = ', &
        abs ( delta )
      write ( *, '(a,g14.6)' ) 'TOL * ( abs ( BETA ) + 1 ) = ', &
        tol * ( abs ( beta ) + 1 )
      write ( *, '(a,g14.6)' ) '  abs ( EPS ) = ', &
        abs ( eps )
      a0 = beta
      a1 = alpha
      a2 = 1.0D+00
      errval = 0.5D+00 * ( abs ( eps ) + abs ( delta ) )
      return
    end if

    call qfact ( n, polin, polout, recur, alpha, beta, a, b, aa, &
      ab, ba, bb )

    scale = max ( abs ( ab ), abs ( bb ) )
    f1 = ab / scale
    f2 = bb / scale
    delta = ( b * f1 - a * f2 ) / ( aa * f2 - ba * f1 )
    scale = max ( abs ( ba ), abs ( aa ) )
    f1 = ba / scale
    f2 = aa / scale
    eps = ( a * f1 - b * f2 ) / ( bb * f2 - ab * f1 )
    alpha = alpha + delta
    beta = beta + eps
!
!  Test for convergence.
!  Stop if correction is less than 1/TOL times the smallest machine
!  relative error.
!
    if ( abs ( delta ) <= tol * ( abs ( alpha ) + 1.0D+00 ) .and. &
         abs ( eps )   <= tol * ( abs ( beta )  + 1.0D+00 ) ) then
      exit
    end if
!
!  Suppressed:
!
!   if ( abs ( alpha ) + tol * abs ( delta ) == abs ( alpha ) .and. &
!        abs ( beta ) + tol * abs ( eps ) == abs ( beta ) ) then
!     exit
!   end if

  end do
!
!  Final iteration to tidy up result.
!
  call qfact ( n, polin, polout, recur, alpha, beta, a, b, aa, ab, ba, bb )

  scale = max ( abs ( ab ), abs ( bb ) )
  f1 = ab / scale
  f2 = bb / scale
  delta = ( b * f1 - a * f2 ) / ( aa * f2 - ba * f1 )
  scale = max ( abs ( ba ), abs ( aa ) )
  f1 = ba / scale
  f2 = aa / scale
  eps = ( a * f1 - b * f2 ) / ( bb * f2 - ab * f1 )
  alpha = alpha + delta
  beta = beta + eps

  a0 = beta
  a1 = alpha
  a2 = 1.0D+00
  errval = 0.5D+00 * ( abs ( eps ) + abs ( delta ) )

  return
end
subroutine check ( n, pnodes, wt, k, h0, recur, test )

!*****************************************************************************80
!
!! CHECK tests a computed quadrature rule.
!
!  Purpose:
!
!    This function carries out a test of the given quadrature rule by
!    computing the appropriate integral of
!
!      W(X) * P(K,X) * P(K,X)
!
!    over the region associated with the weight function W(X) and the
!    orthogonal polynomials P(K,X).
!
!  Modified:
!
!    04 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the quadrature rule.
!    1 <= N.
!
!    Input, real ( kind = 8 ) PNODES(N), the nodes.
!
!    Input, real ( kind = 8 ) WT(N), the weights.
!
!    Input, integer ( kind = 4 ) K, the index of the orthogonal polynomial whose
!    weighted square is to be integrated.  0 <= K.
!
!    Input, real ( kind = 8 ) H0, the integral of the orthogonality weight
!    function over the interval of integration.  Note that
!    P(0,X) is arbitrarily taken to be 1.0.
!
!    Input, external RECUR (), the function which defines the
!    orthogonal polynomials.  See EXTEND for a full description.
!
!    Output, real ( kind = 8 ) TEST, the approximate value of the test integral
!    normalized to unity.  Thus, abs(TEST-1) gives a measure of the
!    quality of the calculated rule.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) cj
  real    ( kind = 8 ) dj
  real    ( kind = 8 ) ej
  real    ( kind = 8 ) h0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) p
  real    ( kind = 8 ) p0
  real    ( kind = 8 ) p1
  real    ( kind = 8 ) pnodes(n)
  external             recur
  real    ( kind = 8 ) test
  real    ( kind = 8 ) wt(n)
  real    ( kind = 8 ) x

  if ( k < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHECK - Fatal error!'
    write ( *, '(a)' ) '  Input value K < 0.'
    stop
  end if

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHECK - Fatal error!'
    write ( *, '(a)' ) '  Input value N < 1.'
    stop
  end if

  if ( h0 <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHECK - Fatal error!'
    write ( *, '(a)' ) '  Input value H0 <= 0.0.'
    stop
  end if

  test = 0.0D+00

  do i = 1, n

    p0 = 0.0D+00
    p1 = 1.0D+00

    x = pnodes(i)

    do j = 0, k - 1
      call recur ( j, cj, dj, ej )
      p = ( cj * x + dj ) * p1 + ej * p0
      p0 = p1
      p1 = p
    end do

    test = test + p1 * p1 * wt(i)

  end do
!
!  Normalize.
!
  test = test / h0

  if ( k == 0 ) then
    return
  end if

  j = 0
  call recur ( j, p, p0, p1 )

  do j = 1, k
    call recur ( j, cj, dj, ej )
    p = - p * ej
  end do

  test = test * cj / p

  return
end
subroutine daxpy ( n, da, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DAXPY computes constant times a vector plus a vector.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
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
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in DX and DY.
!
!    Input, real ( kind = 8 ) DA, the multiplier of DX.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of DX.
!
!    Input/output, real ( kind = 8 ) DY(*), the second vector.
!    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    entries of DY.
!
  implicit none

  real    ( kind = 8 ) da
  real    ( kind = 8 ) dx(*)
  real    ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    return
  end if

  if ( da == 0.0D+00 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dy(iy) = dy(iy) + da * dx(ix)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 4 )

    dy(1:m) = dy(1:m) + da * dx(1:m)

    do i = m+1, n, 4
      dy(i  ) = dy(i  ) + da * dx(i  )
      dy(i+1) = dy(i+1) + da * dx(i+1)
      dy(i+2) = dy(i+2) + da * dx(i+2)
      dy(i+3) = dy(i+3) + da * dx(i+3)
    end do

  end if

  return
end
function ddot ( n, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DDOT forms the dot product of two vectors.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
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
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries in DX.
!
!    Input, real ( kind = 8 ) DY(*), the second vector.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    entries in DY.
!
!    Output, real ( kind = 8 ) DDOT, the sum of the product of the
!    corresponding entries of DX and DY.
!
  implicit none

  real    ( kind = 8 ) ddot
  real    ( kind = 8 ) dtemp
  real    ( kind = 8 ) dx(*)
  real    ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  ddot = 0.0D+00
  dtemp = 0.0D+00

  if ( n <= 0 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dtemp = dtemp + dx(ix) * dy(iy)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 5 )

    do i = 1, m
      dtemp = dtemp + dx(i) * dy(i)
    end do

    do i = m+1, n, 5

      dtemp = dtemp + dx(i  ) * dy(i  ) &
                    + dx(i+1) * dy(i+1) &
                    + dx(i+2) * dy(i+2) &
                    + dx(i+3) * dy(i+3) &
                    + dx(i+4) * dy(i+4)
    end do

  end if

  ddot = dtemp

  return
end
subroutine dgefa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! DGEFA factors a real general matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
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
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    On intput, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers used to obtain
!    it.  The factorization can be written A=L*U, where L is a product of
!    permutation and unit lower triangular matrices, and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity indicator.
!    0, normal value.
!    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
!    but it does indicate that DGESL or DGEDI will divide by zero if called.
!    Use RCOND in DGECO for a reliable indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real    ( kind = 8 ) t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find L = pivot index.
!
    l = idamax ( n-k+1, a(k,k), 1 ) + k - 1
    ipvt(k) = l
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      cycle
    end if
!
!  Interchange if necessary.
!
    if ( l /= k ) then
      t = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Compute multipliers.
!
    t = -1.0D+00 / a(k,k)
    call dscal ( n-k, t, a(k+1,k), 1 )
!
!  Row elimination with column indexing.
!
    do j = k+1, n
      t = a(l,j)
      if ( l /= k ) then
        a(l,j) = a(k,j)
        a(k,j) = t
      end if
      call daxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
    end do

  end do

  ipvt(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
  end if

  return
end
subroutine dgesl ( a, lda, n, ipvt, b, job )

!*****************************************************************************80
!
!! DGESL solves a real general linear system A * X = B.
!
!  Discussion:
!
!    DGESL can solve either of the systems A * X = B or A' * X = B.
!
!    The system matrix must have been factored by DGECO or DGEFA.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if DGECO has set 0.0 < RCOND
!    or DGEFA has set INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
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
!    Input, real ( kind = 8 ) A(LDA,N), the output from DGECO or DGEFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGECO or DGEFA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB.
!    0, solve A * X = B;
!    nonzero, solve A' * X = B.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(lda,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real    ( kind = 8 ) ddot
  real    ( kind = 8 ) t
!
!  Solve A * X = B.
!
  if ( job == 0 ) then

    do k = 1, n-1

      l = ipvt(k)
      t = b(l)

      if ( l /= k ) then
        b(l) = b(k)
        b(k) = t
      end if

      call daxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )

    end do

    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      t = -b(k)
      call daxpy ( k-1, t, a(1,k), 1, b(1), 1 )
    end do

  else
!
!  Solve A' * X = B.
!
    do k = 1, n
      t = ddot ( k-1, a(1,k), 1, b(1), 1 )
      b(k) = ( b(k) - t ) / a(k,k)
    end do

    do k = n-1, 1, -1

      b(k) = b(k) + ddot ( n-k, a(k+1,k), 1, b(k+1), 1 )
      l = ipvt(k)

      if ( l /= k ) then
        t = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

  return
end
subroutine dscal ( n, sa, x, incx )

!*****************************************************************************80
!
!! DSCAL scales a vector by a constant.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
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
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) SA, the multiplier.
!
!    Input/output, real ( kind = 8 ) X(*), the vector to be scaled.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real    ( kind = 8 ) sa
  real    ( kind = 8 ) x(*)

  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 5 )

    x(1:m) = sa * x(1:m)

    do i = m+1, n, 5
      x(i)   = sa * x(i)
      x(i+1) = sa * x(i+1)
      x(i+2) = sa * x(i+2)
      x(i+3) = sa * x(i+3)
      x(i+4) = sa * x(i+4)
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      x(ix) = sa * x(ix)
      ix = ix + incx
    end do

  end if

  return
end
subroutine eprod ( n, j, coeff, work, lw, recur, ifail )

!*****************************************************************************80
!
!! EPROD expands a product of two orthogonal polynomials.
!
!  Discussion:
!
!    This function calculates the expansion of a product of two orthogonal
!    polynomials:
!
!      P(N,X) * P(J,X) = SUM (I = N-J to N+J ) COEFF(I) * P(I,X)
!
!    where J must not exceed N.  The orthogonal polynomials are defined
!    by the recurrence relation calculated by RECUR.
!
!    For proper initialization, the function must first be called
!    with J = 0 and the required value of N.  Subsequent calls must be in
!    the order J = 1,2,,,,,N with the appropriate expansion being
!    generated from previous values and returned in COEFF(*).  The
!    coefficients of P(N-J,X),...., P(N+J,X) are stored in the array
!    COEFF(1),...,COEFF(2*J+1).
!
!  Modified:
!
!    02 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest polynomial degree.  Note that
!    after the initial call with J = 0 the value of N in this argument is
!    ignored.
!
!    Input, integer ( kind = 4 ) J, the current product of P(J,X) with P(N,X)
!    to be calculated.  Note that this function must be first called with J = 0
!    and the required largest N.  Subsequent calls must be
!    in the order J = 1,2,..,N.
!
!    Output, real ( kind = 8 ) COEFF(2*J+1), the coefficients of the expansion.
!
!    Workspace, real ( kind = 8 ) WORK(2*J+1,2).
!    The contents of this work area must not be altered between
!    calls by the calling program.
!
!    Input, integer ( kind = 4 ) LW, the leading dimension of WORK.
!
!    Input, external RECUR ( ), the function which defines the
!    orthogonal polynomials.  See EXTEND for a full description.
!
!    Output, integer ( kind = 4 ) IFAIL, error flag.
!    * 0, Result OK.
!    * 1, J exceeds N.
!    * 2, J has not been called sequentially.
!
  implicit none

  integer ( kind = 4 ) j
  integer ( kind = 4 ) lw

  real    ( kind = 8 ) ci
  real    ( kind = 8 ) cj1
  real    ( kind = 8 ) coeff(2*j+1)
  real    ( kind = 8 ) di
  real    ( kind = 8 ) dj1
  real    ( kind = 8 ) ei
  real    ( kind = 8 ) ej1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ifail
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) itop
  integer ( kind = 4 ), save :: ix(2)
  integer ( kind = 4 ) j2
  integer ( kind = 4 ), save :: last
  integer ( kind = 4 ) n
  external             recur
  integer ( kind = 4 ) s
  integer ( kind = 4 ), save :: ss
  real    ( kind = 8 ) work(lw,2)

  ifail = 0
!
!  Initialize.
!
  if ( j == 0 ) then
    ix(1) = 1
    ix(2) = 2
    coeff(1) = 1.0D+00
    work(1,2) = 1.0D+00
    last = 0
    ss = n
    return
  end if

  s = ss
!
!  Check that J does not exceed S value.
!
  if ( s < j ) then
    ifail = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPROD - Fatal error!'
    return
  end if
!
!  Check that J is used sequentially.
!
  if ( last /= j - 1 ) then
    ifail = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPROD - Fatal error!'
    return
  end if

  last = j
  j2 = j + j
  call recur ( j - 1, cj1, dj1, ej1 )

  if ( j == 1 ) then

    coeff(1:j2+1) = 0.0D+00

  else

    do i = 1, j2 - 3
      coeff(i+2) = work(i,ix(1)) * ej1
    end do

    coeff(1)    = 0.0D+00
    coeff(2)    = 0.0D+00
    coeff(j2)   = 0.0D+00
    coeff(j2+1) = 0.0D+00

  end if

  ibot = s - j + 1
  itop = s + j - 1

  do ii = ibot, itop
    i = ii - ibot + 1
    call recur ( ii, ci, di, ei )
    coeff(i+2) = coeff(i+2) + ( work(i,ix(2) ) / ci ) * cj1
    coeff(i+1) = coeff(i+1) + work(i,ix(2)) * ( dj1 - ( cj1 / ci ) * di )
    coeff(i) = coeff(i) - ( work(i,ix(2)) / ci ) * cj1 * ei
  end do

  ii = ix(1)
  ix(1) = ix(2)
  ix(2) = ii

  do i = 1, j2 + 1
    work(i,ix(2)) = coeff(i)
  end do

  return
end
subroutine extend ( n, m, m0, t, recur, symmet, start, pnodes, h0, nexp, &
  idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda, workb, ldb, iflag )

!*****************************************************************************80
!
!! EXTEND extends a quadrature rule by adding new nodes.
!
!  Discussion:
!
!    This function calculates the N+M node quadrature rule composed of
!    N preassigned nodes together with M nodes chosen optimally to
!    achieve algebraic degree of precision of at least N+2*M-1.
!
!    The orthogonal system of polynomials associated with the quadrature
!    weight is defined generally by the recurrence relation specified in the
!    user supplied function RECUR.
!
!  Modified:
!
!    02 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N; on input, the number of preassigned
!    nodes, and the upper limit to the expansion.  On output, if the
!    computation was successful, N is reset to N+M, which is the appropriate
!    next input value in cases where this routine is called iteratively.
!
!    Input, integer ( kind = 4 ) M, the number of nodes to be optimally added.
!
!    Input/ouput, integer ( kind = 4 ) M0; on input, the lower limit to the
!    expansion of T.  This is ignored if START is TRUE.  If the computation is
!    successful, the output value of M0 is reset to M, which is the appropriate
!    next input value in cases where this routine is called iteratively.
!    On output, the lower limit defining the new orthogonal
!    expansion T.  (Set to M).
!
!    Input/output, real ( kind = 8 ) T(0:max(N-M0,M)); on input, the
!    coefficients TI of the polynomial whose roots define the N preassigned
!    nodes of the quadrature rule and expressed as:
!      SUM (I = M0 to N) (TI/HI) * P(I,X)
!    where HI is the integral of W(X) * P(I,X)**2 over the
!    interval for which orthogonality with respect the weight
!    W(X) is defined (moment integrals) and P(I,X) is the
!    orthogonal polynomial of degree I.  Element T(I-M0) holds the
!    value of TI.
!    Note that T is either,
!    (1) provided explicitly,
!    (2) generated automatically from the N preassigned nodes
!        given in PNODES(*) (if START is TRUE.)
!    or,
!    (3) generated from a previous call to the function.
!    This array should be declared to have at least
!    max(N-M0+1,M+1) elements in the calling program.
!    The service function TRANSF can be used to transform
!    the expansion to the required input form if desired
!    with the parameter IFLAG set to 1.
!    On output, the coefficients TI of the new orthogonal  expansion whose
!    roots are the nodes of the extended quadrature (that is, the preassigned
!    nodes plus the extended nodes) and expressed as:
!      SUM (I = M to N+M) (TI/HI) * P(I,X)
!    T(I-M) holds the value of TI.
!    (For definitions see description of input argument T).  This polynomial
!    can be used as input for further extensions.  The service function
!    TRANSF can be used to remove the moment factors from the expansion if
!    desired with the parameter IFLAG set to 0.
!
!    Input, external RECUR ( ), the user supplied function which defines the
!    orthogonal polynomials.  Given K, CALL RECUR(K,C,D,E) gives
!    the coefficients C,D and E such that,
!      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
!    The parameters are defined as follows:
!      K = Index
!      C, D, E = Parameters in the recurrence relation (functions of K)
!
!    Input, logical SYMMET =
!    * FALSE, if no advantage is to be taken of symmetry, if any,
!      about x = 0 in the interval of integration and the
!      orthogonality  weight function. Note that if symmetry in
!      fact does exist setting this parameter to zero will still
!      produce correct results - only efficiency is effected.
!    * TRUE, if the extended rule computations should
!      exploit symmetry about x = 0 in the interval of
!      integration and the orthogonality  weight function.
!      This reduces the size of the system of linear equations
!      determining EXT by a factor of about 2 (see WORKA). If
!      symmetry does not in fact exist erroneous results will be
!      produced.
!
!    Input, logical START,
!    * TRUE, then the polynomial T is generated to have
!      the preassigned nodes (PNODES) as its roots.
!    * FALSE. then the supplied values of the coefficients
!      of T are used directly.
!
!    Input, real ( kind = 8 ) PNODES(N+M), on input, the preassigned nodes.
!    On output, the nodes of the extended quadrature rule made up from the
!    original preassigned nodes and the new optimally extended nodes.  These
!    values can be used in subsequent iterative use.
!
!    Input, real ( kind = 8 ) H0, the integral of the orthogonality weight
!    function over the interval of integration. Zero moment integral.
!
!    Input, integer ( kind = 4 ) NEXP, the largest negative decimal exponent
!    supported on the computer. (Positive number - typical value 38).
!    Weights less than approximately 10**(-NEXP) are set to zero
!    when the Christoffel-Darboux identity is used (N = M).
!    This may be set to INT(LOG10(X1MACH(2))) where X is set to
!    correspond to the appropriate precision in the PORT library.
!
!    Input, integer ( kind = 4 ) IDIGIT, the node convergence parameter,
!    an integer greater than 0.  An attempt is made to calculate the nodes to
!    the maximum accuracy possible by the machine precision available.
!    IDIGIT controls the assessment procedure to take account of
!    round-off errors and specifies the number of least significan
!    decimal digits that can be ignored (i.e. attributed
!    to round-off) in the computed relative error.  A typical
!    value is 5.
!
!    Output, real ( kind = 8 ) WT(N+M), the quadrature weights for
!    the extended rule associated with the nodes in PNODES.
!
!    Output, integer ( kind = 4 ) NODES, the number of extended nodes found.
!    Normally equals M, but see IFAIL.
!
!    Output, real ( kind = 8 ) QR(M), the real parts of the extended nodes.
!
!    Output, real ( kind = 8 ) QI(M), the imaginary parts of the extended
!    nodes (1,..,NODES).  Hopefully these values are zero!
!
!    Output, real ( kind = 8 ) ERR(M), a measure of the relative error in the
!    nodes.  This may be inspected if the convergence error flag has been raised
!    (IFLAG = 3) to decide if the nodes in question are acceptable.  (ERR(*)
!    actually gives the mean last correction to the quadratic factor in the
!    generalized Bairstow root finder (see BAIR).
!
!    Output, real ( kind = 8 ) EXT(M+1), the coefficients of the polynomial
!    whose roots are the  extended nodes (QRS(*),QIS(*)) and expressed as:
!      EXT = SUM (0 <= I <= M) EXT(I) * P(I,X).
!
!    Output, integer ( kind = 4 ) IWORK(max(M,N)), node convergence flags.
!    Elements 1 to NODES give information on the convergence of the roots of the
!    polynomial EXT corresponding to each extended node.
!    * 0, Convergence of I th root satisfactory;
!    * 1, Convergence of I th root unsatisfactory.
!
!    Workspace, real ( kind = 8 ) WORKA(max(M+1,N+1),max(M+1,N+1)).
!    If SYMMET = TRUE, the dimensions can be reduced to max(M/2+1,N) by
!    max(M/2+1,N+1).
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of WORKA.
!
!    Input, real ( kind = 8 ) WORKB(2*M+1,3).
!
!    Input, integer ( kind = 4 ) LDB, the leading dimension of WORKB.
!
!    Output, integer ( kind = 4 ) IFLAG  = 0, No error detected
!    * 1, The linear system of equations defining the polynomial
!      whose roots are the extended nodes became singular or
!      very  ill-conditioned.   (FATAL).
!    * 2, The linear system of equations used to generate the
!      polynomial T when START is TRUE became singular
!      or very ill-conditioned. (FATAL).
!    * 3, Poor convergence has been detected in the calculation
!      of the roots of EXT corresponding to the new
!      nodes or all nodes have not been found (M not equal
!      to NODES). See also ERR(*).
!    * 4, Possible imaginary nodes detected.
!    * 5, Value of N and M incompatible for SYMMET = TRUE.
!      Both cannot be odd. (FATAL)
!    * 6, Test of new quadrature rule has failed.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldb
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) err(m)
  real    ( kind = 8 ) ext(0:m)
  real    ( kind = 8 ) h0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ideg
  integer ( kind = 4 ) idigit
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iwork(lda)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) nexp
  integer ( kind = 4 ) nlast
  integer ( kind = 4 ) nodes
  integer ( kind = 4 ) num
  real    ( kind = 8 ) pnodes(n+m)
  real    ( kind = 8 ) qi(m)
  real    ( kind = 8 ) qr(m)
  external             recur
  integer ( kind = 4 ) s
  logical              start
  logical              symmet
  real    ( kind = 8 ) t(0:n)
  real    ( kind = 8 ) test
  real    ( kind = 8 ), parameter :: tol = 0.0000001D+00
  real    ( kind = 8 ) worka(0:lda-1,0:*)
  real    ( kind = 8 ) workb(0:ldb-1,3)
  real    ( kind = 8 ) wt(n+m)
!
!  Check the size of LDA.
!
  if ( lda < max ( n + 1, m + 1 ) ) then
     iflag = 7
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'EXTEND - Fatal error!'
     write ( *, '(a)' ) '  IFLAG = 7.'
     write ( *, '(a)' ) '  LDA < max ( N + 1, M + 1 ).'
     stop
   end if
!
!  Check the size of LDB.
!
  if ( ldb < 2 * m + 1 ) then
    iflag = 8
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXTEND - Fatal error!'
    write ( *, '(a)' ) '  IFLAG = 8.'
    write ( *, '(a)' ) '  LDB < 2  M + 1.'
    stop
  end if

  iflag = 0
  nodes = 0
  ideg = n + 2 * m - 1
!
!  N and M cannot both be odd for a symmetric rule.
!
  if ( symmet ) then
    if ( mod ( n, 2 ) == 1 .and. mod ( m, 2 ) == 1 ) then
      iflag = 5
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXTEND - Fatal error!'
      write ( *, '(a)' ) '  IFLAG = 5.'
      write ( *, '(a)' ) '  N and M cannot both be odd for a symmetric rule.'
      return
    end if
  end if
!
!  If required, generate the initial T polynomial corresponding to
!  prescribed preassigned nodes.
!
  if ( start .and. n /= 0 ) then
    call assign ( n, pnodes, recur, t, ierr )
    m0 = 0
    if ( ierr /= 0 ) then
      iflag = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXTEND - Fatal error!'
      write ( *, '(a)' ) '  IFLAG = 2.'
      write ( *, '(a)' ) '  Unable to generate the initial T polynomial'
      write ( *, '(a)' ) '  for the preassigned nodes.'
      stop
    end if
  end if

  nlast = n
!
!  Generate extended expansion coefficients and overwrite T.
!
  call gener ( t, m0, n, m, recur, symmet, ext, iwork, worka, lda, &
    workb, ldb, ierr )

  if ( ierr /= 0 ) then
    iflag = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXTEND - Fatal error!'
    write ( *, '(a)' ) '  IFLAG = 1.'
    return
  end if
!
!  Find extended nodes as roots of EXT(*).
!
  call solve ( ext, m, symmet, recur, idigit, qr, qi, nodes, err, iwork, &
    workb, ldb, ierr )

  if ( ierr /= 0 ) then
    iflag = ierr + 2
  end if

  if ( iflag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXTEND - Fatal error!'
    write ( *, '(a,i8)' ) '  IFLAG = ', iflag
    return
  end if
!
!  Accumulate nodes for extended rule.
!
  pnodes(nlast+1:nlast+m) = qr(1:m)
!
!  Reorder.
!
  call rsort ( pnodes, n, 1 )
!
!  Compute weights (only for positive nodes if symmetric).
!
  if ( symmet ) then
    num = ( n + 1 ) / 2
  else
    num = n
  end if

  do i = 1, num
    call weight ( t, m0, n, pnodes(i), recur, h0, nexp, wt(i) )
    if ( symmet ) then
      wt(n-i+1) = wt(i)
    end if
  end do
!
!  Test the new rule.
!
  do k = 0, min ( 4, ideg / 2 )

    call check ( n, pnodes, wt, k, h0, recur, test )

    if ( tol .lt. abs ( test - 1.0D+00 ) ) then
      iflag = 6
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXTEND - Warning!'
      write ( *, '(a)' ) '  Error test not passed.'
      write ( *, '(a)' ) '  IFLAG = 6.'
      return
    end if

  end do

  return
end
subroutine gener ( t, m0, n, m, recur, symmet, ext, iwork, worka, lda, &
  workb, ldb, iflag )

!*****************************************************************************80
!
!! GENER calculates the polynomial defining the optimal new nodes.
!
!  Discussion:
!
!    Given N preassigned quadrature nodes defined as the roots of the
!    polynomial expansion:
!
!      SUM (I = M0 to N) (TI/HI) * P(I,X)
!
!    calculate the polynomial expansion:
!
!      SUM (I = 0 to M) SI * P(I,X)
!
!    whose roots are the M optimal nodes and the new expansion:
!
!      SUM (I = M to N+M) (RI/HI) * P(I,X)
!
!    whose roots are to the (N + M) nodes of the full extended quadrature rule.
!
!  Modified:
!
!    02 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) T(0:N); on input, the coefficients TI of
!    the polynomial whose roots define the N preassigned nodes of the
!    quadrature rule and expressed as:
!      SUM (I = M0 to N) (TI/HI) * P(I,X)
!    where HI is the integral of W(X) * P(I,X)**2 over the
!    interval for which orthogonality with respect the weight
!    W(X) is defined (moment integrals) and P(I,X) is the
!    orthogonal polynomial of degree I. T(I-M0) holds the
!    value of TI. This array should be declared to have at least
!    max(N-M0+1,M+1) elements in the calling program.
!    On output, the coefficients of the new orthogonal
!    expansion whose roots are the nodes of the extended quadrature rule
!    (that is the preassigned nodes plus the extended nodes).
!    It is expressed as:
!      SUM (I = M to N+M) (TI/HI) * P(I,X)
!    where N and M have their original values. T(I-M) holds
!    the value of TI. See input argument of T for definitions.
!
!    Input/output, integer ( kind = 4 ) M0, the lower limit to the expansion of
!    T.  On output, this is updated to correspond with the output value of T.
!
!    Input/output, integer ( kind = 4 ) N, the upper limit to expansion of T.
!    On output, this is updated to correspond with the output value of T.
!
!    Input, integer ( kind = 4 ) M, the number of nodes to be optimally added.
!
!    Input, external RECUR ( ), the user supplied function which defines the
!    orthogonal polynomials.  Given K, CALL RECUR(K,C,D,E) gives
!    the coefficients C,D and E such that,
!      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
!    The parameters are defined as follows:
!      K = Index
!      C, D, E = parameters in the recurrence relation (functions of K).
!
!    Input, SYMMET
!    * FALSE, if no advantage is to be taken of symmetry, if any,
!      about x = 0 in the interval of integration and the
!      orthogonality  weight function. Note that if symmetry in
!      fact does exist setting this parameter to zero will still
!      produce correct results - only efficiency is effected.
!    * TRUE, if the extended rule computations should
!      exploit symmetry about x = 0 in the interval of
!      integration and the orthogonality  weight function.
!      This reduces the size of the system of linear equations
!      determining EXT by a factor of about 2 (see WORKA). If
!      symmetry does not in fact exist erroneous results will be
!      produced.
!
!    Output, real ( kind = 8 ) EXT(0:M), the coefficients of the polynomial
!    whose roots are the  new extended nodes and expressed as:
!      EXT = SUM (I = 0 to M) EXT(I) * P(I,X)
!
!    Workspace, integer ( kind = 4 ) IWORK(M).
!
!    Workspace, real ( kind = 8 ) WORKA(M+1,max(M+1,N+1)).
!    If SYMMET = TRUE, the dimension can be reduced to
!    M/2+1 by max(M/2+1,N/2+1).
!
!    Input, integer ( kind = 4 ) LDA, the number of elements in the leading
!    dimension of WORKA declared in the calling program.
!
!    Input, real ( kind = 8 ) WORKB(2*M+1,3).
!
!    Input, integer ( kind = 4 ) LDB, the number of elements in the leading
!    dimension of WORKB declared in the calling program.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    * 0, No error detected
!    * 1, The linear system of equations defining the polynomial whose roots
!      are the extended nodes became singular or very ill-conditioned.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldb
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) ext(0:m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ifail
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) iref
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) iwork(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m0
  logical              miss
  logical              msodd
  integer ( kind = 4 ) neq
  logical              neven
  integer ( kind = 4 ) nm
  real    ( kind = 8 ) pmax
  external             recur
  integer ( kind = 4 ) s
  logical              symmet
  real    ( kind = 8 ) t(0:n)
  real    ( kind = 8 ) total
  real    ( kind = 8 ) worka(0:lda-1,0:*)
  real    ( kind = 8 ) workb(0:ldb-1,*)

  iflag = 0
!
!  Look for trivial case.
!
  if ( n == 0 ) then

    ext(0:m-1) = 0.0D+00
    ext(m) = 1.0D+00
    t(0) = 1.0D+00
    n = m
    m0 = m
    return

  end if
!
!  General case.
!
  neven = mod ( n, 2 ) == 0
  nm = n + m
!
!  Form matrix.
!
  do s = 0, m

    msodd = mod ( m + s, 2 ) == 1

    if ( neven .and. msodd .and. symmet ) then
      cycle
    end if

    do j = 0, s

      call eprod ( s, j, workb(0,1), workb(0,2), ldb, recur, ifail )

      if ( mod ( n + s + j, 2 ) /= 1 .or. .not. symmet ) then

        iref = s - j
        itop = min ( n, j + s )
        ibot = max ( m0, iref )

        total = 0.0D+00
        do i = ibot, itop
          total = total + t(i-m0) * workb(i-iref,1)
        end do

        if ( .not. symmet ) then

          worka(s,j) = total
          worka(j,s) = total

        else

          if ( neven ) then
            worka(s/2,j/2) = total
            worka(j/2,s/2) = total
          else if ( msodd ) then
            worka(s/2,j/2) = total
          else
            worka(j/2,s/2) = total
          end if

        end if

      end if

    end do

  end do

  if ( symmet ) then
    neq = m / 2
  else
    neq = m
  end if
!
!  Solve for expansion coefficients.
!
  if ( 0 < neq ) then

    call dgefa ( worka, lda, neq, iwork, info )

    if ( info /= 0 ) then
      iflag = 1
      return
    end if

    call dgesl ( worka, lda, neq, iwork, worka(0,neq), 0 )
!
!  Store expansion coefficients.
!
    ext(0:neq-1) = - worka(0:neq-1,neq)
    ext(neq) = 1.0D+00

  end if
!
!  Calculate new T polynomial.
!
  if ( .not. symmet ) then
!
!  Non-symmetric case.
!
    do s = m, nm

      if ( s /= m ) then

        do j = 0, m

          call eprod  (s, j, workb(0,1), workb(0,2), ldb, recur, ifail )

          iref = s - j
          itop = min ( n, j + s )
          ibot = max ( m0, iref )

          total = 0.0D+00
          do i = ibot, itop
            ir = i - iref
            total = total + t(i-m0) * workb(i-iref,1)
          end do

          worka(m,j) = total

        end do

      end if

      worka(m-1,s-m) = dot_product ( worka(m,0:m), ext(0:m) )

    end do
!
!  Overwrite old values of T.
!
    t(0:n) = worka(m-1,0:n)
!
!  Symmetric case.
!
  else

    do s = m, nm

      if ( mod ( n + m + s, 2 ) /= 1 ) then

        do j = 0, m

          call eprod ( s, j, workb(0,1), workb(0,2), ldb, recur, ifail )

          if ( mod ( n + s + j, 2 ) /= 1 ) then

            iref = s - j
            itop = min ( n, j + s )
            ibot = max ( m0, iref )

            total = 0.0D+00
            do i = ibot, itop
              ir = i - iref
              total = total + t(i-m0) * workb(i-iref,1)
            end do

            worka(neq,j/2) = total

         end if

        end do

        total = dot_product ( ext(0:neq), worka(neq,0:neq) )

        if ( 0 <= neq - 1 ) then
          worka(neq-1,(s-m)/2) = total
        end if

      end if

    end do
!
!  Overwrite old values of T in full unsymmetric form.
!
    ic = n / 2
    miss = .true.

    do j = n, 0, -1
      miss = .not. miss
      if ( miss ) then
        t(j) = 0.0D+00
      else
        if ( 0 <= neq - 1 ) then
          t(j) = worka(neq-1,ic)
          ic = ic - 1
        end if
      end if
    end do
!
!  Convert EXT to full unsymmetric form.
!
    workb(m,1) = 1.0D+00
    ic = neq - 1
    miss = .false.

    do j = m - 1, 0, -1
      miss = .not. miss
      if ( miss ) then
        workb(j,1) = 0.0D+00
      else
        workb(j,1) = ext(ic)
        ic = ic - 1
      end if
    end do

    ext(0:m) = workb(0:m,1)

  end if
!
!  Scale the new T polynomial.
!
  pmax = maxval ( abs ( t(0:n) ) )

  t(0:n) = t(0:n) / pmax

  n = nm
  m0 = m

  return
end
function idamax ( n, dx, incx )

!*****************************************************************************80
!
!! IDAMAX indexes the array element of maximum absolute value.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
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
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector to be examined.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive
!    entries of SX.
!
!    Output, integer ( kind = 4 ) IDAMAX, the index of the element of SX of
!    maximum absolute value.
!
  implicit none

  real    ( kind = 8 ) dmax
  real    ( kind = 8 ) dx(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n

  idamax = 0

  if ( n < 1 .or. incx <= 0 ) then
    return
  end if

  idamax = 1

  if ( n == 1 ) then
    return
  end if

  if ( incx == 1 ) then

    dmax = abs ( dx(1) )

    do i = 2, n
      if ( dmax < abs ( dx(i) ) ) then
        idamax = i
        dmax = abs ( dx(i) )
      end if
    end do

  else

    ix = 1
    dmax = abs ( dx(1) )
    ix = ix + incx

    do i = 2, n
      if ( dmax < abs ( dx(ix) ) ) then
        idamax = i
        dmax = abs ( dx(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
subroutine lfact ( gamma, delta, n, xnode, recur, r, dr )

!*****************************************************************************80
!
!! LFACT removes a linear factor from a polynomial expansion.
!
!  Discussion:
!
!    This function removes the linear factor (X-XNODE) from the polynomial
!    expansion:
!
!      SUM ( 0 <= I <= N ) GAMMA(I) * P(I,X)
!
!    to give the quotient:
!
!      SUM ( 0 <= I <= N-1 ) DELTA(I) * P(I,X).
!
!    and the remainder and its derivative at XNODE.
!
!  Modified:
!
!    02 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) GAMMA(0:N), the polynomial from which the factor
!    is to be removed and expressed as:
!      GAMMA = SUM (0 <= I <= N) GAMMA(I) * P(I,X)
!
!    Output, real ( kind = 8 ) DELTA(0:N-1), the quotient polynomial expressed
!    as:
!      DELTA = SUM (0 <= I <= N-1) DELTA(I) * P(I,X)
!
!    Input, integer ( kind = 4 ) N, the degree of GAMMA.
!
!    Input, real ( kind = 8 ) XNODE, the node to be removed.
!
!    Input, external RECUR ( ), the function which defines the orthogonal
!    polynomials.  See EXTEND for a full description.
!
!    Output, real ( kind = 8 ) R, the remainder from division.
!
!    Output, real ( kind = 8 ) DR, the derivative of R with respect to XNODE.
!    (-R/DR is the Newton correction).
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) bk1
  real    ( kind = 8 ) bk2
  real    ( kind = 8 ) ck
  real    ( kind = 8 ) ck1
  real    ( kind = 8 ) ckm1
  real    ( kind = 8 ) dbk1
  real    ( kind = 8 ) dbk2
  real    ( kind = 8 ) delta(0:n-1)
  real    ( kind = 8 ) dk
  real    ( kind = 8 ) dk1
  real    ( kind = 8 ) dkm1
  real    ( kind = 8 ) dr
  real    ( kind = 8 ) ek
  real    ( kind = 8 ) ek1
  real    ( kind = 8 ) ekm1
  real    ( kind = 8 ) gamma(0:n)
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r
  external             recur
  real    ( kind = 8 ) xnode

  bk1 = 0.0D+00
  bk2 = 0.0D+00
  dbk1 = 0.0D+00
  dbk2 = 0.0D+00
  call recur ( n, ck, dk, ek )
  call recur ( n + 1, ck1, dk1, ek1 )

  do k = n, 0, -1

    r = gamma(k) + ( dk + xnode * ck ) * bk1 + ek1 * bk2
    dr = ( dk + xnode * ck ) * dbk1 + ek1 * dbk2 + ck * bk1
    bk2 = bk1
    bk1 = r
    dbk2 = dbk1
    dbk1 = dr

    if ( k /= 0 ) then
      call recur ( k - 1, ckm1, dkm1, ekm1 )
      delta(k-1) = r * ckm1
    end if

    ek1 = ek
    ck = ckm1
    dk = dkm1
    ek = ekm1

  end do

  return
end
subroutine newton ( t, n, xnode, recur, idigit, delta, errval, ifail )

!*****************************************************************************80
!
!! NEWTON applies Newton's method for a single root of a polynomial.
!
!  Discussion:
!
!    This function applies Newton's method to find a single root of the
!    polynomial T:
!
!      T = SUM (0 <= I <= N) T(I) * P(I,X)
!
!    where P(I,X) are the orthogonal polymonials whose recurrence
!    relation is defined by RECUR.
!
!    The value of T is found from the remainder when T is divided
!    by (X-XNODE).  The derivative (of the remainder) is
!    calculated simultaneously.
!
!    The deflated polynomial:
!
!      DELTA = SUM (I = 0 to N-1) DELTA(I) * P(I,X)
!
!    is also computed.
!
!  Modified:
!
!    02 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(0:N), the polynomial whose roots define the
!    nodes of the quadrature rule and expressed as:
!      T = SUM (0 <= I <= N) T(I) * P(I,X)
!
!    Input, integer ( kind = 4 ) N, the degree of the expansion of T.
!
!    Input/output, real ( kind = 8 ) XNODE, on input, a rough estimate for
!    the root.  On output, this estimate has been improved using Newton's
!    method.
!
!    Input, external RECUR ( ), the function which defines the orthogonal
!    polynomials.  See EXTEND for a full description.
!
!    Input, integer ( kind = 4 ) IDIGIT, the node convergence parameter, an
!    integer greater than 0.  An attempt is made to calculate the nodes to the
!    maximum accuracy possible by the machine precision available.
!    IDIGIT controls the assessment procedure to take account of
!    round-off errors and specifies the number of least significan
!    decimal digits that can be ignored, that is, attributed
!    to round-off, in the computed relative error.  A typical value is 5.
!
!    Output, real ( kind = 8 ) DELTA(0:N-1), the coefficients of the deflated
!    polynomial.
!
!    Output, real ( kind = 8 ) ERRVAL, the value of the correction.  May be used
!    as a measure of the root accuracy when convergence is not achieved.
!
!    Output, integer ( kind = 4 ) IFAIL, error flag.
!    * 0, Convergence OK.
!    * 1, Unsatisfactory convergence after 50 iterations.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) delta(0:n-1)
  real    ( kind = 8 ) dr
  real    ( kind = 8 ) eps
  real    ( kind = 8 ) errval
  integer ( kind = 4 ) idigit
  integer ( kind = 4 ) ifail
  integer ( kind = 4 ) iter
  real    ( kind = 8 ) r
  external             recur
  real    ( kind = 8 ) t(0:n)
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) xnode

  iter = 50
  tol = 10.0D+00**( - max ( 1, idigit ) )

  do

    iter = iter - 1

    if ( iter < 0 ) then
      ifail = 1
      errval = abs ( eps )
      return
    end if

    call lfact ( t, delta, n, xnode, recur, r, dr )
    eps = - r / dr
    xnode = xnode + eps

    if ( abs ( xnode ) + tol * abs ( eps ) == abs ( xnode ) ) then
      exit
    end if

  end do
!
!  Final iteration.
!
  call lfact ( t, delta, n, xnode, recur, r, dr )
  eps = - r / dr
  xnode = xnode + eps
  ifail = 0
  errval = abs ( eps )

  return
end
subroutine qfact ( n, gamma, delta, recur, alpha, beta, a, b, aa, ab, ba, bb )

!*****************************************************************************80
!
!! QFACT divides a polynomial by a quadratic factor.
!
!  Discussion:
!
!    This function divides the polynomial:
!
!      sum ( 0 <= I <= N ) GAMMA(I) * P(I,X)
!
!    by the quadratic factor:
!
!      P(2,X) - ALPHA * P(1,X) - BETA * P(0,X)
!
!    giving the quotient:
!
!      sum ( 0 <= I <= N-2 ) DELTA(I) * P(I,X)
!
!    and remainder:
!
!      A * P(1,X) + B * P(0,X)
!
!    where P(I,X) is the orthogonal polynomial of degree I defined by RECUR.
!
!  Modified:
!
!    02 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of GAMMA.
!
!    Input, real ( kind = 8 ) GAMMA(0:N), a polynomial to be divided by
!    a quadratic factor.
!
!    Output, real ( kind = 8 ) DELTA(0:N-2), the quotient polynomial of
!    degree N-2.
!
!    Input, external RECUR(), the function which defines the
!    orthogonal polynomials.  See EXTEND for a full description.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the coefficients of the quadratic
!    factor.
!
!    Output, real ( kind = 8 ) A, B, the remainder coefficients.
!
!    Output, real ( kind = 8 ) AA, the partial of A with respect to ALPHA.
!
!    Output, real ( kind = 8 ) AB, the partial of A with respect to BETA.
!
!    Output, real ( kind = 8 ) BA, the partial of B with respect to ALPHA.
!
!    Output, real ( kind = 8 ) BB, the partial of B with respect to BETA.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a
  real    ( kind = 8 ) aa
  real    ( kind = 8 ) ab
  real    ( kind = 8 ) adn
  real    ( kind = 8 ) adnm1
  real    ( kind = 8 ) adnp1
  real    ( kind = 8 ) adnp2
  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) b
  real    ( kind = 8 ) ba
  real    ( kind = 8 ) bb
  real    ( kind = 8 ) bdn
  real    ( kind = 8 ) bdnm1
  real    ( kind = 8 ) bdnp1
  real    ( kind = 8 ) bdnp2
  real    ( kind = 8 ) beta
  real    ( kind = 8 ) c0
  real    ( kind = 8 ) c1
  real    ( kind = 8 ) c2
  real    ( kind = 8 ) c3
  real    ( kind = 8 ) cf1
  real    ( kind = 8 ) cf2
  real    ( kind = 8 ) cf3
  real    ( kind = 8 ) cf4
  real    ( kind = 8 ) cf5
  real    ( kind = 8 ) ck
  real    ( kind = 8 ) ck1
  real    ( kind = 8 ) ck2
  real    ( kind = 8 ) ck3
  real    ( kind = 8 ) ck4
  real    ( kind = 8 ) d
  real    ( kind = 8 ) d0
  real    ( kind = 8 ) d1
  real    ( kind = 8 ) d2
  real    ( kind = 8 ) d3
  real    ( kind = 8 ) da
  real    ( kind = 8 ) db
  real    ( kind = 8 ) delta(0:n-2)
  real    ( kind = 8 ) dk
  real    ( kind = 8 ) dk1
  real    ( kind = 8 ) dk2
  real    ( kind = 8 ) dk3
  real    ( kind = 8 ) dk4
  real    ( kind = 8 ) dn
  real    ( kind = 8 ) dnm1
  real    ( kind = 8 ) dnp1
  real    ( kind = 8 ) dnp2
  real    ( kind = 8 ) e0
  real    ( kind = 8 ) e1
  real    ( kind = 8 ) e2
  real    ( kind = 8 ) e3
  real    ( kind = 8 ) ek
  real    ( kind = 8 ) ek1
  real    ( kind = 8 ) ek2
  real    ( kind = 8 ) ek3
  real    ( kind = 8 ) ek4
  real    ( kind = 8 ) gamma(0:n)
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r0
  real    ( kind = 8 ) r1
  real    ( kind = 8 ) r2
  external             recur
  real    ( kind = 8 ) rk1
  real    ( kind = 8 ) rk2
  real    ( kind = 8 ) rk3
  real    ( kind = 8 ) rs
  real    ( kind = 8 ) rs0
  real    ( kind = 8 ) rs1
  real    ( kind = 8 ) sk2
  real    ( kind = 8 ) sn1
  real    ( kind = 8 ) sn2
  real    ( kind = 8 ) sn3
  real    ( kind = 8 ) sn4
  real    ( kind = 8 ) tk2
  real    ( kind = 8 ) uk2
  real    ( kind = 8 ) v1
  real    ( kind = 8 ) vk2
  real    ( kind = 8 ) vlk1
  real    ( kind = 8 ) vlk2
  real    ( kind = 8 ) vlk3
  real    ( kind = 8 ) vlk4
  real    ( kind = 8 ) vm1
  real    ( kind = 8 ) vm2
  real    ( kind = 8 ) vmk1
  real    ( kind = 8 ) vmk2
  real    ( kind = 8 ) vmk3
  real    ( kind = 8 ) w0
  real    ( kind = 8 ) w1
  real    ( kind = 8 ) wk2
!
!  Initialize coefficients.
!
  dnp2 = 0.0D+00
  dnp1 = 0.0D+00
  dn   = 0.0D+00
  dnm1 = 0.0D+00
!
!  Partial coefficients wrt ALPHA.
!
  adnp2 = 0.0D+00
  adnp1 = 0.0D+00
  adn   = 0.0D+00
  adnm1 = 0.0D+00
!
!  Partial coefficients wrt BETA.
!
  bdnp2 = 0.0D+00
  bdnp1 = 0.0D+00
  bdn   = 0.0D+00
  bdnm1 = 0.0D+00
!
!  Scaling parameters.
!
  sn1 = 1.0D+00
  sn2 = 1.0D+00
  sn3 = 1.0D+00
  sn4 = 1.0D+00

  call recur ( 0, c0, d0, e0 )
  call recur ( 1, c1, d1, e1 )
  call recur ( 2, c2, d2, e2 )
  call recur ( 3, c3, d3, e3 )

  r0 = - c0 * e1 / c1
  r1 = - c0 * e2 / c2
  r2 = - c0 * e3 / c3

  vm1 = d0 - c0 * d1 / c1
  vm2 = d0 - c0 * d2 / c2

  w0 = - r1 * e1
  w1 = - c1 * r2 * e2 / c2
  v1 = d1 * r1 - c1 * vm2 * e2 / c2 - c1 * r1 * d1 / c1
  k = n - 2

  call recur ( k + 4, ck4, dk4, ek4 )
  call recur ( k + 3, ck3, dk3, ek3 )
  call recur ( k + 2, ck2, dk2, ek2 )
  call recur ( k + 1, ck1, dk1, ek1 )

  vlk4 = c0 / ck3
  vlk3 = c0 / ck2
  vlk2 = c0 / ck1
  rk3 = - c0 * ek4 / ck4
  rk2 = - c0 * ek3 / ck3
  vmk3 = d0 - dk3 * vlk4
  vmk2 = d0 - dk2 * vlk3
!
!  Extract quadratic factor and find partial derivatives
!
  do k = n - 2, 0, - 1

    call recur ( k, ck, dk, ek )
    vlk1 = c0 / ck
    rk1 = - c0 * ek2 / ck2
    vmk1 = d0 - dk1 * vlk2
    sk2 = c1 * vlk1 * vlk2 / c0
    tk2 = vlk2 * ( d1 - c1 * dk2 / ck2 ) + c1 * vmk1 / ck1
    uk2 = d1 * vmk2 + e1 - c1 * vlk3 * ek3 / ck3 &
      - c1 * vmk2 * dk2 / ck2 + c1 * rk1 / ck1
    vk2 = d1 * rk2 - c1 * vmk3 * ek3 / ck3 - c1 * rk2 * dk2 / ck2
    wk2 = - c1 * rk3 * ek3 / ck3
    cf1 = ( alpha * vlk2 - tk2 ) / sn1
    cf2 = ( beta + alpha*  vmk2 - uk2 ) / sn2
    cf3 = ( alpha * rk2 - vk2 ) / sn3
    cf4 = - wk2 / sn4
    rs = gamma(k+2)
    d = rs + cf1 * dnm1 + cf2 * dn + cf3 * dnp1 + cf4 * dnp2
    delta(k) = d / sk2
    da = vlk2 * dnm1 / sn1 + vmk2 * dn / sn2 + rk2 * dnp1 / sn3 &
      + cf1 * adnm1 + cf2 * adn + cf3 * adnp1 + cf4 * adnp2
    db = dn / sn2 + cf1 * bdnm1 + cf2 * bdn + cf3 * bdnp1 + cf4 * bdnp2
!
!  Recycle old values.
!
    sn4 = sn3
    sn3 = sn2
    sn2 = sn1
    sn1 = sk2
    dnp2 = dnp1
    dnp1 = dn
    dn = dnm1
    dnm1 = d
    adnp2 = adnp1
    adnp1 = adn
    adn = adnm1
    adnm1 = da
    bdnp2 = bdnp1
    bdnp1 = bdn
    bdn = bdnm1
    bdnm1 = db
    ck4 = ck3
    ck3 = ck2
    ck2 = ck1
    ck1 = ck
    dk4 = dk3
    dk3 = dk2
    dk2 = dk1
    dk1 = dk
    ek4 = ek3
    ek3 = ek2
    ek2 = ek1
    ek1 = ek
    vlk4 = vlk3
    vlk3 = vlk2
    vlk2 = vlk1
    rk3 = rk2
    rk2 = rk1
    vmk3 = vmk2
    vmk2 = vmk1

  end do

  cf1 = alpha
  cf2 = beta + alpha * vm1 - r1
  cf3 = alpha * r1 - v1
  cf4 = - w1
  cf5 = alpha * r0
  rs0 = gamma(0)
  rs1 = gamma(1)
  dnm1 = dnm1 / sn1
  dn = dn / sn2
  dnp1 = dnp1 / sn3
  dnp2 = dnp2 / sn4
  adnm1 = adnm1 / sn1
  adn = adn / sn2
  adnp1 = adnp1 / sn3
  adnp2 = adnp2 / sn4
  bdnm1 = bdnm1 / sn1
  bdn = bdn / sn2
  bdnp1 = bdnp1 / sn3
  bdnp2 = bdnp2 / sn4
!
!  Remainder.
!
  a = rs1 + cf1 * dnm1 + cf2 * dn + cf3 * dnp1 + cf4 * dnp2
  b = rs0 + beta * dnm1 + cf5 * dn - w0 * dnp1
!
!  Partials.
!
  aa = dnm1 + vm1 * dn + r1 * dnp1 + cf1 * adnm1 + cf2 * adn &
    + cf3 * adnp1 + cf4 * adnp2
  ab = dn + cf1 * bdnm1 + cf2 * bdn + cf3 * bdnp1 + cf4 * bdnp2
  ba = r0 * dn + beta * adnm1 + cf5 * adn - w0 * adnp1
  bb = dnm1 + beta * bdnm1 + cf5 * bdn - w0 * bdnp1

  return
end
subroutine roots ( a0, a1, a2, zreal1, zimag1, zreal2, zimag2, recur, info )

!*****************************************************************************80
!
!! ROOTS calculates roots of a quadratic factor.
!
!  Discussion:
!
!    This function calculates the roots corresponding to the quadratic factor
!
!      A2 * P(2,X) - A1 * P(1,X) - A0 * P(0,X)
!
!    where P(I,X) is a general orthogonal polynomial of degree I
!    defined by the recurrence calculated by RECUR.
!
!  Modified:
!
!    02 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A0, A1, A2, the coefficients of the quadratic
!    factor.
!
!    Output, real ( kind = 8 ) ZREAL1, ZIMAG1, the real and imaginary parts
!    of root 1.
!
!    Output, real ( kind = 8 ) ZREAL2, ZIMAG2, the real and Imaginary parts
!    of root 2.
!
!    Input, external RECUR ( ), the function which defines the orthogonal
!    polynomials.  See EXTEND for full description.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    * 0, two roots found.
!    * 1, one root only (A2 = 0).
!
  implicit none

  real    ( kind = 8 ) a0
  real    ( kind = 8 ) a1
  real    ( kind = 8 ) a2
  real    ( kind = 8 ) aa
  real    ( kind = 8 ) bb
  real    ( kind = 8 ) c0
  real    ( kind = 8 ) c1
  real    ( kind = 8 ) cc
  real    ( kind = 8 ) d0
  real    ( kind = 8 ) d1
  real    ( kind = 8 ) e0
  real    ( kind = 8 ) e1
  integer ( kind = 4 ) info
  integer ( kind = 4 ) k
  external             recur
  real    ( kind = 8 ) z
  real    ( kind = 8 ) zimag1
  real    ( kind = 8 ) zimag2
  real    ( kind = 8 ) zr
  real    ( kind = 8 ) zreal1
  real    ( kind = 8 ) zreal2

  info = 0

  k = 0
  call recur ( k, c0, d0, e0 )

  if ( a2 == 0.0D+00 ) then
    zreal1 = - ( a0 + a1 * d0 ) / a1 / c0
    zreal2 = 0.0D+00
    zimag1 = 0.0D+00
    zimag2 = 0.0D+00
    info = 1
    return
  end if

  k = 1
  call recur ( k, c1, d1, e1 )

  aa = - c0 * c1 * a2
  bb = - a2 * ( c0 * d1 + d0 * c1 ) + c0 * a1
  cc = - d0 * d1 * a2 - e1 * a2 + a0 + a1 * d0
  z = bb * bb - 4.0D+00 * aa * cc
  zr = sqrt ( abs ( z ) )

  if ( 0.0D+00 <= z ) then
    zimag1 = 0.0D+00
    zimag2 = 0.0D+00
    zreal1 = 0.5D+00 * ( - bb - sign ( zr, bb ) ) / aa
    zreal2 = cc / aa / zreal1
  else
    zreal1 = - 0.5D+00 * bb / aa
    zreal2 = zreal1
    zimag1 = 0.5D+00 * zr / aa
    zimag2 = - zimag1
  end if

  return
end
subroutine rsort ( a, n, iflag )

!*****************************************************************************80
!
!! RSORT carries out a simple ripple sort.
!
!  Modified:
!
!    16 February 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(N), the array to be sorted.
!
!    Input, integer ( kind = 4 ) N, the number of elements to be sorted
!
!    Input, integer ( kind = 4 ) IFLAG, determines the sort direction.
!    * 0, for ascending sort;
!    * 1, for descending sort.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  logical              ascend
  logical              done
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real    ( kind = 8 ) val

  ascend = iflag == 0
!
!  Begin scans.
!
  do j = n - 1, 1, - 1

    done = .true.

    do k = 1, j

      if ( ascend ) then
        k1 = k
        k2 = k + 1
      else
        k1 = k + 1
        k2 = k
      end if
!
!  Exchange elements.
!
      if ( a(k2) < a(k1) ) then
        val = a(k1)
        a(k1) = a(k2)
        a(k2) = val
        done = .false.
      end if

    end do

    if ( done ) then
      return
    end if

  end do

  return
end
subroutine solve ( ext, m, symmet, recur, idigit, qr, qi, nodes, err, &
  icheck, work, ldw, ierr )

!*****************************************************************************80
!
!! SOLVE calculates roots of an orthogonal polynomial expansion.
!
!  Discussion:
!
!    This function calculates the roots of the orthogonal polynomial expansion:
!
!      SUM ( 0 <= I <= M ) EXT(I) * P(I,X)
!
!  Modified:
!
!    16 February 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) EXT(0:M), the coefficients of the polynomial whose
!    roots are required (nodes of the quadrature rule) and expressed as:
!      SUM (I = 0 to M) EXT(I) * P(I,X)
!    The recurrence relation for the orthogonal polynomials
!    P(I,X) is defined by RECUR.
!
!    Input, integer ( kind = 4 ) M, the upper limit to expansion EXT
!    (polynomial degree).
!
!    Input, logical SYMMET
!    * FALSE, if no advantage can be taken of symmetry
!      about x = 0 in the interval of integration and the
!      orthogonality weight function.
!    * TRUE, if symmetry exists about x = 0 in the interval of
!      integration and the orthogonality weight function.
!
!    Input, external RECUR ( ), the user supplied function which defines the
!    orthogonal polynomials.  Given K, RECUR ( K, C, D, E ) gives
!    the coefficients C,D and E such that,
!      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
!    The parameters are defined as follows:
!      K = Index
!      C, D, E = parameters in the recurrence relation (functions of K).
!
!    Input, integer ( kind = 4 ) IDIGIT, the node convergence parameter, an
!    integer greater than 0.  An attempt is made to calculate the nodes to the
!    maximum accuracy possible by the machine precision available.
!    IDIGIT controls the assessment procedure to take account of
!    round-off errors and specifies the number of least significan
!    decimal digits that can be ignored (i.e. attributed
!    to round-off) in the computed relative error.  A typical
!    value is 5.
!
!    Output, real ( kind = 8 ) QR(M), the real parts of the roots fo EXT.
!
!    Output, real ( kind = 8 ) QI(M), the imaginary parts of the roots
!    of EXT.  (Hopefully these values are zero!).
!
!    Output, integer ( kind = 4 ) NODES, the number of extended nodes found.
!    Normally equals M but see IERR.
!
!    Output, real ( kind = 8 ) ERR(M), a measure of the relative error in the
!    roots.  This may be inspected if the convergence error flag has been
!    raised (IERR = 2) to decide if the roots in question are acceptable.
!    (ERR(*) actually gives the mean last correction to the quadratic factor
!    in the generalized Bairstow root finder (see BAIR).
!
!    Output, integer ( kind = 4 ) ICHECK(M), root convergence flags.  Elements
!    1 to NODES give information on the convergence of the roots of the
!    polynomial EXT.
!    * 0, Convergence of I th root satisfactory;
!    * 1, Convergence of I th root unsatisfactory.
!
!    Worskpace, real ( kind = 8 ) WORK(M+1,2).
!
!    Input, integer ( kind = 4 ) LDW, the leading dimension of WORK (which must
!    be at least M+1).
!
!    Output, integer ( kind = 4 ) IERR
!    * 0, No error detected
!    * 1, Poor convergence has been detected in the calculation
!      of the roots of EXT or all roots have not
!      been found (M not equal to NODES).  See also ERR.
!    * 2, Possible imaginary nodes detected.
!
  implicit none

  integer ( kind = 4 ) ldw
  integer ( kind = 4 ) m

  real    ( kind = 8 ) a0
  real    ( kind = 8 ) a1
  real    ( kind = 8 ) a2
  real    ( kind = 8 ) c0
  real    ( kind = 8 ) c1
  real    ( kind = 8 ) d0
  real    ( kind = 8 ) d1
  real    ( kind = 8 ) det
  real    ( kind = 8 ) e0
  real    ( kind = 8 ) e1
  real    ( kind = 8 ) err(m)
  real    ( kind = 8 ) errval
  real    ( kind = 8 ) ext(0:m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icheck(m)
  integer ( kind = 4 ) idigit
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifail
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nodes
  integer ( kind = 4 ) nroot
  real    ( kind = 8 ) p1a
  real    ( kind = 8 ) p1b
  real    ( kind = 8 ) p2a
  real    ( kind = 8 ) p2b
  real    ( kind = 8 ) pmax
  real    ( kind = 8 ) qi(m)
  real    ( kind = 8 ) qr(m)
  external             recur
  logical              reset
  real    ( kind = 8 ) rt1
  real    ( kind = 8 ) rt2
  real    ( kind = 8 ) sa0
  real    ( kind = 8 ) sa1
  real    ( kind = 8 ) sfa0
  real    ( kind = 8 ) sfa1
  logical              symmet
  real    ( kind = 8 ), parameter :: vrt1 = 0.0000001D+00
  real    ( kind = 8 ) work(0:ldw-1,2)
  real    ( kind = 8 ) zi1
  real    ( kind = 8 ) zi2
  real    ( kind = 8 ) zr1
  real    ( kind = 8 ) zr2

  nodes = 0
  ierr = 0
!
!  If M is odd, find and remove initial real root using Newton iteration.
!  Set WORK(*,1) to polynomial to be processed.
!
  if ( mod ( m, 2 ) == 1 ) then

    zr1 = vrt1
    call newton ( ext, m, zr1, recur, idigit, work(0,1), errval, ifail )
    nodes = nodes + 1
    icheck(nodes) = ifail
    err(nodes) = errval
    qr(nodes) = zr1
    qi(nodes) = 0.0D+00
    nroot = m - 1

  else

    work(0:m,1) = ext(0:m)
    nroot = m

  end if
!
!  Find the remaining root pairs.
!  Calculate seed approximation for quadratic factor.
!
  if ( nroot /= 0 ) then

    call recur ( 0, c0, d0, e0 )
    call recur ( 1, c1, d1, e1 )
    rt1 = vrt1

    if ( symmet ) then
      rt2 = - rt1
    else
      rt2 = 0.0D+00
    end if

    p1a = c0 * rt1 + d0
    p1b = c0 * rt2 + d0
    p2a = ( c1 * rt1 + d1 ) * p1a + e1
    p2b = ( c1 * rt2 + d1 ) * p1b + e1
    det = c0 * ( rt1 - rt2 )
    sa1 = ( p2a - p2b ) / det
    sa0 = ( p1a * p2b - p1b * p2a ) / det
    reset = .true.
!
!  Alternate approximation which introduces a small complex component.
!
    rt1 = vrt1
    rt2 = vrt1
    sfa1 = ( c0 * d1 + d0 * c1 ) / c0 + 2.0D+00 * c1 * rt1
    sfa0 = d0 * d1 + e1 - d0 * sfa1 &
      - c0 * c1 * ( rt1 * rt1 + rt2 * rt2 )
!
!  IP1 points to the current deflated polynomial.
!
    ip1 = 1
    ip2 = 2

    do

      if ( reset ) then
        a2 = 1.0D+00
        a1 = sa1
        a0 = sa0
        reset = .false.
      end if

      call bair ( nroot, work(0,ip1), work(0,ip2), a0, a1, a2, recur, &
        idigit, errval, ifail )
!
!  On failure, try again with complex components introduced.
!
      if ( ifail /= 0 ) then
        a2 = 1.0D+00
        a1 = sfa1
        a0 = sfa0
        reset = .true.
        call bair ( nroot, work(0,ip1), work(0,ip2), a0, a1, a2, &
          recur, idigit, errval, ifail )
      end if
!
!  Apply Bairstow to full expansion to avoid error accumulation.
!
      call bair ( m, ext, work(0,ip2), a0, a1, a2, recur, idigit, errval, &
        ifail )
!
!  Tidy up the quotient polynomial.
!
      call qfact ( nroot, work(0,ip1), work(0,ip2), recur, a1, a0, zr1, zr1, &
        zr1, zr1, zr1, zr1 )

      call roots ( a0, a1, a2, zr1, zi1, zr2, zi2, recur, info )
!
!  Record node information.
!  If IFAIL is nonzero, then the calculation is going to be rejected.
!
      nodes = nodes + 1
      icheck(nodes) = ifail
      err(nodes) = errval
      qr(nodes) = zr1
      qi(nodes) = zi1

      nodes = nodes + 1
      icheck(nodes) = ifail
      err(nodes) = errval
      qr(nodes) = zr2
      qi(nodes) = zi2

      nroot = nroot - 2
!
!  Make the deflated polynomial current.
!
      i = ip1
      ip1 = ip2
      ip2 = i
!
!  Scale the deflated polynomial.
!
      if ( nroot <= 0 ) then
        exit
      end if

      pmax = maxval ( abs ( work(0:nroot,ip1) ) )
      work(0:nroot,ip1) = work(0:nroot,ip1) / pmax

    end do

  end if
!
!  Check for poor convergence.
!
  i = sum ( icheck(1:nodes) )

  if ( i /= 0 ) then
    ierr = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SOLVE - Warning:'
    write ( *, '(a)' ) '  Poor convergence for some roots.'
    return
  end if

  if ( nodes /= m ) then
    ierr = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SOLVE - Warning:'
    write ( *, '(a)' ) '  NODE /= M.'
    return
  end if
!
!  Look for possible imaginary nodes.
!
  do j = 1, nodes
    if ( qi(j) /= 0.0D+00 ) then
      write ( *, * ) 'QI(J) /= 0'
      write ( *, * ) 'QI(J) = ', qi(j)
      ierr = 2
      return
    end if
  end do

  return
end
subroutine transf ( t, m, n, recur, iflag )

!*****************************************************************************80
!
!! TRANSF scales a polynomial expansion with respect to the moments.
!
!  Discussion:
!
!    This function scales the polynomial expansion:
!
!      SUM (M to N) TI * P(I,X)
!
!    with respect to the moments HI of the orthogonality weight function
!    giving the expansion:
!
!      H0 * SUM (M to N) (TI/HI) * P(I,X)
!
!    or
!
!      (1/H0) * SUM (M to N) (TI*HI) * P(I,X)
!
!    depending on the value of IFLAG.
!
!  Modified:
!
!    27 February 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) T(0:N), the coefficients TI of the
!    polynomial expansion to be scaled and expressed as:
!      SUM (I = M to N) TI * P(I,X)
!    T(I-M) holds the value of TI.  On output, the polynomial has been scaled.
!
!    Input, integer ( kind = 4 ) M, the lower limit to the expansion of T.
!
!    Input, integer ( kind = 4 ) N, the upper limit to the expansion of T.
!
!    Input, external RECUR ( ), the function which defines the orthogonal
!    polynomials.  See EXTEND for a full description.
!
!    Input, integer ( kind = 4 ) IFLAG, the operation to be carried out.
!    * 0, if coefficient TI is to be replaced by TI*(H0/HI).
!    * 1, if coefficient TI is to be replaced by TI*(HI/H0).
!
  implicit none

  real    ( kind = 8 ) ck
  real    ( kind = 8 ) ckm1
  real    ( kind = 8 ) dk
  real    ( kind = 8 ) ek
  real    ( kind = 8 ) h
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  external             recur
  real    ( kind = 8 ) t(0:n)

  h = 1.0D+00

  do k = 0, n

    call recur ( k, ck, dk, ek )

    if ( k /= 0 ) then
      h = - ckm1 / ck * ek * h
    end if

    if ( m <= k ) then
      if ( iflag == 0 ) then
        t(k-m) = t(k-m) / h
      else
        t(k-m) = t(k-m) * h
      end if
    end if

    ckm1 = ck

  end do

  return
end
subroutine weight ( t, m, n, xnode, recur, h0, nexp, wt )

!*****************************************************************************80
!
!! WEIGHT calculates quadrature weights.
!
!  Discussion:
!
!    This function calculates the quadrature weight associated with the node
!    XNODE in the rule whose nodes are defined by the roots of polynomial T.
!
!    The weight is calculated by dividing T by (X-XNODE) to give:
!
!      S(X) = T(X)/(X-XNODE) = SUM (0 to N-1) G(I) * P(I,X).
!
!    S(X) is then divided by (X-XNODE) to give the remainder R.
!
!    The weight is finally given by H0*G(0)/R. If N = M the
!    Christoffel-Darboux identity result is used to reduce extreme
!    cancellation effects at high degree.
!
!  Modified:
!
!    27 February 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Patterson,
!    An algorithm for generating interpolatory quadrature rules of the highest
!    degree of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 123-136.
!
!    Thomas Patterson,
!    Algorithm 672:
!    EXTEND: generation of interpolatory quadrature rules of the highest degree
!    of precision with preassigned nodes for general weight functions,
!    Transactions on Mathematical Software,
!    Volume 15, Number 2, June 1989, pages 137-143.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(0:N), the coefficients TI of the polynomial
!    whose roots define the N preassigned nodes of the quadrature
!    rule and expressed as:
!      SUM (I = M to N) (TI/HI) * P(I,X)
!    where HI is the integral of W(X) * P(I,X)**2 over the
!    interval for which orthogonality with respect the weight
!    W(X) is defined (moment integrals) and P(I,X) is the
!    orthogonal polynomial of degree I. T(I-M) holds the
!    value of TI. This array should be declared to have at least
!    N-M+1 elements in the calling program.
!
!    Input, integer ( kind = 4 ) M, the lower limit to the expansion of T.
!
!    input, integer ( kind = 4 ) N, the upper limit to the expansion of T.
!
!    Input, real ( kind = 8 ) XNODE, the node whose weight is required
!
!    Input, external RECUR ( ), the function which defines the orthogonal
!    polynomials.  See EXTEND for a full description.
!
!    Input, real ( kind = 8 ) H0, th eintegral of the orthogonality weight
!    function over the interval of integration. Zero moment integral.  Note that
!    P(0,X) is arbitrarily taken to be 1.0
!
!    Input, integer ( kind = 4 ) NEXP, the largest negative decimal exponent
!    supported on the computer. (Positive number - typical value 38).
!    Weights less than approximately 10**(-NEXP) are set to zero
!    when the Christoffel-Darboux identity is used (N = M).
!
!    Output, real ( kind = 8 ) WT, the weight associated with XNODE.
!
  implicit none

  real    ( kind = 8 ) bb
  real    ( kind = 8 ) bk1
  real    ( kind = 8 ) bk2
  real    ( kind = 8 ) ck
  real    ( kind = 8 ) ck1
  real    ( kind = 8 ) ckm1
  real    ( kind = 8 ) d0
  real    ( kind = 8 ) dd
  real    ( kind = 8 ) dk
  real    ( kind = 8 ) dk1
  real    ( kind = 8 ) dk2
  real    ( kind = 8 ) dkm1
  real    ( kind = 8 ) e0
  real    ( kind = 8 ) ek
  real    ( kind = 8 ) ek1
  real    ( kind = 8 ) ekm1
  real    ( kind = 8 ) h
  real    ( kind = 8 ) h0
  integer ( kind = 4 ) iscale
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nexp
  external             recur
  real    ( kind = 8 ) rk1
  real    ( kind = 8 ) rk2
  real    ( kind = 8 ) rs
  real    ( kind = 8 ) scale
  real    ( kind = 8 ) t(0:n)
  real    ( kind = 8 ) wt
  real    ( kind = 8 ) xnode
!
!  Check for special case.
!
!  Use Christoffel-Darboux result.
!
  if ( m == n ) then

    bk1 = 0.0D+00
    bk2 = 1.0D+00
    dk1 = 0.0D+00
    dk2 = 0.0D+00
    iscale = 0
    k = 0
    call recur ( k, h, d0, e0 )

    do k = 0, n - 1

      call recur ( k, ck, dk, ek )

      if ( 0 < k ) then
        h = - ek * h
      end if

      bb = ( ck * xnode + dk ) * bk2 + ek * bk1
      dd = ( ck * xnode + dk ) * dk2 + ek * dk1 + ck * bk2
      bk1 = bk2
      bk2 = bb
      dk1 = dk2
      dk2 = dd

      if ( bk2 /= 0.0D+00 ) then
        j = int ( log10 ( abs ( bk2 ) ) )
        if ( 2 < abs ( j ) ) then
!
!  Scale to control overflow/underflow.
!
          iscale = iscale - 2 * j
          scale = 10.0D+00**j
          bk2 = bk2 / scale
          bk1 = bk1 / scale
          dk1 = dk1 / scale
          dk2 = dk2 / scale
        end if
      end if

      if ( h /= 0.0D+00 ) then
        j = int ( log10 ( abs ( h ) ) )
        if ( 2 <= abs ( j ) ) then
          iscale = iscale + j
          h = h / 10.0D+00**j
        end if
      end if

    end do

    wt = h0 * h / dk2 / bk1

    if ( wt /= 0.0D+00 ) then
      itest = int ( log10 ( abs ( wt ) ) ) + iscale
      if ( - nexp <= itest ) then
        wt = wt * 10.0D+00**iscale
      else
        wt = 0.0D+00
      end if
    end if
    return
  end if
!
!  General case.
!
  bk2 = 0.0D+00
  bk1 = 0.0D+00
  rk2 = 0.0D+00
  rk1 = 0.0D+00
  call recur ( n, ck, dk, ek )
  call recur ( n + 1, ck1, dk1, ek1 )
  h = 1.0D+00
  iscale = 0

  do k = n, 1, -1

    if ( m <= k ) then

      rs = t(k-m) / h
!
!  Scale and adjust for possible overflow/underflow.
!
      if ( nexp < iscale ) then
        rs = 0.0D+00
      else
        rs = rs / 10.0D+00**iscale
      end if
    else
      rs = 0.0D+00
    end if

    bb = rs + ( dk + xnode * ck ) * bk1 + ek1 * bk2
    bk2 = bk1
    bk1 = bb

    call recur ( k - 1, ckm1, dkm1, ekm1 )

    if ( n /= m ) then
      h = - h * ck / ek / ckm1
    end if

    bb = bb * ckm1
    wt = bb + ( dkm1 + xnode * ckm1 ) * rk1 + ek * rk2
    rk2 = rk1
    rk1 = wt
    ck1 = ck
    dk1 = dk
    ek1 = ek
    ck = ckm1
    dk = dkm1
    ek = ekm1

    if ( bk1 /= 0.0D+00 ) then

      j = int ( log10 ( abs ( bk1 ) ) )
!
!  Scale to control overflow/underflow.
!
      if ( 2 < abs ( j ) ) then
        iscale = iscale + j
        scale = 10.0D+00**j
        bk1 = bk1 / scale
        bk2 = bk2 / scale
        rk1 = rk1 / scale
        rk2 = rk2 / scale
      end if
    end if

  end do

  wt = h0 * bb / wt

  return
end
