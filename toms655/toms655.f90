subroutine cawiq ( nt, t, mlt, nwts, ndx, key, nst, aj, bj, jdf, zemu, wts )

!*****************************************************************************80
!
!! CAWIQ computes quadrature weights for a given set of knots.
!
!  Discussion:
!
!    This routine is given a set of distinct knots, T, their multiplicities MLT,
!    the Jacobi matrix associated with the polynomials orthogonal with respect
!    to the weight function W(X), and the zero-th moment of W(X).
!
!    It computes the weights of the quadrature formula
!
!      sum ( 1 <= J <= NT ) sum ( 0 <= I <= MLT(J) - 1 ) wts(j) d^i/dx^i f(t(j))
!
!    which is to approximate
!
!      integral ( a < x < b ) f(t) w(t) dt
!
!    The routine makes various checks, as indicated below, sets up
!    various vectors and, if necessary, calls for the diagonalization
!    of the Jacobi matrix that is associated with the polynomials
!    orthogonal with respect to W(X) on the interval A, B.
!
!    Then for each knot, the weights of which are required, it calls the
!    routine CWIQD which to compute the weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input/output, integer ( kind = 4 ) NDX(NT), associates with each distinct
!    knot T(J), an integer NDX(J) which is such that the weight to the I-th
!    derivative value of F at the J-th knot, is stored in
!      WTS(abs(NDX(J))+I) for J = 1,2,...,NT, and I = 0,1,2,...,MLT(J)-1.
!    The sign of NDX includes the following information:
!    > 0, weights are wanted for this knot
!    < 0, weights not wanted for this knot but it is included in the quadrature
!    = 0. means ignore this knot completely.
!
!    Input, integer ( kind = 4 ) KEY, indicates structure of WTS and NDX.
!    KEY is an integer with absolute value between 1 and 4.
!    The sign of KEY choosed the form of WTS:
!    0 < KEY, WTS in standard form.
!    0 > KEY, J]WTS(J) required.
!    The absolute value has the following effect:
!    1, set up pointers in NDX for all knots in T array (routine CAWIQ does
!    this).  the contents of NDX are not tested on input and weights are
!    packed sequentially in WTS as indicated above.
!    2, set up pointers only for knots which have nonzero NDX on input.  All
!    knots which have a non-zero flag are allocated space in WTS.
!    3, set up pointers only for knots which have NDX > 0 on input.  Space in
!    WTS allocated only for knots with NDX > 0.
!    4, NDX assumed to be preset as pointer array on input.
!
!    Input, integer ( kind = 4 ) NST, the dimension of the Jacobi matrix.
!    NST should be between (N+1)/2 and N.  The usual choice will be (N+1)/2.
!
!    Input/output, real ( kind = 8 ) AJ(NST), BJ(NST).
!    If JDF = 0 then AJ contains the  diagonal of the Jacobi matrix and
!    BJ(1:NST-1) contains the subdiagonal.
!    If JDF = 1, AJ contains the eigenvalues of the Jacobi matrix and
!    BJ contains the squares of the elements of the first row of U, the
!    orthogonal matrix which diagonalized the Jacobi matrix as U*D*U'.
!
!    Input/output, integer ( kind = 4 ) JDF, indicates whether the Jacobi
!    matrix needs to be diagonalized.
!    0, diagonalization required;
!    1, diagonalization not required.
!
!    Input, real ( kind = 8 ) ZEMU, the zero-th moment of the weight
!    function W(X).
!
!    Output, real ( kind = 8 ) WTS(NWTS), the weights.
!
  implicit none

  integer ( kind = 4 ) nst
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( kind = 8 ) aj(nst)
  real ( kind = 8 ) bj(nst)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jdf
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) k
  integer ( kind = 4 ) key
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) mnm
  integer ( kind = 4 ) mtj
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndx(nt)
  real ( kind = 8 ) p
  real ( kind = 8 ) prec
  real ( kind = 8 ), allocatable :: r(:)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) tmp
  real ( kind = 8 ) wts(nwts)
  real ( kind = 8 ), allocatable :: xk(:)
  real ( kind = 8 ), allocatable :: z(:)
  real ( kind = 8 ) zemu

  prec = epsilon ( prec )

  if ( nt < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CAWIQ - Fatal error!'
    write ( *, '(a)' ) '  NT < 1.'
    stop
  end if
!
!  Check for indistinct knots.
!
  if ( 1 < nt ) then

    k = nt - 1

    do i = 1, k
      tmp = t(i)
      l = i + 1
      do j = l, nt
        if ( abs ( tmp - t(j) ) <= prec ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CAWIQ - Fatal error!'
          write ( *, '(a)' ) '  Knots too close.'
          stop
        end if
      end do

    end do

  end if
!
!  Check multiplicities,
!  Set up various useful parameters and
!  set up or check pointers to WTS array.
!
  l = abs ( key )

  if ( l < 1 .or. 4 < l ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CAWIQ - Fatal error!'
    write ( *, '(a)' ) '  Magnitude of KEY not between 1 and 4.'
    stop
  end if

  k = 1

  if ( l == 1 ) then

    do i = 1, nt
      ndx(i) = k
      if ( mlt(i) < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CAWIQ - Fatal error!'
        write ( *, '(a)' ) '  MLT(I) < 1.'
        stop
      end if
      k = k + mlt(i)
    end do

    n = k - 1

  else if ( l == 2 .or. l == 3 ) then

    n = 0

    do i = 1, nt

      if ( ndx(i) == 0 ) then
        cycle
      end if

      if ( mlt(i) < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CAWIQ - Fatal error!'
        write ( *, '(a)' ) '  MLT(I) < 1.'
        stop
      end if

      n = n + mlt(i)

      if ( ndx(i) < 0 .and. l == 3 ) then
        cycle
      end if

      ndx(i) = sign ( k, ndx(i) )
      k = k + mlt(i)

    end do

    if ( nwts + 1 < k ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CAWIQ - Fatal error!'
      write ( *, '(a)' ) '  NWTS + 1 < K.'
      stop
    end if

  else if ( l == 4 ) then

    do i = 1, nt

      ip = abs ( ndx(i) )

      if ( ip == 0 ) then
        cycle
      end if

      if ( nwts < ip + mlt(i) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CAWIQ - Fatal error!'
        write ( *, '(a)' ) '  NWTS < IPM.'
        stop
      end if

      if ( i == nt ) then
        exit
      end if

      l = i + 1
      do j = l, nt
        jp = abs ( ndx(j) )
        if ( jp /= 0 ) then
          if ( jp <= ip + mlt(i) .and. ip <= jp + mlt(j) ) then
            exit
          end if
        end if
      end do

    end do

  end if
!
!  Test some parameters.
!
  if ( nst < ( n + 1 ) / 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CAWIQ - Fatal error!'
    write ( *, '(a)' ) '  NST < ( N + 1 ) / 2.'
    stop
  end if

  if ( zemu <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CAWIQ - Fatal error!'
    write ( *, '(a)' ) '  ZEMU <= 0.'
    stop
  end if
!
!  Treat a quadrature formula with 1 simple knot first.
!
  if ( n <= 1 ) then

    do i = 1, nt
      if ( 0 < ndx(i) ) then
        wts( abs ( ndx(i) ) ) = zemu
        return
      end if
    end do

  end if
!
!  Carry out diagonalization if not already done.
!
  if ( jdf == 0 ) then
!
!  Set unit vector in work field to get back first row of Q.
!
    allocate ( z(1:nst) )

    do i = 1, nst
      z(i) = 0.0D+00
    end do
    z(1) = 1.0D+00
!
!  Diagonalize the Jacobi matrix.
!
    call imtqlx ( nst, aj, bj, z )
!
!  Signal Jacobi matrix now diagonalized successfully.
!
    jdf = 1
!
!  Save squares of first row of U in subdiagonal array.
!
    do i = 1, nst
      bj(i) = z(i) * z(i)
    end do

    deallocate ( z )

  end if
!
!  Find all the weights for each knot flagged.
!
  do i = 1, nt

    if ( ndx(i) <= 0 ) then
      cycle
    end if

    m = mlt(i)
    mnm = max ( n - m, 1 )
    l = min ( m, n - m + 1 )
!
!  Set up K-hat matrix for CWIQD with knots according to their multiplicities.
!
    allocate ( xk(1:mnm) )

    k = 1
    do j = 1, nt
      if ( ndx(j) /= 0 ) then
        if ( j /= i ) then
          do jj = 1, mlt(j)
            xk(k) = t(j)
            k = k + 1
          end do
        end if
      end if
    end do
!
!  Set up the right principal vector.
!
    allocate ( r(1:l) )

    r(1) = 1.0D+00 / zemu
    r(2:l) = 0.0D+00
!
!  Pick up pointer for the location of the weights to be output.
!
    k = ndx(i)
!
!  Find all the weights for this knot.
!
    call cwiqd ( m, mnm, l, t(i), xk, nst, aj, bj, r, wts(k) )

    deallocate ( r )
    deallocate ( xk )

    if ( key < 0 ) then
      cycle
    end if
!
!  Divide by factorials for weights in standard form.
!
    tmp = 1.0D+00
    do j = 2, m - 1
      p = j
      tmp = tmp * p
      wts(k+j) = wts(k+j) / tmp
    end do

  end do

  return
end
subroutine cdgqf ( nt, kind, alpha, beta, t, wts )

!*****************************************************************************80
!
!! CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with a classical weight function with default values for A and B,
!    and only simple knots.
!
!    There are no moments checks and no printing is done.
!
!    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) aj(nt)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) kind
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) zemu

  call parchk ( kind, 2 * nt, alpha, beta )
!
!  Get the Jacobi matrix and zero-th moment.
!
  call class_matrix ( kind, nt, alpha, beta, aj, bj, zemu )
!
!  Compute the knots and weights.
!
  call sgqf ( nt, aj, bj, zemu, t, wts )

  return
end
subroutine cegqf ( nt, kind, alpha, beta, a, b, f, qfsum )

!*****************************************************************************80
!
!! CEGQF computes a quadrature formula and applies it to a function.
!
!  Discussion:
!
!    The user chooses the quadrature formula to be used, as well as the
!    interval (A,B) in which it is applied.
!
!    Note that the knots and weights of the quadrature formula are not
!    returned to the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
!    Input, real ( kind = 8 ), external F, the name of a routine which
!    evaluates the function and some of its derivatives.  The routine
!    must have the form
!      function f ( x, i )
!      real ( kind = 8 ) f
!      integer ( kind = 4 ) i
!      real ( kind = 8 ) x
!    and return in F the value of the I-th derivative of the function
!    at X.  The value I will always be 0.  The value X will always be a knot.
!
!    Output, real ( kind = 8 ) QFSUM, the value of the quadrature formula
!    applied to F.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) kind
  real ( kind = 8 ) qfsum
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)

  call cgqf ( nt, kind, alpha, beta, a, b, t, wts )
!
!  Evaluate the quadrature sum.
!
  call eiqfs ( nt, t, wts, f, qfsum )

  return
end
subroutine cegqfs ( nt, kind, alpha, beta, f, qfsum )

!*****************************************************************************80
!
!! CEGQFS estimates an integral using a standard quadrature formula.
!
!  Discussion:
!
!    The user chooses one of the standard quadrature rules
!    with the default values of A and B.  This routine determines
!    the corresponding weights and evaluates the quadrature formula
!    on a given function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ), external F, the name of a routine which
!    evaluates the function and some of its derivatives.  The routine
!    must have the form
!      function f ( x, i )
!      real ( kind = 8 ) f
!      integer ( kind = 4 ) i
!      real ( kind = 8 ) x
!    and return in F the value of the I-th derivative of the function
!    at X.  The value  I will always be 0.  The value X will always be a knot.
!
!    Output, real ( kind = 8 ) QFSUM, the value of the quadrature formula
!    applied to F.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lu
  real ( kind = 8 ) qfsum
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
!
!  Assign workspace for knots and weights.
!
  lu = 0

  call cgqfs ( nt, kind, alpha, beta, lu, t, wts )
!
!  Evaluate the quadrature sum.
!
  call eiqfs ( nt, t, wts, f, qfsum )

  return
end
subroutine ceiqf ( nt, t, mlt, kind, alpha, beta, a, b, f, qfsum )

!*****************************************************************************80
!
!! CEIQF constructs and applies a quadrature formula based on user knots.
!
!  Discussion:
!
!    The knots may have multiplicity.  The quadrature interval is over
!    any valid A, B.  A classical weight function is selected by the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
!    Input, real ( kind = 8 ), external F, the name of a routine which
!    evaluates the function and some of its derivatives.  The routine
!    must have the form
!      function f ( x, i )
!      real ( kind = 8 ) f
!      integer ( kind = 4 ) i
!      real ( kind = 8 ) x
!    and return in F the value of the I-th derivative of the function
!    at X.  The highest value of I will be the maximum value in MLT minus
!    one.  The value X will always be a knot.
!
!    Output, real ( kind = 8 ) QFSUM, the value of the quadrature formula
!    applied to F.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) key
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lu
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ), allocatable :: ndx(:)
  integer ( kind = 4 ) nwts
  real ( kind = 8 ) qfsum
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ), allocatable :: wts(:)

  lu = 0
  nwts = sum ( mlt(1:nt) )

  key = 1
  allocate ( ndx(1:nt) )
  allocate ( wts(1:nwts) )

  call ciqf ( nt, t, mlt, nwts, ndx, key, kind, alpha, beta, a, b, lu, wts )

  call eiqf ( nt, t, mlt, wts, nwts, ndx, key, f, qfsum )

  deallocate ( ndx )
  deallocate ( wts )

  return
end
subroutine ceiqfs ( nt, t, mlt, kind, alpha, beta, f, qfsum )

!*****************************************************************************80
!
!! CEIQFS computes and applies a quadrature formula based on user knots.
!
!  Discussion:
!
!    The knots may have multiplicity.  The quadrature interval is over
!    the standard interval A, B for the classical weight function selected
!    by the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ), external F, the name of a routine which
!    evaluates the function and some of its derivatives.  The routine
!    must have the form
!      function f ( x, i )
!      real ( kind = 8 ) f
!      integer ( kind = 4 ) i
!      real ( kind = 8 ) x
!    and return in F the value of the I-th derivative of the function
!    at X.  The highest value of I will be the maximum value in MLT minus
!    one.  The value X will always be a knot.
!
!    Output, real ( kind = 8 ) QFSUM, the value of the quadrature formula
!    applied to F.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) key
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lu
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) n
  integer ( kind = 4 ), allocatable :: ndx(:)
  real ( kind = 8 ) qfsum
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ), allocatable :: wts(:)

  lu = 0
  n = sum ( mlt(1:nt) )

  allocate ( ndx(1:nt) )
  key = 1
  allocate ( wts(1:n) )

  call ciqfs ( nt, t, mlt, n, ndx, key, kind, alpha, beta, lu, wts )

  call eiqf ( nt, t, mlt, wts, n, ndx, key, f, qfsum )

  deallocate ( ndx )
  deallocate ( wts )

  return
end
subroutine cgqf ( nt, kind, alpha, beta, a, b, t, wts )

!*****************************************************************************80
!
!! CGQF computes knots and weights of a Gauss quadrature formula.
!
!  Discussion:
!
!    The user may specify the interval (A,B).
!
!    Only simple knots are produced.
!
!    Use routine EIQFS to evaluate this quadrature formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints, or
!    other parameters.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ), allocatable :: ndx(:)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
!
!  Compute the Gauss quadrature formula for default values of A and B.
!
  call cdgqf ( nt, kind, alpha, beta, t, wts )
!
!  Prepare to scale the quadrature formula to other weight function with
!  valid A and B.
!
  allocate ( mlt(1:nt) )

  mlt(1:nt) = 1

  allocate ( ndx(1:nt) )

  do i = 1, nt
    ndx(i) = i
  end do

  call scqf ( nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b )

  deallocate ( mlt )
  deallocate ( ndx )

  return
end
subroutine cgqfs ( nt, kind, alpha, beta, lo, t, wts )

!*****************************************************************************80
!
!! CGQFS computes knots and weights of a Gauss quadrature formula.
!
!  Discussion:
!
!    This routine computes the knots and weights of a Gauss quadrature
!    formula with:
!
!    * a classical weight function with default values for A and B;
!    * only simple knots
!    * optionally print knots and weights and a check of the moments
!
!    Use routine EIQFS to evaluate a quadrature formula computed by
!    this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, integer ( kind = 4 ) LO, selects the action.
!    > 0, compute and print knots and weights.  Print moments check.
!    = 0, compute knots and weights.
!    < 0, compute and print knots and weights.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) key
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) m
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ) mmex
  integer ( kind = 4 ) mop
  integer ( kind = 4 ), allocatable :: ndx(:)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ) wts(nt)
!
!  Check there is enough workfield and assign workfield
!
  key = 1
  mop = 2 * nt
  m = mop + 1
  mmex = max ( m + 2, 1 )
!
!  Compute the Gauss quadrature formula for default values of A and B.
!
  call cdgqf ( nt, kind, alpha, beta, t, wts )
!
!  Exit if no print required.
!
  if ( lo /= 0 ) then

    allocate ( mlt(1:nt) )

    mlt(1:nt) = 1

    allocate ( ndx(1:nt) )

    do i = 1, nt
      ndx(i) = i
    end do

    allocate ( w(1:mmex) )

    call chkqfs ( t, wts, mlt, nt, nt, ndx, key, w, mop, mmex, kind, alpha, &
      beta, lo )

    deallocate ( mlt )
    deallocate ( ndx )
    deallocate ( w )

  end if

  return
end
subroutine chkqf ( t, wts, mlt, nt, nwts, ndx, key, mop, mex, kind, alpha, &
  beta, lo, a, b )

!*****************************************************************************80
!
!! CHKQF computes and prints the moments of a quadrature formula.
!
!  Discussion:
!
!    The quadrature formula is based on a clasical weight function with
!    any valid A, B.
!
!    No check can be made for non-classical weight functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(NT), the knots.
!
!    Input, real ( kind = 8 ) WTS(NWTS), the weights.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input, integer ( kind = 4 ) NDX(NT), used to index the array WTS.
!    If KEY = 1, then NDX need not be preset.  For more details see the
!    comments in CAWIQ.
!
!    Input, integer ( kind = 4 ) KEY, indicates the structure of the WTS
!    array.  It will normally be set to 1.  This will cause the weights to be
!    packed sequentially in array WTS.  For more details see the comments
!    in CAWIQ.
!
!    Input, integer ( kind = 4 ) MOP, the expected order of precision of the
!    quadrature formula.
!
!    Input, integer ( kind = 4 ) MEX, the number of moments required to be
!    tested.  Set MEX = 1 and LO < 0 for no moments check.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, integer ( kind = 4 ) LO, selects the action to carry out.
!     > 0, print weights and moment tests.
!     = 0, print nothing. compute moment test.
!     < 0, print weights only. don't compute moment tests.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
  implicit none

  integer ( kind = 4 ) mex
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) izero
  integer ( kind = 4 ) key
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mop
  integer ( kind = 4 ) neg
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) ndx(nt)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ), allocatable :: t2(:)
  real ( kind = 8 ) tmp
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ) wts(nwts)

  allocate ( w(1:mex) )

  call parchk ( kind, mex, alpha, beta )

  if ( lo /= 0 ) then

    izero = 0

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Interpolatory quadrature formula'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  Type  Interval       Weight function               Name'
    write ( *, '(a)' ) ' '
    if ( kind == 1 ) then
      write ( *, '(a)' ) &
        '    1  (a,b)              1.0                    Legendre'
    else if ( kind == 2 ) then
      write ( *, '(a)' ) &
        '    2  (a,b)      ((b-x)*(x-a))^(-0.5)          Chebyshev Type 1'
    else if ( kind == 3 ) then
      write ( *, '(a)' ) &
        '    3  (a,b)      ((b-x)*(x-a))^alpha           Gegenbauer'
    else if ( kind == 4 ) then
      write ( *, '(a)' ) &
        '    4  (a,b)    (b-x)^alpha*(x-a)^beta          Jacobi'
    else if ( kind == 5 ) then
      write ( *, '(a)' ) &
        '    5  (a,+oo)  (x-a)^alpha*exp(-b*(x-a))      Gen Laguerre'
    else if ( kind == 6 ) then
      write ( *, '(a)' ) &
        '    6  (-oo,+oo) |x-a|^alpha*exp(-b*(x-a)^2)    Gen Hermite'
    else if ( kind == 7 ) then
      write ( *, '(a)' ) &
        '    7  (a,b)      |x-(a+b)/2.0|^alpha        Exponential'
    else if ( kind == 8 ) then
      write ( *, '(a)' ) &
        '    8  (a,+oo)   (x-a)^alpha*(x+b)^beta         Rational'
    else if ( kind == 9 ) then
      write ( *, '(a)' ) &
        '    9  (a,b)      ((b-x)*(x-a))^(+0.5)          Chebyshev Type 2'
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,f12.5)' ) '     Parameters   A          ', a
    write ( *, '(a,f12.5)' ) '                  B          ', b
    if ( 3 <= kind .and. kind <= 8 ) then
      write ( *, '(a,f12.5)' ) '                  alpha      ', alpha
    end if

    if ( kind == 4 .or. kind == 8 ) then
      write ( *, '(a,f12.5)' ) '                  beta       ', beta
    end if

    call chkqfs ( t, wts, mlt, nt, nwts, ndx, key, w, mop, mex, izero, &
      alpha, beta, - abs ( lo ) )

  end if

  if ( 0 <= lo ) then
!
!  Compute the moments in W.
!
    call scmm ( mex, kind, alpha, beta, a, b, w )

    if ( kind == 1 .or. kind == 2 .or. kind == 3 .or. kind == 4 &
      .or. kind == 7 .or. kind == 9)  then
      tmp = ( b + a ) / 2.0D+00
    else if ( kind == 5 .or. kind == 6 .or. kind == 8 ) then
      tmp = a
    end if

    allocate ( t2(1:nt))

    t2(1:nt) = t(1:nt) - tmp

    neg = -1
!
!  Check moments.
!
    call chkqfs ( t2, wts, mlt, nt, nwts, ndx, key, w, mop, mex, neg, &
      alpha, beta, lo  )

    deallocate ( t2 )

  end if

  deallocate ( w )

  return
end
subroutine chkqfs ( t, wts, mlt, nt, nwts, ndx, key, w, mop, mex, kind, &
  alpha, beta, lo )

!*****************************************************************************80
!
!! CHKQFS checks the polynomial accuracy of a quadrature formula.
!
!  Discussion:
!
!    This routine will optionally print weights, and results of a moments test.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(NT), the knots.
!
!    Input, real ( kind = 8 ) WTS(NWTS), the weights.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input, integer ( kind = 4 ) NDX(NT), used to index the array WTS.
!    If KEY = 1, then NDX need not be preset.  For more details see the
!    comments in CAWIQ.
!
!    Input, integer ( kind = 4 ) KEY, indicates the structure of the WTS
!    array.  It will normally be set to 1.  This will cause the weights to be
!    packed sequentially in array WTS.  For more details see the comments
!    in CAWIQ.
!
!    Input/output, real ( kind = 8 ) W(MEX), the moments array.
!    This is input only if KIND = 0.
!
!    Input, integer ( kind = 4 ) MOP, the expected order of precision of the
!    quadrature formula.
!
!    Input, integer ( kind = 4 ) MEX, the number of moments to be tested.
!    MEX must be at least 1.  Set MEX = 1 and LO < 0 for no moment check.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    0, unknown weight function (the user must set the first MEX moments in
!       array W in this case.)
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, integer ( kind = 4 ) LO, selects the action to carry out.
!     > 0, print weights and moment tests.
!     = 0, print nothing. compute moment test.
!     < 0, print weights only. don't compute moment tests.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) E(MEX), ER(MEX), the absolute and relative
!    errors of the quadrature formula applied to (X-DEL)^n.
!
!    Local, real ( kind = 8 ) QM(MEX), the values of the quadrature formula
!    applied to (X-DEL)^N.
!
  implicit none

  integer ( kind = 4 ) mex
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) e(mex)
  real ( kind = 8 ) ek
  real ( kind = 8 ) emn
  real ( kind = 8 ) emx
  real ( kind = 8 ) erest
  real ( kind = 8 ) ern
  real ( kind = 8 ) erx
  real ( kind = 8 ) er(mex)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jl
  integer ( kind = 4 ) k
  integer ( kind = 4 ) key
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) kindp
  integer ( kind = 4 ) kjl
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) mop
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) ndx(nt)
  real ( kind = 8 ) px
  real ( kind = 8 ) tmp
  real ( kind = 8 ) tmpx
  real ( kind = 8 ) prec
  real ( kind = 8 ) qm(mex)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) w(mex)
  real ( kind = 8 ) wts(nwts)
!
!  KIND may be set to -1 to allow printing of moments only.
!
!  This feature is only used internally, by CHKQF.
!
  kindp = max ( 0, kind )

  if ( lo /= 0 .and. kind /= -1 ) then

    if ( kindp /= 0 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Interpolatory quadrature formula'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) &
        '  Type  Interval       Weight function               Name'
      write ( *, '(a)' ) ' '
      if ( kindp == 1 ) then
        write ( *, '(a)' ) &
          '    1  (-1,+1)            1.0                    Legendre'
      else if ( kindp == 2 ) then
        write ( *, '(a)' ) &
          '    2  (-1,+1)    ((b-x)*(x-a))^(-0.5)          Chebyshev Type 1'
      else if ( kindp == 3 ) then
        write ( *, '(a)' ) &
          '    3  (-1,+1)    ((b-x)*(x-a))^alpha           Gegenbauer'
      else if ( kindp == 4 ) then
        write ( *, '(a)' ) &
          '    4  (-1,+1)  (b-x)^alpha*(x-a)^beta          Jacobi'
      else if ( kindp == 5 ) then
        write ( *, '(a)' ) &
          '    5  (a,+oo)   (x-a)^alpha*exp(-b*(x-a))     Gen Laguerre'
      else if ( kindp == 6 ) then
        write ( *, '(a)' ) &
          '    6  (-oo,+oo) |x-a|^alpha*exp(-b*(x-a)^2)  Gen Hermite'
      else if ( kindp == 7 ) then
        write ( *, '(a)' ) &
          '    7  (-1,+1)    |x-(a+b)/2.0|^alpha        Exponential'
      else if ( kindp == 8 ) then
        write ( *, '(a)' ) &
          '    8  (0,+oo)    (x-a)^alpha*(x+b)^beta         Rational'
      else if ( kindp == 9 ) then
        write ( *, '(a)' ) &
          '    9  (-1,+1)    ((b-x)*(x-a))^(+0.5)          Chebyshev Type 2'
      end if

      if ( 3 <= kindp .and. kindp <= 8 ) then
        write ( *, '(a,f12.5)' ) '                  alpha      ', alpha
      end if

      if ( kindp == 4 .or. kindp == 8 ) then
        write ( *, '(a,f12.5)' ) '                  beta       ', beta
      end if

    end if

    if ( kind /= -1 ) then
      prec = epsilon ( prec )
      write ( *, '(a)' ) ' '
      write ( *, '(a,d13.1)' ) '  Machine precision = ', prec
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '           Knots               Mult                Weights'
    write ( *, '(a)' ) ' '

    do i = 1, nt
      k = abs ( ndx(i) )
      if ( k /= 0 ) then
        write ( *, '(i4,d26.17,i4,d26.17)' ) i, t(i), mlt(i), wts(k)
        do j = k + 1, k + mlt(i) - 1
          write ( *, '(34x,d26.17)' ) wts(j)
        end do
      end if
    end do

  end if

  if ( 0 <= lo ) then
!
!  Compute the moments in W.
!
    if ( kindp /= 0 ) then

      call wm ( mex, kindp, alpha, beta, w )

    end if

    qm(1:mex) = 0.0D+00
    erest = 0.0D+00

    do k = 1, nt

      tmp = 1.0D+00
      l = abs ( ndx(k) )
      if ( l == 0 ) then
        cycle
      end if

      erest = erest + abs ( wts(l) )
      do j = 1, mex
        qm(j) = qm(j) + tmp * wts(l)
        tmpx = tmp
        px = 1.0D+00
        do jl = 2, min ( mlt(k), mex - j + 1 )
          kjl = j + jl - 1
          tmpx = tmpx * ( kjl - 1 )
          qm(kjl) = qm(kjl) + tmpx * wts(l+jl-1) / px
          if ( key <= 0 ) then
            px = px * jl
          end if
        end do
        tmp = tmp * t(k)
      end do

    end do

    e(1:mex) = w(1:mex) - qm(1:mex)
    er(1:mex) = e(1:mex) / ( abs ( w(1:mex) ) + 1.0D+00 )
    erest = erest / ( abs ( w(1) ) + 1.0D+00 )

  end if

  if ( 0 < lo ) then

    m = mop + 1
    mx = min ( mop, mex )

    emx = maxval ( abs ( e(1:mx) ) )
    emn = minval ( abs ( e(1:mx) ) )
    erx = maxval ( abs ( er(1:mx) ) )
    ern = minval ( abs ( er(1:mx) ) )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Comparison of moments'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Order of precision ', mop
    write ( *, '(a)' ) '  Errors :    Absolute    Relative'
    write ( *, '(a)' ) '  ---------+-------------------------'
    write ( *, '(a,2d12.3)' ) '  Minimum :', emn, ern
    write ( *, '(a,2d12.3)' ) '  Maximum :', emx, erx
    write ( *, '(a)' ) ' '
    write ( *, '(a,d13.3)' ) '  Weights ratio       ', erest

    if ( m <= mex ) then

      ek = e(m)
      do j = 1, mop
        ek = ek / real ( j, kind = 8 )
      end do

      write ( *, '(a,i2,a,d13.3)' ) '  Error in ', mop, 'th power ', e(m)
      write ( *, '(a,d13.3)' ) '  Error constant      ', ek

    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Moments:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '            True             from QF            Error      Relative'
    write ( *, '(a)' ) ' '
    do j = 1, mx
      write ( *, '(i4,2d19.10,2d12.3)' ) j, w(j), qm(j), e(j), er(j)
    end do
    write ( *, '(a)' ) ' '
    do j = m, mex
      write ( *, '(i4,2d19.10,2d12.3)' ) j, w(j), qm(j), e(j), er(j)
    end do

  end if

  return
end
subroutine ciqf ( nt, t, mlt, nwts, ndx, key, kind, alpha, beta, a, b, lo, &
  wts )

!*****************************************************************************80
!
!! CIQF computes weights for a classical weight function and any interval.
!
!  Discussion:
!
!    This routine compute somes or all the weights of a quadrature formula
!    for a classical weight function with any valid A, B and a given set of
!    knots and multiplicities.
!
!    The weights may be packed into the output array WTS according to a
!    user-defined pattern or sequentially.
!
!    The routine will also optionally print knots and weights and a check
!    of the moments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input/output, integer ( kind = 4 ) NDX(NT), used to index the output
!    array WTS.  If KEY = 1, then NDX need not be preset.  For more
!    details see the comments in CAWIQ.
!
!    Input, integer ( kind = 4 ) KEY, indicates the structure of the WTS
!    array.  It will normally be set to 1.  This will cause the weights to be
!    packed sequentially in array WTS.  For more details see the comments
!    in CAWIQ.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
!    Input, integer ( kind = 4 ) LO, selects the actions to perform.
!     > 0, compute and print weights.  Print moments check.
!     = 0, compute weights.
!     < 0, compute and print weights.
!
!    Output, real ( kind = 8 ) WTS(NWTS), the weights.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) key
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) lu
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mex
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) mop
  integer ( kind = 4 ) ndx(nt)
  real ( kind = 8 ) st(nt)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nwts)

  m = 1
  l = abs ( key )

  do j = 1, nt
    if ( l == 1 .or. abs ( ndx(j) ) /= 0 ) then
      m = m + mlt(j)
    end if
  end do

  if ( nwts + 1 < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIQF - Fatal error!'
    write ( *, '(a)' ) '  NWTS + 1 < M.'
    stop
  end if

  mex = 2 + m
!
!  Scale the knots to default A, B.
!
  call sct ( nt, t, kind, a, b, st )

  lu = 0

  call ciqfs ( nt, st, mlt, nwts, ndx, key, kind, alpha, beta, lu, wts )
!
!  Don't scale user's knots - only scale weights.
!
  call scqf ( nt, st, mlt, wts, nwts, ndx, wts, st, kind, alpha, &
    beta, a, b )

  if ( lo /= 0 ) then

    mop = m - 1

    call chkqf ( t, wts, mlt, nt, nwts, ndx, key, mop, mex, kind, &
      alpha, beta, lo, a, b )

  end if

  return
end
subroutine ciqfs ( nt, t, mlt, nwts, ndx, key, kind, alpha, beta, lo, wts )

!*****************************************************************************80
!
!! CIQFS computes some weights of a quadrature formula in the default interval.
!
!  Discussion:
!
!    This routine computes some or all the weights of a quadrature formula
!    for a classical weight function with default values of A and B,
!    and a given set of knots and multiplicities.
!
!    The weights may be packed into the output array WTS according to a
!    user-defined pattern or sequentially.
!
!    The routine will also optionally print knots and weights and a check of
!    the moments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input/output, integer ( kind = 4 ) NDX(NT),  used to index the output
!    array WTS.  If KEY = 1, then NDX need not be preset.  For more
!    details see the comments in CAWIQ.
!
!    Input, integer ( kind = 4 ) KEY, indicates the structure of the WTS
!    array.  It will normally be set to 1.  For more details see
!    the comments in CAWIQ.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, integer ( kind = 4 ) LO, selects the actions to perform.
!     > 0, compute and print weights.  Print moments check.
!     = 0, compute weights.
!     < 0, compute and print weights.
!
!    Output, real ( kind = 8 ) WTS(NWTS), the weights.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( kind = 8 ), allocatable :: aj(:)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ), allocatable :: bj(:)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jdf
  integer ( kind = 4 ) key
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mex
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) mmex
  integer ( kind = 4 ) mop
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndx(nt)
  integer ( kind = 4 ) nst
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ) wts(nwts)
  real ( kind = 8 ) zemu

  jdf = 0
  n = 0
  l = abs ( key )

  do j = 1, nt

    if ( l == 1 .or. abs ( ndx(j) ) /= 0 ) then
      n = n + mlt(j)
    end if

  end do
!
!  N knots when counted according to multiplicity.
!
  if ( nwts < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIQFS - Fatal error!'
    write ( *, '(a)' ) '  NWTS < N.'
    stop
  end if

  mex = n + 3
  nst = ( n + 1 ) / 2
!
!  Get the Jacobi matrix.
!
  allocate ( aj(1:nst) )
  allocate ( bj(1:nst) )

  call class_matrix ( kind, nst, alpha, beta, aj, bj, zemu )
!
!  Call weights routine.
!
  call cawiq ( nt, t, mlt, n, ndx, key, nst, aj, bj, jdf, zemu, wts )

  deallocate ( aj )
  deallocate ( bj )
!
!  Return if no printing or checking required.
!
  if ( lo /= 0 ) then

    mop = n

    allocate ( w(1:mex) )

    call chkqfs ( t, wts, mlt, nt, n, ndx, key, w, mop, mex, kind, &
      alpha, beta, lo )

    deallocate ( w )

  end if

  return
end
subroutine class_matrix ( kind, m, alpha, beta, aj, bj, zemu )

!*****************************************************************************80
!
!! CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
!
!  Discussion:
!
!    This routine computes the diagonal AJ and sub-diagonal BJ
!    elements of the order M tridiagonal symmetric Jacobi matrix
!    associated with the polynomials orthogonal with respect to
!    the weight function specified by KIND.
!
!    For weight functions 1-7, M elements are defined in BJ even
!    though only M-1 are needed.  For weight function 8, BJ(M) is
!    set to zero.
!
!    The zero-th moment of the weight function is returned in ZEMU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, integer ( kind = 4 ) M, the order of the Jacobi matrix.
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) AJ(M), BJ(M), the diagonal and subdiagonal
!    of the Jacobi matrix.
!
!    Output, real ( kind = 8 ) ZEMU, the zero-th moment.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a2b2
  real ( kind = 8 ) ab
  real ( kind = 8 ) aba
  real ( kind = 8 ) abi
  real ( kind = 8 ) abj
  real ( kind = 8 ) abti
  real ( kind = 8 ) aj(m)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) apone
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp2
  real ( kind = 8 ) zemu

  temp = epsilon ( temp )

  call parchk ( kind, 2 * m - 1, alpha, beta )

  temp2 = r8_gamma ( 0.5D+00 )

  if ( 500.0D+00 * temp < abs ( temp2 * temp2 - pi ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLASS_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  Gamma function does not match machine parameters.'
    stop
  end if

  if ( kind == 1 ) then

    ab = 0.0D+00

    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) = sqrt ( bj(1:m) )

  else if ( kind == 2 ) then

    zemu = pi

    aj(1:m) = 0.0D+00

    bj(1) =  sqrt ( 0.5D+00 )
    bj(2:m) = 0.5D+00

  else if ( kind == 3 ) then

    ab = alpha * 2.0D+00
    zemu = 2.0D+00**( ab + 1.0D+00 ) * ( r8_gamma ( alpha + 1.0D+00 ) )**2 &
      / r8_gamma ( ab + 2.0D+00 )

    aj(1:m) = 0.0D+00
    bj(1) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
    do i = 2, m
      bj(i) = i * ( i + ab ) &
        / ( 4.0D+00 * ( i + alpha ) * ( i + alpha ) - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 4 ) then

    ab = alpha + beta
    abi = 2.0D+00 + ab
    zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 ) &
      * r8_gamma ( beta + 1.0D+00 ) / r8_gamma ( abi )
    aj(1) = ( beta - alpha ) / abi
    bj(1) = 4.0D+00 * ( 1.0 + alpha ) * ( 1.0D+00 + beta ) &
      / ( ( abi + 1.0D+00 ) * abi * abi )
    a2b2 = beta * beta - alpha * alpha

    do i = 2, m
      abi = 2.0D+00 * i + ab
      aj(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
      abi = abi * abi
      bj(i) = 4.0D+00 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) &
        / ( ( abi - 1.0D+00 ) * abi )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 5 ) then

    zemu = r8_gamma ( alpha + 1.0D+00 )

    do i = 1, m
      aj(i) = 2.0D+00 * i - 1.0D+00 + alpha
      bj(i) = i * ( i + alpha )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 6 ) then

    zemu = r8_gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      bj(i) = ( i + alpha * mod ( i, 2 ) ) / 2.0D+00
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 7 ) then

    ab = alpha
    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 8 ) then

    ab = alpha + beta
    zemu = r8_gamma ( alpha + 1.0D+00 ) * r8_gamma ( - ( ab + 1.0D+00 ) ) &
      / r8_gamma ( - beta )
    apone = alpha + 1.0D+00
    aba = ab * apone
    aj(1) = - apone / ( ab + 2.0D+00 )
    bj(1) = - aj(1) * ( beta + 1.0D+00 ) / ( ab + 2.0D+00 ) / ( ab + 3.0D+00 )
    do i = 2, m
      abti = ab + 2.0D+00 * i
      aj(i) = aba + 2.0D+00 * ( ab + i ) * ( i - 1 )
      aj(i) = - aj(i) / abti / ( abti - 2.0D+00 )
    end do

    do i = 2, m - 1
      abti = ab + 2.0D+00 * i
      bj(i) = i * ( alpha + i ) / ( abti - 1.0D+00 ) * ( beta + i ) &
        / ( abti * abti ) * ( ab + i ) / ( abti + 1.0D+00 )
    end do

    bj(m) = 0.0D+00
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 9 ) then

    zemu = pi / 2.0D+00
    aj(1:m) = 0.0D+00
    bj(1:m) = 0.5D+00

  end if

  return
end
subroutine cliqf ( nt, t, kind, alpha, beta, a, b, lo, wts )

!*****************************************************************************80
!
!! CLIQF computes a classical quadrature formula, with optional printing.
!
!  Discussion:
!
!    This routine computes all the weights of an interpolatory
!    quadrature formula with
!    1. only simple knots and
!    2. a classical weight function with any valid A and B, and
!    3. optionally prints the knots and weights and a check of the moments.
!
!    To evaluate this quadrature formula for a given function F,
!    call routine EIQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
!    Input, integer ( kind = 4 ) LO, indicates what is to be done.
!    > 0, compute and print weights and moments check.
!    = 0, compute weights.
!    < 0, compute and print weights.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) key
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lo
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ), allocatable :: ndx(:)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)

  key = 1

  allocate ( mlt(1:nt) )
  mlt(1:nt) = 1

  allocate ( ndx(1:nt) )

  call ciqf ( nt, t, mlt, nt, ndx, key, kind, alpha, beta, a, b, lo, wts )

  deallocate ( mlt )
  deallocate ( ndx )

  return
end
subroutine cliqfs ( nt, t, kind, alpha, beta, lo, wts )

!*****************************************************************************80
!
!! CLIQFS computes the weights of a quadrature formula in the default interval.
!
!  Discussion:
!
!    This routine computes the weights of an interpolatory quadrature formula
!    with a classical weight function, in the default interval A, B,
!    using only simple knots.
!
!    It can optionally print knots and weights and a check of the moments.
!
!    To evaluate a quadrature computed by CLIQFS, call EIQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, integer ( kind = 4 ) LO, chooses the printing option.
!     > 0, compute weights, print them, print the moment check results.
!     0, compute weights.
!     < 0, compute weights and print them.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) key
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) lo
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ), allocatable :: ndx(:)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)

  key = 1

  allocate ( mlt(1:nt) )
  mlt(1:nt) = 1

  allocate ( ndx(1:nt) )

  call ciqfs ( nt, t, mlt, nt, ndx, key, kind, alpha, beta, lo, wts )

  deallocate ( mlt )
  deallocate ( ndx )

  return
end
subroutine cwiqd ( m, nm, l, v, xk, nstar, phi, a, r, d )

!*****************************************************************************80
!
!! CWIQD computes all the weights for a given knot.
!
!  Discussion:
!
!    The variable names correspond to the 1982 reference, and explanations of
!    some of the terminology may be found there.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Jaroslav Kautsky, Sylvan Elhay,
!    Calculation of the Weights of Interpolatory Quadratures,
!    Numerische Mathematik,
!    Volume 40, 1982, pages 407-422.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the multiplicity of the knot in question.
!
!    Input, integer ( kind = 4 ) NM, is equal to max ( N - M, 1 ), where N is
!    the number of knots used, counted according to multiplicity.
!
!    Input, integer ( kind = 4 ) L, min ( M, N - M + 1), where N is the number
!    of knots used, counted according to multiplicity.
!
!    Input, real ( kind = 8 ) V,  the knot in question.
!
!    Input, real ( kind = 8 ) XK(NM), all but the last M entries in the
!    diagonal of K-hat.
!
!    Input, integer ( kind = 4 ) NSTAR, the dimension of the Jacobi matrix.
!
!    Input, real ( kind = 8 ) PHI(NSTAR), the eigenvalues of the Jacobi matrix.
!
!    Input, real ( kind = 8 ) A(NSTAR), the square of the first row of the
!    orthogonal matrix that diagonalizes the Jacobi matrix.
!
!    Input, real ( kind = 8 ) R(L), used to compute the right
!    principal vectors.
!
!    Output, real ( kind = 8 ) D(M), the weights.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nm
  integer ( kind = 4 ) nstar

  real ( kind = 8 ) a(nstar)
  real ( kind = 8 ) d(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) last
  integer ( kind = 4 ) minil
  real ( kind = 8 ) phi(nstar)
  real ( kind = 8 ) r(l)
  real ( kind = 8 ) sum
  real ( kind = 8 ) tmp
  real ( kind = 8 ) v
  real ( kind = 8 ) wf(nstar)
  real ( kind = 8 ) xk(nm)
  real ( kind = 8 ) y(m)
  real ( kind = 8 ) z(m)
!
!  Compute products required for Y-hat.
!
  do j = 1, nstar
    wf(j) = a(j)
    do i = 1, nm
      wf(j) = wf(j) * ( phi(j) - xk(i) )
    end do
  end do
!
!  Compute Y-hat.
!
  do i = 1, m
    sum = 0.0D+00
    do j = 1, nstar
      sum = sum + wf(j)
      wf(j) = wf(j) * ( phi(j) - v )
    end do
    y(i) = sum
  end do
!
!  If N = 1 the right principal vector is already in R.
!  Otherwise compute the R-principal vector of grade M-1.
!
  do i = 1, nm

    tmp = v - xk(i)

    last = min ( l, i + 1 )
    do j = last, 2, -1
      r(j) = tmp * r(j) + r(j-1)
    end do

    r(1) = tmp * r(1)

  end do
!
!  Compute left principal vector(s) and weight for highest derivative.
!  The following statement contains the only division in this
!  routine.  Any test for overflow should be made after it.
!
  d(m) = y(m) / r(1)

  if ( m == 1 ) then
    return
  end if
!
!  Compute left principal vector.
!
  z(1) = 1.0D+00 / r(1)
  do i = 2, m
    sum = 0.0D+00
    minil = min ( i, l )
    do j = 2, minil
      k = i - j + 1
      sum = sum + r(j) * z(k)
    end do
    z(i) = - sum * z(1)
  end do
!
!  Accumulate weights.
!
  do i = 2, m
    sum = 0.0D+00
    do j = 1, i
      k = m - i + j
      sum = sum + z(j) * y(k)
    end do
    k = m - i + 1
    d(k) = sum
  end do

  return
end
subroutine eiqf ( nt, t, mlt, wts, nwts, ndx, key, f, qfsum )

!*****************************************************************************80
!
!! EIQF evaluates an interpolatory quadrature formula.
!
!  Discussion:
!
!   The knots, weights and integrand are supplied.
!
!   All knots with nonzero NDX are used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2008
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, real ( kind = 8 ) WTS(NWTS), the weights.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input, integer ( kind = 4 ) NDX(NT), used to index the array WTS.
!    If KEY = 1, then NDX need not be preset.  For more details see the
!    comments in CAWIQ.
!
!    Input, integer ( kind = 4 ) KEY, indicates the structure of the WTS
!    array.  It will normally be set to 1.  This will cause the weights to be
!    packed sequentially in array WTS.  For more details see the comments
!    in CAWIQ.
!
!    Input, real ( kind = 8 ), external F, the name of a routine which
!    evaluates the function and some of its derivatives.  The routine
!    must have the form
!      function f ( x, i )
!      real ( kind = 8 ) f
!      integer ( kind = 4 ) i
!      real ( kind = 8 ) x
!    and return in F the value of the I-th derivative of the function
!    at X.  The highest value of I will be the maximum value in MLT minus
!    one.  The value X will always be a knot.
!
!    Output, real ( kind = 8 ) QFSUM, the value of the quadrature formula
!    applied to F.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) key
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) ndx(nt)
  real ( kind = 8 ) p
  real ( kind = 8 ) qfsum
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nwts)

  l = abs ( key )

  if ( l < 1 .or. 4 < l ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EIQF - Fatal error!'
    write ( *, '(a)' ) '  Magnitude of KEY must be between 1 and 4.'
    stop
  end if

  qfsum = 0.0D+00
  do j = 1, nt
    l = abs ( ndx(j) )
    if ( l /= 0 ) then
      p = 1.0D+00
      do i = 1, mlt(j)
        qfsum = qfsum + wts(l+i-1) * f ( t(j), i - 1 ) / p
        if ( key <= 0 ) then
          p = p * i
        end if
      end do
    end if
  end do

  return
end
subroutine eiqfs ( nt, t, wts, f, qfsum )

!*****************************************************************************80
!
!! EIQFS evaluates a quadrature formula defined by CLIQF or CLIQFS.
!
!  Discussion:
!
!    This routine evaluates an interpolatory quadrature formula with all knots
!    simple and all knots included in the quadrature.  This routine will be used
!    typically after CLIQF or CLIQFS has been called.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the knots.
!
!    Input, real ( kind = 8 ) WTS(NT), the weights.
!
!    Input, real ( kind = 8 ), external F, the name of a routine which
!    evaluates the function and some of its derivatives.  The routine
!    must have the form
!      function f ( x, i )
!      real ( kind = 8 ) f
!      integer ( kind = 4 ) i
!      real ( kind = 8 ) x
!    and return in F the value of the I-th derivative of the function
!    at X.  The value of I will always be 0.  The value X will always be a knot.
!
!    Output, real ( kind = 8 ) QFSUM, the value of the quadrature formula
!    applied to F.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) j
  real ( kind = 8 ) qfsum
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)

  qfsum = 0.0D+00
  do j = 1, nt
    qfsum = qfsum + wts(j) * f ( t(j), 0 )
  end do

  return
end
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal 
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: itn = 30
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) prec
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end
subroutine parchk ( kind, m, alpha, beta )

!*****************************************************************************80
!
!! PARCHK checks parameters ALPHA and BETA for classical weight functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, integer ( kind = 4 ) M, the order of the highest moment to
!    be calculated.  This value is only needed when KIND = 8.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters, if required
!    by the value of KIND.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) m
  real ( kind = 8 ) tmp

  if ( kind <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  KIND <= 0.'
    stop
  end if
!
!  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
!
  if ( 3 <= kind .and. kind <= 8 .and. alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  3 <= KIND and ALPHA <= -1.'
    stop
  end if
!
!  Check BETA for Jacobi.
!
  if ( kind == 4 .and. beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  KIND == 4 and BETA <= -1.0.'
    stop
  end if
!
!  Check ALPHA and BETA for rational.
!
  if ( kind == 8 ) then
    tmp = alpha + beta + m + 1.0D+00
    if ( 0.0D+00 <= tmp .or. tmp <= beta ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PARCHK - Fatal error!'
      write ( *, '(a)' ) '  KIND == 8 but condition on ALPHA and BETA fails.'
      stop
    end if
  end if

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none

  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = 8 ), parameter :: eps = 2.22D-16
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 171.624D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.0D+30
  real ( kind = 8 ), parameter :: xminin = 2.23D-308
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end
subroutine scmm ( m, kind, alpha, beta, a, b, w )

!*****************************************************************************80
!
!! SCMM computes moments of a classical weight function scaled to [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of moments.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
!    Output, real ( kind = 8 ) W(M), the scaled moments.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a
  real ( kind = 8 ) al
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) be
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) temp
  real ( kind = 8 ) tmp
  real ( kind = 8 ) w(m)

  temp = epsilon ( temp )

  if ( kind == 1 ) then

    al = 0.0D+00
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCMM - Fatal error!'
      write ( *, '(a)' ) '  B - A too small!'
      stop
    end if

    q = ( b - a ) / 2.0D+00
    p = q**( al + be + 1.0D+00 )

  else if ( kind == 2 ) then

    al = -0.5D+00
    be = -0.5D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCMM - Fatal error!'
      write ( *, '(a)' ) '  B - A too small!'
      stop
    end if

    q = ( b - a ) / 2.0D+00
    p = q**( al + be + 1.0D+00 )

  else if ( kind == 3 ) then

    al = alpha
    be = alpha

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCMM - Fatal error!'
      write ( *, '(a)' ) '  B - A too small!'
      stop
    end if

    q = ( b - a ) / 2.0D+00
    p = q**( al + be + 1.0D+00 )

  else if ( kind == 4 ) then

    al = alpha
    be = beta

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCMM - Fatal error!'
      write ( *, '(a)' ) '  B - A too small!'
      stop
    end if

    q = ( b - a ) / 2.0D+00
    p = q**( al + be + 1.0D+00 )

  else if ( kind == 5 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCMM - Fatal error!'
      write ( *, '(a)' ) '  B <= 0!'
      stop
    end if

    q = 1.0D+00 / b
    p = q**( alpha + 1.0D+00 )

  else if ( kind == 6 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCMM - Fatal error!'
      write ( *, '(a)' ) '  B <= 0!'
      stop
    end if

    q = 1.0D+00 / sqrt ( b )
    p = q**( alpha + 1.0D+00 )

  else if ( kind == 7 ) then

    al = alpha
    be = 0.0D+00

    if ( abs ( b - a ) <= temp) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCMM - Fatal error!'
      write ( *, '(a)' ) '  B - A too small!'
      stop
    end if

    q = ( b - a ) / 2.0D+00
    p = q**( al + be + 1.0D+00 )

  else if ( kind == 8 ) then

    if ( a + b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCMM - Fatal error!'
      write ( *, '(a)' ) '  A + B <= 0'
      stop
    end if

    q = a + b
    p = q**( alpha + beta + 1.0D+00 )

  else if ( kind == 9 ) then

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCMM - Fatal error!'
      write ( *, '(a)' ) '  B - A too small!'
      stop
    end if

    q = ( b - a ) / 2.0D+00
    p = q * q

  end if
!
!  Compute the moments in W.
!
  call wm ( m, kind, alpha, beta, w )

  tmp = p

  do i = 1, m
    w(i) = w(i) * tmp
    tmp = tmp * q
  end do

  return
end
subroutine scqf ( nt, t, mlt, wts, nwts, ndx, swts, st, kind, alpha, beta, a, &
  b )

!*****************************************************************************80
!
!! SCQF scales a quadrature formula to a nonstandard interval.
!
!  Discussion:
!
!    The arrays WTS and SWTS may coincide.
!
!    The arrays T and ST may coincide.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the original knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, real ( kind = 8 ) WTS(NWTS), the weights.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input, integer ( kind = 4 ) NDX(NT), used to index the array WTS.
!    For more details see the comments in CAWIQ.
!
!    Output, real ( kind = 8 ) SWTS(NWTS), the scaled weights.
!
!    Output, real ( kind = 8 ) ST(NT), the scaled knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( kind = 8 ) a
  real ( kind = 8 ) al
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) be
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) ndx(nt)
  real ( kind = 8 ) p
  real ( kind = 8 ) shft
  real ( kind = 8 ) slp
  real ( kind = 8 ) st(nt)
  real ( kind = 8 ) swts(nwts)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) temp
  real ( kind = 8 ) tmp
  real ( kind = 8 ) wts(nwts)

  temp = epsilon ( temp )

  call parchk ( kind, 1, alpha, beta )

  if ( kind == 1 ) then

    al = 0.0D+00
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 2 ) then

    al = -0.5D+00
    be = -0.5D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 3 ) then

    al = alpha
    be = alpha

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 4 ) then

    al = alpha
    be = beta

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 5 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  B <= 0'
      stop
    end if

    shft = a
    slp = 1.0D+00 / b
    al = alpha
    be = 0.0D+00

  else if ( kind == 6 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  B <= 0.'
      stop
    end if

    shft = a
    slp = 1.0D+00 / sqrt ( b )
    al = alpha
    be = 0.0D+00

  else if ( kind == 7 ) then

    al = alpha
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 8 ) then

    if ( a + b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  A + B <= 0.'
      stop
    end if

    shft = a
    slp = a + b
    al = alpha
    be = beta

  else if ( kind == 9 ) then

    al = 0.5D+00
    be = 0.5D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  end if

  p = slp**( al + be + 1.0D+00 )

  do k = 1, nt

    st(k) = shft + slp * t(k)
    l = abs ( ndx(k) )

    if ( l /= 0 ) then
      tmp = p
      do i = l, l + mlt(k) - 1
        swts(i) = wts(i) * tmp
        tmp = tmp * slp
      end do
    end if

  end do

  return
end
subroutine sct ( nt, t, kind, a, b, st )

!*****************************************************************************80
!
!! SCT rescales distinct knots to an interval [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the original knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints for which the
!    knots ST should be scaled.
!
!    Output, real ( kind = 8 ) ST(NT), the scaled knots.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  real ( kind = 8 ) shft
  real ( kind = 8 ) slp
  real ( kind = 8 ) st(nt)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) tmp

  if ( kind < 1 .or. 9 < kind ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCT - Fatal error!'
    write ( *, '(a)' ) '  KIND falls outside range of 1 to 9.'
    stop
  end if

  if ( kind == 1 .or. kind == 2 .or. kind == 3 .or. kind == 4 &
    .or. kind == 7 .or. kind == 9 ) then

    tmp = epsilon ( tmp )
    bma = b - a

    if ( bma <= tmp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCT - Fatal error!'
      write ( *, '(a)' ) '  B - A too small.'
      stop
    end if

    slp = 2.0D+00 / bma
    shft = - ( a + b ) / bma

  else if ( kind == 5 ) then

    if ( b < 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCT - Fatal error!'
      write ( *, '(a)' ) '  B < 0.'
      stop
    end if

    slp = b
    shft = - a * b

  else if ( kind == 6 ) then

    if ( b < 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCT - Fatal error!'
      write ( *, '(a)' ) '  B < 0.'
      stop
    end if

    slp = sqrt ( b )
    shft = - a * slp

  else if ( kind == 8 ) then

    slp = 1.0D+00 / ( a + b )

    if ( slp <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCT - Fatal error.'
      write ( *, '(a)' ) '  1 / ( A + B ) <= 0.'
      stop
    end if

    shft = - a * slp

  end if

  st(1:nt) = shft + slp * t(1:nt)

  return
end
subroutine sgqf ( nt, aj, bj, zemu, t, wts )

!*****************************************************************************80
!
!! SGQF computes knots and weights of a Gauss Quadrature formula.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with simple knots from the Jacobi matrix and the zero-th
!    moment of the weight function, using the Golub-Welsch technique.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) AJ(NT), the diagonal of the Jacobi matrix.
!
!    Input/output, real ( kind = 8 ) BJ(NT), the subdiagonal of the Jacobi
!    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
!
!    Input, real ( kind = 8 ) ZEMU, the zero-th moment of the weight function.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) aj(nt)
  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) zemu
!
!  Exit if the zero-th moment is not positive.
!
  if ( zemu <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGQF - Fatal error!'
    write ( *, '(a)' ) '  ZEMU <= 0.'
    stop
  end if
!
!  Set up vectors for IMTQLX.
!
  t(1:nt) = aj(1:nt)

  wts(1) = sqrt ( zemu )
  wts(2:nt) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt)**2

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
subroutine wm ( m, kind, alpha, beta, w )

!*****************************************************************************80
!
!! WM evaluates the first M moments of classical weight functions.
!
!  Discussion:
!
!    W(K) = Integral ( A <= X <= B ) X**(K-1) * W(X) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of moments to evaluate.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) W(M), the first M moments.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) alpha
  real ( kind = 8 ) als
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ja
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kind
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) rk
  real ( kind = 8 ) sum
  real ( kind = 8 ) tmpa
  real ( kind = 8 ) tmpb
  real ( kind = 8 ) trm
  real ( kind = 8 ) w(m)

  call parchk ( kind, m, alpha, beta )

  do k = 2, m, 2
    w(k) = 0.0D+00
  end do

  if ( kind == 1 ) then

    do k = 1, m, 2
      rk = real ( k, kind = 8 )
      w(k) = 2.0D+00 / rk
    end do

  else if ( kind == 2 ) then

    w(1) = pi
    do k = 3, m, 2
      rk = real ( k, kind = 8 )
      w(k) = w(k-2) * ( rk - 2.0D+00 ) / ( rk - 1.0D+00 )
    end do

  else if ( kind == 3 ) then

    w(1) = sqrt ( pi ) * r8_gamma ( alpha + 1.0D+00 ) &
      / r8_gamma ( alpha + 3.0D+00 / 2.0D+00 )

    do k = 3, m, 2
      rk = real ( k, kind = 8 )
      w(k) = w(k-2) * ( rk - 2.0D+00 ) / ( 2.0D+00 * alpha + rk )
    end do

  else if ( kind == 4 ) then

    als = alpha + beta + 1.0D+00
    w(1) = 2.0D+00**als * r8_gamma ( alpha + 1.0D+00 ) &
      / r8_gamma ( als + 1.0D+00 ) * r8_gamma ( beta + 1.0D+00 )

    do k = 2, m

      sum = 0.0D+00
      trm = 1.0D+00
      rk = real ( k, kind = 8 )

      do i = 0, ( k - 2 ) / 2

        tmpa = trm
        do ja = 1, 2 * i
          tmpa = tmpa * ( alpha + ja ) / ( als + ja )
        end do

        do jb = 1, k - 2 * i - 1
          tmpa = tmpa * ( beta + jb ) / ( als + 2 * i + jb )
        end do

        tmpa = tmpa / ( 2 * i + 1.0D+00 ) * &
          ( 2 * i * ( beta + alpha ) + beta - ( rk - 1.0D+00 ) * alpha ) &
          / ( beta + rk - 2 * i - 1.0D+00 )
        sum = sum + tmpa

        trm = trm * ( rk - 2 * i - 1.0D+00 ) &
          / ( 2 * i + 1.0D+00 ) * ( rk - 2 * i - 2.0D+00 ) / ( 2 * i + 2.0D+00 )

      end do

      if ( mod ( k, 2 ) /= 0 ) then
        tmpb = 1.0D+00
        do i = 1, k - 1
          tmpb = tmpb * ( alpha + i ) / ( als + i )
        end do
        sum = sum + tmpb
      end if

      w(k) = sum * w(1)

    end do

  else if ( kind == 5 ) then

    w(1) = r8_gamma ( alpha + 1.0D+00 )

    do k = 2, m
      rk = real ( k, kind = 8 )
      w(k) = ( alpha + rk - 1.0D+00 ) * w(k-1)
    end do

  else if ( kind == 6 ) then

    w(1) = r8_gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )

    do k = 3, m, 2
      rk = real ( k, kind = 8 )
      w(k) = w(k-2) * ( alpha + rk - 2.0D+00 ) / 2.0D+00
    end do

  else if ( kind == 7 ) then

    als = alpha
    do k = 1, m, 2
      rk = real ( k, kind = 8 )
      w(k) = 2.0D+00 / ( rk + als )
    end do

  else if ( kind == 8 ) then

    w(1) = r8_gamma ( alpha + 1.0D+00 ) &
      * r8_gamma ( - alpha - beta - 1.0D+00 ) &
      / r8_gamma ( - beta )

    do k = 2, m
      rk = real ( k, kind = 8 )
      w(k) = - w(k-1) * ( alpha + rk - 1.0D+00 ) / ( alpha + beta + rk )
    end do

  else if ( kind == 9 ) then

    w(1) = pi / 2.0D+00

    do k = 3, m, 2
      rk = real ( k, kind = 8 )
      w(k) = w(k-2) * ( rk - 2.0D+00 ) / ( rk + 1.0D+00 )
    end do

  end if

  return
end
subroutine wtfn ( t, nt, kind, alpha, beta, w )

!*****************************************************************************80
!
!! WTFN evaluates the classical weight functions at given points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(NT), the points where the weight function
!    is to be evaluated.
!
!    Input, integer ( kind = 4 ) NT, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) W(NT), the value of the weight function.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) kind
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) w(nt)

  call parchk ( kind, 1, alpha, beta )

  if ( kind == 1 ) then

    w(1:nt) = 1.0D+00

  else if ( kind == 2 ) then

    w(1:nt) = 1.0D+00 / sqrt ( ( 1.0D+00 - t(1:nt) ) * ( 1.0D+00 + t(1:nt) ) )

  else if ( kind == 3 ) then

    if ( alpha == 0.0D+00 ) then
      w(1:nt) = 1.0D+00
    else
      w(1:nt) = ( ( 1.0D+00 - t(1:nt) ) * ( 1.0D+00 + t(1:nt) ) )**alpha
    end if

  else if ( kind == 4 ) then

    if ( alpha == 0.0D+00 ) then
      w(1:nt) = 1.0D+00
    else
      w(1:nt) = ( 1.0D+00 - t(1:nt) )**alpha
    end if

    if ( beta /= 0.0D+00 ) then
      w(1:nt) = w(1:nt) * ( 1.0D+00 + t(1:nt) )**beta
    end if

  else if ( kind == 5 ) then

    if ( alpha == 0.0D+00 ) then
      w(1:nt) = exp ( - t(1:nt) )
    else
      w(1:nt) = exp ( - t(1:nt) ) * t(1:nt)**alpha
    end if

  else if ( kind == 6 ) then

    if ( alpha == 0.0D+00 ) then
      w(1:nt) = exp ( - t(1:nt)**2 )
    else
      w(1:nt) = exp ( - t(1:nt)**2 ) * abs ( t(1:nt) )**alpha
    end if

  else if ( kind == 7 ) then

    if ( alpha /= 0.0D+00 ) then
      w(1:nt) = abs ( t(1:nt) )**alpha
    else
      w(1:nt) = 1.0D+00
    end if

  else if ( kind == 8 ) then

    if ( alpha == 0.0D+00 ) then
      w(1:nt) = 1.0D+00
    else
      w(1:nt) = t(1:nt)**alpha
    end if

    if ( beta /= 0.0D+00 ) then
      w(1:nt) = w(1:nt) * ( 1.0D+00 + t(1:nt) )**beta
    end if

  else if ( kind == 9 ) then

    w(1:nt) = sqrt ( ( 1.0D+00 - t(1:nt) ) * ( 1.0D+00 + t(1:nt) ) )

  end if

  return
end
