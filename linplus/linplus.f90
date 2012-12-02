function c8_le_l2 ( x, y )

!*****************************************************************************80
!
!! C8_LE_L2 := X <= Y for the L2 norm on C8 values.
!
!  Discussion:
!
!    The L2 norm can be defined here as:
!
!      C8_NORM2(X) = sqrt ( ( real (X) )^2 + ( imag (X) )^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, Y, the values to be compared.
!
!    Output, logical C8_LE_L2, is TRUE if X <= Y.
!
  implicit none

  complex ( kind = 8 ) x
  complex ( kind = 8 ) y
  logical c8_le_l2

  if ( ( real ( x, kind = 8 ) )**2 + ( imag ( x ) )**2 <= &
       ( real ( y, kind = 8 ) )**2 + ( imag ( y ) )**2 ) then
    c8_le_l2 = .true.
  else
    c8_le_l2 = .false.
  end if

  return
end
subroutine c8_swap ( x, y )

!*****************************************************************************80
!
!! C8_SWAP swaps two C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, complex ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  complex ( kind = 8 ) x
  complex ( kind = 8 ) y
  complex ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine c83_cr_fa ( n, a, a_cr )

!*****************************************************************************80
!
!! C83_CR_FA decomposes a C83 matrix using cyclic reduction.
!
!  Discussion:
!
!    The D83 storage format is used for a real tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    Once C83_CR_FA has decomposed a matrix A, then C83_CR_SL may be used 
!    to solve linear systems A * x = b.
!
!    C83_CR_FA does not employ pivoting.  Hence, the results can be more
!    sensitive to ill-conditioning than standard Gauss elimination.  In
!    particular, C83_CR_FA will fail if any diagonal element of the matrix
!    is zero.  Other matrices may also cause C83_CR_FA to fail.
!
!    C83_CR_FA can be guaranteed to work properly if the matrix is strictly
!    diagonally dominant, that is, if the absolute value of the diagonal
!    element is strictly greater than the sum of the absolute values of
!    the offdiagonal elements, for each equation.
!
!    The algorithm may be illustrated by the following figures:
!
!    The initial matrix is given by:
!
!          D1 U1
!          L1 D2 U2
!             L2 D3 U3
!                L3 D4 U4
!                   L4 D5 U5
!                      L5 D6
!
!    Rows and columns are permuted in an odd/even way to yield:
!
!          D1       U1
!             D3    L2 U3
!                D5    L4 U5
!          L1 U2    D2
!             L3 U4    D4
!                L5       D6
!
!    A block LU decomposition is performed to yield:
!
!          D1      |U1
!             D3   |L2 U3
!                D5|   L4 U5
!          --------+--------
!                  |D2'F3
!                  |F1 D4'F4
!                  |   F2 D6'
!
!    For large systems, this reduction is repeated on the lower right hand
!    tridiagonal subsystem until a completely upper triangular system
!    is obtained.  The system has now been factored into the product of a
!    lower triangular system and an upper triangular one, and the information
!    defining this factorization may be used by C83_CR_SL to solve linear
!    systems.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex ( kind = 8 ) A(3,N), the matrix.
!
!    Output, complex ( kind = 8 ) A_CR(3,0:2*N), factorization information 
!    needed by C83_CR_SL.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) a_cr(3,0:2*n)
  integer ( kind = 4 ) iful
  integer ( kind = 4 ) ifulp
  integer ( kind = 4 ) ihaf
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ilp
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) ipnt
  integer ( kind = 4 ) ipntp

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C83_CR_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive N = ', n
    stop
  end if

  if ( n == 1 ) then
    a_cr(1,0:2) = 0.0D+00
    a_cr(2,0) = 0.0D+00
    a_cr(2,1) = 1.0D+00 / a(2,1)
    a_cr(2,2) = 0.0D+00
    a_cr(3,0:2) = 0.0D+00
    return
  end if
!
!  Zero out the workspace entries.
!
  a_cr(1,0) = 0.0D+00
  a_cr(1,1:n-1) = a(1,2:n)
  a_cr(1,n:2*n) = 0.0D+00

  a_cr(2,0) = 0.0D+00
  a_cr(2,1:n) = a(2,1:n)
  a_cr(2,n+1:2*n) = 0.0D+00

  a_cr(3,0) = 0.0D+00
  a_cr(3,1:n-1) = a(3,1:n-1)
  a_cr(3,n:2*n) = 0.0D+00

  il = n
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    if ( mod ( il, 2 ) == 1 ) then
      inc = il + 1
    else
      inc = il
    end if

    incr = inc / 2
    il = il / 2
    ihaf = ipntp + incr + 1
    ifulp = ipnt + inc + 2

!dir$ ivdep
    do ilp = incr, 1, -1
      ifulp = ifulp - 2
      iful = ifulp - 1
      ihaf = ihaf - 1
      a_cr(2,iful) = 1.0D+00 / a_cr(2,iful)
      a_cr(3,iful)  = a_cr(3,iful)  * a_cr(2,iful)
      a_cr(1,ifulp) = a_cr(1,ifulp) * a_cr(2,ifulp+1)
      a_cr(2,ihaf)  = a_cr(2,ifulp) - a_cr(1,iful)  * a_cr(3,iful) &
                                  - a_cr(1,ifulp) * a_cr(3,ifulp)
      a_cr(3,ihaf) = -a_cr(3,ifulp) * a_cr(3,ifulp+1)
      a_cr(1,ihaf) = -a_cr(1,ifulp) * a_cr(1,ifulp+1)
    end do

  end do

  a_cr(2,ipntp+1) = 1.0D+00 / a_cr(2,ipntp+1)

  return
end
subroutine c83_cr_sls ( n, a_cr, nb, b, x )

!*****************************************************************************80
!
!! C83_CR_SLS solves several linear systems factored by C83_CR_FA.
!
!  Discussion:
!
!    The matrix A must be tridiagonal.  C83_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.  If
!    the factors computed by C83_CR_FA are passed to C83_CR_SL, then one or many
!    linear systems involving the matrix A may be solved.
!
!    Note that C83_CR_FA does not perform pivoting, and so the solution 
!    produced by C83_CR_SL may be less accurate than a solution produced 
!    by a standard Gauss algorithm.  However, such problems can be 
!    guaranteed not to occur if the matrix A is strictly diagonally 
!    dominant, that is, if the absolute value of the diagonal coefficient 
!    is greater than the sum of the absolute values of the two off diagonal 
!    coefficients, for each row of the matrix.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex ( kind = 8 ) A_CR(3,0:2*N), factorization information 
!    computed by C83_CR_FA.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input, real ( kind = 8 ) B(N,NB), the right hand sides.
!
!    Output, real ( kind = 8 ) X(N,NB), the solutions of the linear systems.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  complex ( kind = 8 ) a_cr(3,0:2*n)
  complex ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) iful
  integer ( kind = 4 ) ifulm
  integer ( kind = 4 ) ihaf
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ipnt
  integer ( kind = 4 ) ipntp
  integer ( kind = 4 ) ndiv
  complex ( kind = 8 ) rhs(0:2*n,nb)
  complex ( kind = 8 ) x(n,nb)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C83_CR_SLS - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive N = ', n
    stop
  end if

  if ( n == 1 ) then
    x(1,1:nb) = a_cr(2,1) * b(1,1:nb)
    return
  end if
!
!  Set up RHS.
!
  rhs(0,1:nb) = 0.0D+00
  rhs(1:n,1:nb) = b(1:n,1:nb)
  rhs(n+1:2*n,1:nb) = 0.0D+00

  il = n
  ndiv = 1
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    il = il / 2
    ndiv = ndiv * 2
    ihaf = ipntp

!dir$ ivdep
    do iful = ipnt+2, ipntp, 2
      ihaf = ihaf + 1
      rhs(ihaf,1:nb) = rhs(iful,1:nb) &
        - a_cr(3,iful-1) * rhs(iful-1,1:nb) &
        - a_cr(1,iful)   * rhs(iful+1,1:nb)
    end do

  end do

  rhs(ihaf,1:nb) = rhs(ihaf,1:nb) * a_cr(2,ihaf)
  ipnt = ipntp

  do while ( 0 < ipnt )

    ipntp = ipnt
    ndiv = ndiv / 2
    il = n / ndiv
    ipnt = ipnt - il
    ihaf = ipntp

!dir$ ivdep
    do ifulm = ipnt+1, ipntp, 2
      iful = ifulm + 1
      ihaf = ihaf + 1
      rhs(iful,1:nb) = rhs(ihaf,1:nb)
      rhs(ifulm,1:nb) = a_cr(2,ifulm) &
        * (                     rhs(ifulm,1:nb) &
            - a_cr(3,ifulm-1) * rhs(ifulm-1,1:nb) &
            - a_cr(1,ifulm)   * rhs(iful,1:nb) )
    end do

  end do

  x(1:n,1:nb) = rhs(1:n,1:nb)

  return
end
subroutine c83_indicator ( n, a )

!*****************************************************************************80
!
!! C83_INDICATOR sets up a C83 indicator matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Output, complex ( kind = 8 ) A(3,N), the indicator matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  a(1,1) = 0.0D+00
  do j = 2, n
    i = j - 1
    a(1,j) = cmplx ( i, j, kind = 8 )
  end do

  do j = 1, n
    i = j
    a(2,j) = cmplx ( i, j, kind = 8 )
  end do

  do j = 1, n-1
    i = j + 1
    a(3,j) = cmplx ( i, j, kind = 8 )
  end do
  a(3,n) = 0.0D+00

  return
end
subroutine c83_jac_sl ( n, a, b, x, tol, it_max, job, it, diff )

!*****************************************************************************80
!
!! C83_JAC_SL solves a C83 system using Jacobi iteration.
!
!  Discussion:
!
!    The C83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex ( kind = 8 ) A(3,N), the matrix.
!
!    Input, complex ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Input/output, complex ( kind = 8 ) X(N), an approximate solution 
!    to the system.
!
!    Input, real ( kind = 8 ) TOL, a tolerance.  If the maximum change in
!    the solution is less than TOL, the iteration is terminated early.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations taken.
!
!    Output, real ( kind = 8 ) DIFF, the maximum change in the solution
!    on the last iteration.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) b(n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) job
  real ( kind = 8 ) tol
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) x_new(n)
  real ( kind = 8 ) x_norm
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C83_JAC_SL - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero diagonal entry, index = ', i
      stop
    end if
  end do

  if ( job == 0 ) then

    do it_num = 1, it_max

      it = it_num

      x_new(1) =   b(1)                   - a(3,1) * x(2)
      do i = 2, n - 1
        x_new(i) = b(i) - a(1,i) * x(i-1) - a(3,i) * x(i+1)
      end do
      x_new(n) =   b(n) - a(1,n) * x(n-1)
!
!  Divide by diagonal terms.
!
      x_new(1:n) = x_new(1:n) / a(2,1:n)
!
!  Measure norms of solution and change:
!
      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( x_new(1:n) - x(1:n) ) )
!
!  Update.
!
      x(1:n) = x_new(1:n)
!
!  Test for early termination.
!
      if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
        exit
      end if

    end do

  else

    do it_num = 1, it_max

      it = it_num

      x_new(1) =   b(1)                     - a(1,2) * x(2)
      do i = 2, n - 1
        x_new(i) = b(i) - a(3,i-1) * x(i-1) - a(1,i+1) * x(i+1)
      end do
      x_new(n) =   b(n) - a(3,n-1) * x(n-1)
!
!  Divide by diagonal terms.
!
      x_new(1:n) = x_new(1:n) / a(2,1:n)
!
!  Measure norms of solution and change:
!
      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( x_new(1:n) - x(1:n) ) )
!
!  Update:
!
      x(1:n) = x_new(1:n)
!
!  Test for early termination.
!
      if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
        exit
      end if

    end do

  end if

  return
end
subroutine c83_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! C83_MXV multiplies a C83 matrix times a C8VEC.
!
!  Discussion:
!
!    The C83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, complex ( kind = 8 ) A(3,N), the matrix.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, complex ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) x(n)

  b(1:n)   =            a(2,1:n)   * x(1:n)
  b(1:n-1) = b(1:n-1) + a(1,2:n)   * x(2:n)
  b(2:n)   = b(2:n)   + a(3,1:n-1) * x(1:n-1)

  return
end
subroutine c83_np_det ( n, a_lu, det )

!*****************************************************************************80
!
!! C83_NP_DET returns the determinant of a C83 system factored by C83_NP_FA.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex ( kind = 8 ) A_LU(3,N), the LU factors computed 
!    by C83_NP_FA.
!
!    Output, complex ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a_lu(3,n)
  complex ( kind = 8 ) det

  det = product ( a_lu(2,1:n) )

  return
end
subroutine c83_np_fa ( n, a, info )

!*****************************************************************************80
!
!! C83_NP_FA factors a C83 matrix without pivoting.
!
!  Discussion:
!
!    The C83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    Because this routine does not use pivoting, it can fail even when
!    the matrix is not singular, and it is liable to make larger
!    errors.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input/output, complex ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.  On output, factorization information.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  info = 0

  do i = 1, n - 1

    if ( a(2,i) == 0.0D+00 ) then
      info = i
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C83_NP_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Store the multiplier in L.
!
    a(3,i) = a(3,i) / a(2,i)
!
!  Modify the diagonal entry in the next column.
!
    a(2,i+1) = a(2,i+1) - a(3,i) * a(1,i+1)

  end do

  if ( a(2,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C83_NP_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop
  end if

  return
end
subroutine c83_np_fs ( n, a, b, x )

!*****************************************************************************80
!
!! C83_NP_FS factors and solves a C83 system.
!
!  Discussion:
!
!    The C83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, complex ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, complex ( kind = 8 ) B(N), the right hand side.
!
!    Output, complex ( kind = 8 ) X(N), the solution.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C83_NP_FS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
      stop
    end if
  end do

  x(1:n) = b(1:n)

  do i = 2, n
    xmult = a(3,i-1) / a(2,i-1)
    a(2,i) = a(2,i) - xmult * a(1,i)
    x(i)   = x(i)   - xmult * x(i-1)
  end do

  x(n) = x(n) / a(2,n)
  do i = n-1, 1, -1
    x(i) = ( x(i) - a(1,i+1) * x(i+1) ) / a(2,i)
  end do

  return
end
subroutine c83_np_ml ( n, a_lu, x, b, job )

!*****************************************************************************80
!
!! C83_NP_ML computes A * x or x * A, where A has been factored by C83_NP_FA.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex ( kind = 8 ) A(3,N), the LU factors from C83_FA.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, complex ( kind = 8 ) B(N), the product.
!
!    Input, integer ( kind = 4 ) JOB, specifies the product to find.
!    0, compute A * x.
!    nonzero, compute A' * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a_lu(3,n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  complex ( kind = 8 ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Compute X := U * X
!
    do i = 1, n

      b(i) = a_lu(2,i) * b(i)

      if ( i < n ) then
        b(i) = b(i) + a_lu(1,i+1) * b(i+1)
      end if

    end do
!
!  Compute X: = L * X.
!
    do i = n, 2, -1
      b(i) = b(i) + a_lu(3,i-1) * b(i-1)
    end do

  else
!
!  Compute X: = L' * X.
!
    do i = 1, n-1
      b(i) = b(i) + a_lu(3,i) * b(i+1)
    end do
!
!  Compute X: = U' * X.
!
    do i = n, 2, -1
      b(i) = a_lu(2,i) * b(i)
      b(i) = b(i) + a_lu(1,i) * b(i-1)
    end do
    b(1) = a_lu(2,1) * b(1)

  end if

  return
end
subroutine c83_np_sl ( n, a_lu, b, job )

!*****************************************************************************80
!
!! C83_NP_SL solves a C83 system factored by C83_NP_FA.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex ( kind = 8 ) A(3,N), the factor information from C83_NP_FA.
!
!    Input/output, complex ( kind = 8 ) B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution of the linear system.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a_lu(3,n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job

  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    do i = 2, n
      b(i) = b(i) - a_lu(3,i-1) * b(i-1)
    end do
!
!  Solve U * X = Y.
!
    do i = n, 1, -1
      b(i) = b(i) / a_lu(2,i)
      if ( 1 < i ) then
        b(i-1) = b(i-1) - a_lu(1,i) * b(i)
      end if
    end do

  else
!
!  Solve U' * Y = B
!
    do i = 1, n
      b(i) = b(i) / a_lu(2,i)
      if ( i < n ) then
        b(i+1) = b(i+1) - a_lu(1,i+1) * b(i)
      end if
    end do
!
!  Solve L' * X = Y.
!
    do i = n-1, 1, -1
      b(i) = b(i) - a_lu(3,i) * b(i+1)
    end do

  end if

  return
end
subroutine c83_print ( n, a, title )

!*****************************************************************************80
!
!! C83_PRINT prints a C83 matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
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
!    Input, complex ( kind = 8 ) A(3,N), the C83 matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  character ( len = * ) title

  call c83_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine c83_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C83_PRINT_SOME prints some of a C83 matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
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
!    Input, complex ( kind = 8 ) A(3,N), the C83 matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column, to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 3
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  character ( len = 12 ) citemp(incx)
  character ( len = 12 ) crtemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  real ( kind = 8 ) xi
  real ( kind = 8 ) xr
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( crtemp(j2), '(i8,6x)' ) j
      write ( citemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col:  '',6a12)' ) ( crtemp(j2), citemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - 1 )
    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + 1 )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( 1 < i-j .or. 1 < j-i ) then

          crtemp(j2) = ' '
          citemp(j2) = ' '

        else

          if ( j == i-1 ) then
            xr = real ( a(1,i), kind = 8 )
            xi = imag ( a(1,i) )
          else if ( j == i ) then
            xr = real ( a(2,i), kind = 8 )
            xi = imag ( a(2,i) )
          else if ( j == i+1 ) then
            xr = real ( a(3,i), kind = 8 )
            xi = imag ( a(3,i) )
          end if

          if ( xr == 0.0D+00 .and. xi == 0.0D+00 ) then
            crtemp(j2) = '    0.0'
            citemp(j2) = ' '
          else if ( xr == 0.0D+00 .and. xi /= 0.0D+00 ) then
            crtemp(j2) = ' '
            write ( citemp(j2), '(g12.5)' ) xi
          else if ( xr /= 0.0D+00 .and. xi == 0.0D+00 ) then
            write ( crtemp(j2), '(g12.5)' ) xr
            citemp(j2) = ' '
          else
            write ( crtemp(j2), '(g12.5)' ) xr
            write ( citemp(j2), '(g12.5)' ) xi
          end if

        end if

      end do

      write ( *, '(i5,1x,6a12)' ) i, ( crtemp(j2), citemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c83_random ( n, seed, a )

!*****************************************************************************80
!
!! C83_RANDOM randomizes a C83 matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, complex ( kind = 8 ) A(3,N), the C83 matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) seed

  a(1,1) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  call c8vec_uniform_01 ( n - 1, seed, a(1,2:n) )

  call c8vec_uniform_01 ( n,     seed, a(2,1:n) )

  call c8vec_uniform_01 ( n - 1, seed, a(3,1:n-1) )
  a(3,n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  return
end
subroutine c83_to_c8ge ( n, a1, a2 )

!*****************************************************************************80
!
!! C83_TO_C8GE copies a C83 matrix into a C8GE matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    The C8GE storage format is used for a complex general M by N matrix.  
!    A storage space is made for each logical entry.  The two 
!    dimensional logical array is mapped to a vector, in which storage 
!    is by columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex ( kind = 8 ) A1(3,N), the C83 matrix.
!
!    Output, complex ( kind = 8 ) A2(N,N), the C8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a1(3,n)
  complex ( kind = 8 ) a2(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n

      if ( j == i - 1 ) then
        a2(i,j) = a1(1,i)
      else if ( i == j ) then
        a2(i,j) = a1(2,i)
      else if ( j == i + 1 ) then
        a2(i,j) = a1(3,i)
      else
        a2(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine c83_vxm ( n, a, x, b )

!*****************************************************************************80
!
!! C83_VXM multiplies a C8VEC by a C83 matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, complex ( kind = 8 ) A(3,N), the C83 matrix.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A'.
!
!    Output, complex ( kind = 8 ) B(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) x(n)

  b(1:n)   =            a(2,1:n)   * x(1:n)
  b(1:n-1) = b(1:n-1) + a(3,1:n-1) * x(2:n) 
  b(2:n)   = b(2:n)   + a(1,2:n)   * x(1:n-1)

  return
end
subroutine c8ci_eval ( n, a, lambda )

!*****************************************************************************80
!
!! C8CI_EVAL returns the eigenvalues of a C8CI matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis,
!    Circulant Matrices,
!    Second Edition,
!    Chelsea, 1994,
!    ISBN: 0828403384,
!    LC: QA188.D37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Output, complex ( kind = 8 ) LAMBDA(N), the eigenvalues.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  complex ( kind = 8 ) lambda(n)
  complex ( kind = 8 ) w(n)

  call c8vec_unity ( n, w )

  lambda(1:n) = a(n)
  do i = n-1, 1, -1
    lambda(1:n) = lambda(1:n) * w(1:n) + a(i)
  end do

  call c8vec_sort_a2 ( n, lambda )

  return
end
subroutine c8ci_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! C8CI_MXV multiplies a C8CI matrix times a C8VEC.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, complex ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  do i = 1, n
    b(i) = zero
    do j = 1, i-1
      b(i) = b(i) + a(n+j+1-i) * x(j)
    end do
    do j = i, n
      b(i) = b(i) + a(j+1-i) * x(j)
    end do
  end do

  return
end
subroutine c8ci_print ( n, a, title )

!*****************************************************************************80
!
!! C8CI_PRINT prints a C8CI matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  character ( len = * ) title

  call c8ci_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine c8ci_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8CI_PRINT_SOME prints some of a C8CI matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
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
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) aij
  character ( len = 20 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  complex ( kind = 8 ) zero
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(j+1-i)
        else
          aij = a(n+j+1-i)
        end if

        if ( aij == zero ) then
          ctemp(j2) = '     0.0            '
        else if ( imag ( aij ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( aij, kind = 8 )
        else
          write ( ctemp(j2), '(2g10.3)' ) aij
        end if

      end do

      write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8ci_random ( n, seed, a )

!*****************************************************************************80
!
!! C8CI_RANDOM randomizes a C8CI matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
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
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, complex ( kind = 8 ) A(N), the C8CI matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) seed

  call c8vec_uniform_01 ( n, seed, a )

  return
end
subroutine c8ci_sl ( n, a, b, x, job )

!*****************************************************************************80
!
!! C8CI_SL solves a C8CI system.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Input, complex ( kind = 8 ) B(N), the right hand side.
!
!    Output, complex ( kind = 8 ) X(N), the solution of the linear system.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  integer ( kind = 4 ) nsub
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2
  complex ( kind = 8 ) r3
  complex ( kind = 8 ) r5
  complex ( kind = 8 ) r6
  complex ( kind = 8 ) work(2*n-2)
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  if ( job == 0 ) then
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = zero
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(n+2-nsub)
      r6 = a(nsub)

      if ( 2 < nsub ) then

        work(nsub-1) = r2

        do i = 1, nsub-2
          r5 = r5 + a(n+1-i) * work(nsub-i)
          r6 = r6 + a(i+1) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( 2 < nsub ) then

        r6 = work(n)
        work(n-1+nsub-1) = zero
        do i = 2, nsub-1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = zero
      do i = 1, nsub-1
        r5 = r5 + a(n+1-i) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      do i = 1, nsub-1
        x(i) = x(i) + work(n-1+i) * r6
      end do

      x(nsub) = r6

    end do

  else
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = zero
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(nsub)
      r6 = a(n+2-nsub)

      if ( 2 < nsub ) then

        work(nsub-1) = r2

        do i = 1, nsub-2
          r5 = r5 + a(i+1) * work(nsub-i)
          r6 = r6 + a(n+1-i) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( 2 < nsub ) then

        r6 = work(n)
        work(n-1+nsub-1) = zero
        do i = 2, nsub-1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = zero
      do i = 1, nsub-1
        r5 = r5 + a(i+1) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      do i = 1, nsub-1
        x(i) = x(i) + work(n-1+i) * r6
      end do

      x(nsub) = r6

    end do

  end if

  return
end
subroutine c8ci_to_c8ge ( n, a, b )

!*****************************************************************************80
!
!! C8CI_TO_C8GE copies a C8CI matrix to a C8GE matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!    The C8GE storage format is used for a complex general M by N matrix.  
!    A storage space is made for each logical entry.  The two 
!    dimensional logical array is mapped to a vector, in which storage 
!    is by columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Output, complex ( kind = 8 ) B(N,N), the C8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, i-1
      b(i,j) = a(n+j+1-i)
    end do
    do j = i, n
      b(i,j) = a(j+1-i)
    end do
  end do

  return
end
subroutine c8ci_vxm ( n, a, x, b )

!*****************************************************************************80
!
!! C8CI_VXM multiplies a C8VEC by a C8CI matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Input, complex X(N), the vector to be multiplied by A.
!
!    Output, complex B(N), the product A' * X.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  b(1:n) = zero

  do i = 1, n
    do j = 1, i
      b(i) = b(i) + a(i+1-j) * x(j)
    end do
    do j = i+1, n
      b(i) = b(i) + a(n+i+1-j) * x(j)
    end do
  end do

  return
end
subroutine c8ge_random ( m, n, seed, a )

!*****************************************************************************80
!
!! C8GE_RANDOM randomizes a C8GE matrix.
!
!  Discussion:
!
!    The C8GE storage format is used for a complex general M by N matrix.  
!    A storage space is made for each logical entry.  The two 
!    dimensional logical array is mapped to a vector, in which storage 
!    is by columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, complex ( kind = 8 ) A(M,N), the C8GE matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) seed

  call c8mat_uniform_01 ( m, n, seed, a )

  return
end
subroutine c8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! C8MAT_PRINT prints a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call c8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8MAT_PRINT_SOME prints some of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  character ( len = 20 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == zero ) then
          ctemp(j2) = '       0.0          '
        else if ( imag ( a(i,j) ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = 8 )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8mat_uniform_01 ( m, n, seed, c )

!*****************************************************************************80
!
!! C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
!
!  Discussion:
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should 
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, double complex C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta

  do j = 1, n
    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      theta = 2.0D+00 * pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

    end do

  end do

  return
end
subroutine c8to_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! C8TO_MXV multiplies a C8TO matrix times a C8VEC.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(2*N-1), the entries of the first row of 
!    the Toeplitz matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, complex ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  do i = 1, n

    b(i) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

    do j = 1, i-1
      b(i) = b(i) + a(n+i-j) * x(j)
    end do

    do j = i, n
      b(i) = b(i) + a(j+1-i) * x(j)
    end do

  end do

  return
end
subroutine c8to_print ( n, a, title )

!*****************************************************************************80
!
!! C8TO_PRINT prints a C8TO matrix.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 March 2001
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
!    Input, complex ( kind = 8 ) A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  character ( len = * ) title

  call c8to_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine c8to_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8TO_PRINT_SOME prints some of a C8TO matrix.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 March 2001
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
!    Input, complex ( kind = 8 ) A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) aij
  character ( len = 20 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  complex ( kind = 8 ) zero
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(j+1-i)
        else
          aij = a(n+i-j)
        end if

        if ( aij == zero ) then
          ctemp(j2) = '    0.0'
        else if ( imag ( aij ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( aij, kind = 8 )
        else
          write ( ctemp(j2), '(2g10.3)' ) aij
        end if

      end do

      write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8to_random ( n, seed, a )

!*****************************************************************************80
!
!! C8TO_RANDOM randomizes a C8TO matrix.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
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
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, complex ( kind = 8 ) A(2*N-1), the randomized matrix, with 
!    entries between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) seed

  n2 = 2 * n - 1

  call c8vec_uniform_01 ( n2, seed, a )

  return
end
subroutine c8to_sl ( n, a, b, x, job )

!*****************************************************************************80
!
!! C8TO_SL solves a C8TO system.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(2*N-1), the first row of the Toeplitz
!    matrix, followed by the first column of the Toeplitz matrix, beginning 
!    with the second element.
!
!    Input, complex ( kind = 8 ) B(N) the right hand side vector.
!
!    Output, complex ( kind = 8 ) X(N), the solution vector.  X and B may 
!    share the same storage.
!
!    Input, integer ( kind = 4 ) JOB,
!    0 to solve A*X=B,
!    nonzero to solve A'*X=B.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) c1(n-1)
  complex ( kind = 8 ) c2(n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  integer ( kind = 4 ) nsub
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2
  complex ( kind = 8 ) r3
  complex ( kind = 8 ) r5
  complex ( kind = 8 ) r6
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  if ( n < 1 ) then
    return
  end if
!
!  Solve the system with the principal minor of order 1.
!
  r1 = a(1)
  x(1) = b(1) / r1

  if ( n == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system with the Toeplitz matrix.
!
  do nsub = 2, n
!
!  Compute multiples of the first and last columns of the inverse of
!  the principal minor of order NSUB.
!
    if ( job == 0 ) then
      r5 = a(n+nsub-1)
      r6 = a(nsub)
    else
      r5 = a(nsub)
      r6 = a(n+nsub-1)
    end if

    if ( 2 < nsub ) then

      c1(nsub-1) = r2

      do i = 1, nsub-2
        if ( job == 0 ) then
          r5 = r5 + a(n+i) * c1(nsub-i)
          r6 = r6 + a(i+1) * c2(i)
        else
          r5 = r5 + a(i+1) * c1(nsub-i)
          r6 = r6 + a(n+i) * c2(i)
        end if
      end do

    end if

    r2 = -r5 / r1
    r3 = -r6 / r1
    r1 = r1 + r5 * r3

    if ( 2 < nsub ) then

      r6 = c2(1)
      c2(nsub-1) = zero
      do i = 2, nsub-1
        r5 = c2(i)
        c2(i) = c1(i) * r3 + r6
        c1(i) = c1(i) + r6 * r2
        r6 = r5
      end do

    end if

    c2(1) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
    r5 = zero

    do i = 1, nsub-1
      if ( job == 0 ) then
        r5 = r5 + a(n+i) * x(nsub-i)
      else
        r5 = r5 + a(i+1) * x(nsub-i)
      end if
    end do

    r6 = ( b(nsub) - r5 ) / r1

    do i = 1, nsub-1
      x(i) = x(i) + c2(i) * r6
    end do

    x(nsub) = r6

  end do

  return
end
subroutine c8to_to_c8ge ( n, a, b )

!*****************************************************************************80
!
!! C8TO_TO_C8GE copies a C8TO matrix to a C8GE matrix.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!    The C8GE storage format is used for a complex general M by N matrix.  
!    A storage space is made for each logical entry.  The two 
!    dimensional logical array is mapped to a vector, in which storage 
!    is by columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(2*N-1), the C8TO matrix.
!
!    Output, complex ( kind = 8 ) B(N,N), the C8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, i-1
      b(i,j) = a(n+i-j)
    end do
    do j = i, n
      b(i,j) = a(j-i+1)
    end do
  end do

  return
end
subroutine c8to_vxm ( n, a, x, b )

!*****************************************************************************80
!
!! C8TO_VXM multiplies a C8VEC by a C8TO matrix.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(2*N-1), the entries of the first row of the
!    Toeplitz matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, complex ( kind = 8 ) B(N), the product A' * X.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  do i = 1, n

    b(i) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

    do j = 1, i
      b(i) = b(i) + a(i+1-j) * x(j)
    end do

    do j = i+1, n
      b(i) = b(i) + a(n+j-i) * x(j)
    end do

  end do

  return
end
subroutine c8vec_indicator ( n, a )

!*****************************************************************************80
!
!! C8VEC_INDICATOR sets a C8VEC to an "indicator" vector.
!
!  Discussion:
!
!    X(1:N) = ( 1-1i, 2-2i, 3-3i, 4-4i, ... )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, complex ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = cmplx ( i, -i, kind = 8 )
  end do

  return
end
subroutine c8vec_print ( n, a, title )

!*****************************************************************************80
!
!! C8VEC_PRINT prints a C8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, complex ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,2g14.6)' ) i, a(i)
  end do

  return
end
subroutine c8vec_print_some ( n, x, max_print, title )

!*****************************************************************************80
!
!! C8VEC_PRINT_SOME prints some of a C8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines 
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title
  complex ( kind = 8 ) x(n)

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i8,2x,2g14.6)' ) i, x(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print-2
      write ( *, '(i8,2x,2g14.6)' ) i, x(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i8,2x,2g14.6)' ) i, x(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i8,2x,2g14.6)' ) i, x(i)
    end do
    i = max_print
    write ( *, '(i8,2x,2g14.6,2x,a)' ) i, x(i), '...more entries...'

  end if

  return
end
subroutine c8vec_sort_a2 ( n, x )

!*****************************************************************************80
!
!! C8VEC_SORT_A2 ascending sorts a C8VEC by L2 norm.
!
!  Discussion:
!
!    The L2 norm of A+Bi is sqrt ( A**2 + B**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input/output, complex ( kind = 8 ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  logical c8_le_l2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call c8_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c8_le_l2 ( x(i), x(j) ) ) then
        isgn = - 1
      else
        isgn = + 1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
!
!  Discussion:
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values to compute.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should 
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, double complex C(N), the pseudorandom complex vector.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    theta = 2.0D+00 * pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

  end do

  return
end
subroutine c8vec_unity ( n, a )

!*****************************************************************************80
!
!! C8VEC_UNITY returns the N roots of unity as a C8VEC.
!
!  Discussion:
!
!    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
!
!    X(1:N)**N = ( (1,0), (1,0), ..., (1,0) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, complex ( kind = 8 ) A(N), the N roots of unity.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta

  do i = 1, n
    theta = pi * real ( 2 * ( i - 1 ), kind = 8 ) / real ( n, kind = 8 )
    a(i) = cmplx ( cos ( theta ), sin ( theta ), kind = 8 )
  end do

  return
end
subroutine daxpy ( n, sa, x, incx, y, incy )

!*****************************************************************************80
!
!! DAXPY adds a constant times one R8VEC to another.
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
!    Original FORTRAN77 version by Lawson, Hanson, Kincaid, Krogh.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Charle Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) SA, the multiplier.
!
!    Input, real ( kind = 8 ) X(*), the vector to be scaled and added to Y.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
!    Input/output, real ( kind = 8 ) Y(*), the vector to which a 
!    multiple of X is to be added.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive 
!    entries of Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) n
  real ( kind = 8 ) sa
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

  else if ( sa == 0.0D+00 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    y(1:n) = y(1:n) + sa * x(1:n)

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy  ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      y(iy) = y(iy) + sa * x(ix)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine file_delete ( file_name )

!*****************************************************************************80
!
!! FILE_DELETE deletes a named file if it exists.
!
!  Discussion:
!
!    You might want to call this routine to get rid of any old copy
!    of a file, before trying to open a new copy with the OPEN argument:
!      status = 'new'.
!
!    It's not always safe to open a file with " STATUS = 'UNKNOWN' ".
!    For instance, on the SGI, the most recent version of the FORTRAN
!    compiler seems to go crazy when I open an unformatted direct
!    access file this way.  It creates an enormous file (of somewhat
!    random size).  The problem goes away if I delete any old copy
!    using this routine, and then open a fresh copy with
!    " STATUS = 'NEW' ".  It's a scary world.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
  implicit none

  logical file_exist
  logical file_is_open
  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical, parameter :: verbose = .false.
!
!  Does the file exist?
!
  if ( .not. file_exist ( file_name ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE - Warning!'
    write ( *, '(a)' ) '  There is no file of the given name.'
    return
  end if
!
!  Is the file open?
!
  if ( file_is_open ( file_name ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE - Warning!'
    write ( *, '(a)' ) '  The file "' // trim ( file_name ) &
      // '" is currently open.'
    write ( *, '(a)' ) '  It must be closed before it can be deleted.'
    return
  end if
!
!  Get a free unit number.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  A free FORTRAN unit could not be found.'
    return
  end if

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE:'
    write ( *, '(a)' ) '  Deleting "' // trim ( file_name ) // '".'
  end if

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_DELETE: Warning!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    return
  end if

  close ( unit = iunit, status = 'delete' )

  return
end
function file_exist ( file_name )

!*****************************************************************************80
!
!! FILE_EXIST reports whether a file exists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_EXIST, is TRUE if the file exists.
!
  implicit none

  character ( len = * ) file_name
  logical file_exist

  inquire ( file = file_name, exist = file_exist )

  return
end
function file_is_open ( file_name )

!*****************************************************************************80
!
!! FILE_IS_OPEN reports whether a file (specified by filename) is open.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_IS_OPEN, is TRUE if the file is open.
!
  implicit none

  character ( len = * ) file_name
  logical file_is_open

  inquire ( file = file_name, opened = file_is_open )

  return
end
subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 ) today
  integer ( kind = 4 ) values(8)
  character ( len = 5 ) zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) / real (  11.0, kind = 8 )
  temp = temp + real ( values(3) - 1, kind = 8 ) / real (  30.0, kind = 8 )
  temp = temp + real ( values(5),     kind = 8 ) / real (  23.0, kind = 8 )
  temp = temp + real ( values(6),     kind = 8 ) / real (  59.0, kind = 8 )
  temp = temp + real ( values(7),     kind = 8 ) / real (  59.0, kind = 8 )
  temp = temp + real ( values(8),     kind = 8 ) / real ( 999.0, kind = 8 )
  temp = temp                                    / real (   6.0, kind = 8 )

  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( huge ( 1 ), kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == huge ( 1 ) ) then
    seed = seed - 1
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine hilbert_inverse ( n, a )

!*****************************************************************************80
!
!! HILBERT_INVERSE returns the inverse of the Hilbert matrix.
!
!  Formula:
!
!    A(I,J) =  (-1)**(I+J) * (N+I-1)! * (N+J-1)! /
!           [ (I+J-1) * ((I-1)!*(J-1)!)**2 * (N-I)! * (N-J)! ]
!
!  Example:
!
!    N = 5
!
!       25    -300     1050    -1400     630
!     -300    4800   -18900    26880  -12600
!     1050  -18900    79380  -117600   56700
!    -1400   26880  -117600   179200  -88200
!      630  -12600    56700   -88200   44100
!
!  Properties:
!
!    A is symmetric.
!
!    Because A is symmetric, it is normal, so diagonalizable.
!
!    A is almost impossible to compute accurately by general routines
!    that compute the inverse.
!
!    A is integral.
!
!    The sum of the entries of A is N**2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) A(N,N), the inverse Hilbert matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!
!  Set the (1,1) entry.
!
  a(1,1) = real ( n * n, kind = 8 )
!
!  Define Row 1, Column J by recursion on Row 1 Column J-1
!
  i = 1
  do j = 2, n
    a(i,j) = - a(i,j-1) &
      * real ( ( n + j - 1 ) * ( i + j - 2 ) * ( n + 1 - j ), kind = 8 ) &
      / real ( ( i + j - 1 ) * ( j - 1 ) * ( j - 1 ), kind = 8 )
  end do
!
!  Define Row I by recursion on row I-1
!
  do i = 2, n
    do j = 1, n

      a(i,j) = - a(i-1,j) &
        * real ( (n+i-1) * (i+j-2) * (n+1-i), kind = 8 ) &
        / real ( (i+j-1) * ( i - 1 ) * ( i - 1 ), kind = 8 )

    end do
  end do

  return
end
function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!
!  Example:
!
!        I  I4_LOG_10
!    -----  --------
!        0    0
!        1    0
!        2    0
!        9    0
!       10    1
!       11    1
!       99    1
!      100    2
!      101    2
!      999    2
!     1000    3
!     1001    3
!     9999    3
!    10000    4
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 10 
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_10, the integer part of the 
!    logarithm base 10 of the absolute value of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) ten_pow

  if ( i == 0 ) then

    i4_log_10 = 0

  else

    i4_log_10 = 0
    ten_pow = 10

    i_abs = abs ( i )

    do while ( ten_pow <= i_abs )
      i4_log_10 = i4_log_10 + 1
      ten_pow = ten_pow * 10
    end do

  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) big
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(i8,1x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(i8,1x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i8,i11)' ) i, a(i)
    end do
  end if

  return
end
subroutine i4vec_search_binary_a ( n, a, b, indx )

!*****************************************************************************80
!
!! I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    Binary search is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.9,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the array to be searched.  A must
!    be sorted in ascending order.
!
!    Input, integer ( kind = 4 ) B, the value to be searched for.
!
!    Output, integer ( kind = 4 ) INDX, the result of the search.
!    -1, B does not occur in A.
!    I, A(I) = B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) high
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid

  indx = -1
  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == b ) then
      indx = mid
      exit
    else if ( a(mid) < b ) then
      low = mid + 1
    else if ( b < a(mid) ) then
      high = mid - 1
    end if

  end do

  return
end
function r8_is_int ( r )

!*****************************************************************************80
!
!! R8_IS_INT determines if an R8 represents an integer value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the number to be checked.
!
!    Output, logical R8_IS_INT, is TRUE if R is an integer value.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  logical r8_is_int

  if ( real ( huge ( i ), kind = 8 ) < r ) then
    r8_is_int = .false.
  else if ( r < - real ( huge ( i ), kind = 8 ) ) then
    r8_is_int = .false.
  else if ( r == real ( int ( r ), kind = 8 ) ) then
    r8_is_int = .true.
  else
    r8_is_int = .false.
  end if

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP switches two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
function r8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, 
!    which should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r8_uniform = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, 
!    which should NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r83_cr_fa ( n, a, a_cr )

!*****************************************************************************80
!
!! R83_CR_FA decomposes an R83 matrix using cyclic reduction.
!
!  Discussion:
!
!    The R83 storage format is used for a real tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    Once R83_CR_FA has decomposed a matrix A, then R83_CR_SL may be used 
!    to solve linear systems A * x = b.
!
!    R83_CR_FA does not employ pivoting.  Hence, the results can be more
!    sensitive to ill-conditioning than standard Gauss elimination.  In
!    particular, R83_CR_FA will fail if any diagonal element of the matrix
!    is zero.  Other matrices may also cause R83_CR_FA to fail.
!
!    R83_CR_FA can be guaranteed to work properly if the matrix is strictly
!    diagonally dominant, that is, if the absolute value of the diagonal
!    element is strictly greater than the sum of the absolute values of
!    the offdiagonal elements, for each equation.
!
!    The algorithm may be illustrated by the following figures:
!
!    The initial matrix is given by:
!
!          D1 U1
!          L1 D2 U2
!             L2 D3 U3
!                L3 D4 U4
!                   L4 D5 U5
!                      L5 D6
!
!    Rows and columns are permuted in an odd/even way to yield:
!
!          D1       U1
!             D3    L2 U3
!                D5    L4 U5
!          L1 U2    D2
!             L3 U4    D4
!                L5       D6
!
!    A block LU decomposition is performed to yield:
!
!          D1      |U1
!             D3   |L2 U3
!                D5|   L4 U5
!          --------+--------
!                  |D2'F3
!                  |F1 D4'F4
!                  |   F2 D6'
!
!    For large systems, this reduction is repeated on the lower right hand
!    tridiagonal subsystem until a completely upper triangular system
!    is obtained.  The system has now been factored into the product of a
!    lower triangular system and an upper triangular one, and the information
!    defining this factorization may be used by R83_CR_SL to solve linear
!    systems.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2004
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(3,N), the R83 matrix.
!
!    Output, real ( kind = 8 ) A_CR(3,0:2*N), factorization information 
!    needed by R83_CR_SL.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) a_cr(3,0:2*n)
  integer ( kind = 4 ) iful
  integer ( kind = 4 ) ifulp
  integer ( kind = 4 ) ihaf
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ilp
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) ipnt
  integer ( kind = 4 ) ipntp

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83_CR_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive N = ', n
    stop
  end if

  if ( n == 1 ) then
    a_cr(1,0:2) = 0.0D+00
    a_cr(2,0) = 0.0D+00
    a_cr(2,1) = 1.0D+00 / a(2,1)
    a_cr(2,2) = 0.0D+00
    a_cr(3,0:2) = 0.0D+00
    return
  end if
!
!  Zero out the workspace entries.
!
  a_cr(1,0) = 0.0D+00
  a_cr(1,1:n-1) = a(1,2:n)
  a_cr(1,n:2*n) = 0.0D+00

  a_cr(2,0) = 0.0D+00
  a_cr(2,1:n) = a(2,1:n)
  a_cr(2,n+1:2*n) = 0.0D+00

  a_cr(3,0) = 0.0D+00
  a_cr(3,1:n-1) = a(3,1:n-1)
  a_cr(3,n:2*n) = 0.0D+00

  il = n
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    if ( mod ( il, 2 ) == 1 ) then
      inc = il + 1
    else
      inc = il
    end if

    incr = inc / 2
    il = il / 2
    ihaf = ipntp + incr + 1
    ifulp = ipnt + inc + 2

!dir$ ivdep
    do ilp = incr, 1, -1
      ifulp = ifulp - 2
      iful = ifulp - 1
      ihaf = ihaf - 1
      a_cr(2,iful) = 1.0D+00 / a_cr(2,iful)
      a_cr(3,iful)  = a_cr(3,iful)  * a_cr(2,iful)
      a_cr(1,ifulp) = a_cr(1,ifulp) * a_cr(2,ifulp+1)
      a_cr(2,ihaf)  = a_cr(2,ifulp) - a_cr(1,iful)  * a_cr(3,iful) &
                                  - a_cr(1,ifulp) * a_cr(3,ifulp)
      a_cr(3,ihaf) = -a_cr(3,ifulp) * a_cr(3,ifulp+1)
      a_cr(1,ihaf) = -a_cr(1,ifulp) * a_cr(1,ifulp+1)
    end do

  end do

  a_cr(2,ipntp+1) = 1.0D+00 / a_cr(2,ipntp+1)

  return
end
subroutine r83_cr_sl ( n, a_cr, b, x )

!*****************************************************************************80
!
!! R83_CR_SL solves a linear systems factored by R83_CR_FA.
!
!  Discussion:
!
!    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.  If
!    the factors computed by R83_CR_FA are passed to R83_CR_SL, then 
!    a linear system involving the matrix A may be solved.
!
!    Note that R83_CR_FA does not perform pivoting, and so the solution 
!    produced by R83_CR_SL may be less accurate than a solution produced 
!    by a standard Gauss algorithm.  However, such problems can be 
!    guaranteed not to occur if the matrix A is strictly diagonally 
!    dominant, that is, if the absolute value of the diagonal coefficient 
!    is greater than the sum of the absolute values of the two off diagonal 
!    coefficients, for each row of the matrix.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A_CR(3,0:2*N), factorization information 
!    computed by R83_CR_FA.
!
!    Input, real ( kind = 8 ) B(N), the right hand side.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear systems.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_cr(3,0:2*n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) iful
  integer ( kind = 4 ) ifulm
  integer ( kind = 4 ) ihaf
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ipnt
  integer ( kind = 4 ) ipntp
  integer ( kind = 4 ) ndiv
  real ( kind = 8 ) rhs(0:2*n)
  real ( kind = 8 ) x(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83_CR_SL - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive N = ', n
    stop
  end if

  if ( n == 1 ) then
    x(1) = a_cr(2,1) * b(1)
    return
  end if
!
!  Set up RHS.
!
  rhs(0) = 0.0D+00
  rhs(1:n) = b(1:n)
  rhs(n+1:2*n) = 0.0D+00

  il = n
  ndiv = 1
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    il = il / 2
    ndiv = ndiv * 2
    ihaf = ipntp

!dir$ ivdep
    do iful = ipnt + 2, ipntp, 2
      ihaf = ihaf + 1
      rhs(ihaf) = rhs(iful) &
        - a_cr(3,iful-1) * rhs(iful-1) &
        - a_cr(1,iful)   * rhs(iful+1)
    end do

  end do

  rhs(ihaf) = a_cr(2,ihaf) * rhs(ihaf)

  ipnt = ipntp

  do while ( 0 < ipnt )

    ipntp = ipnt
    ndiv = ndiv / 2
    il = n / ndiv
    ipnt = ipnt - il
    ihaf = ipntp

!dir$ ivdep
    do ifulm = ipnt + 1, ipntp, 2
      iful = ifulm + 1
      ihaf = ihaf + 1
      rhs(iful) = rhs(ihaf)
      rhs(ifulm) = a_cr(2,ifulm) &
        * (                     rhs(ifulm) &
            - a_cr(3,ifulm-1) * rhs(ifulm-1) &
            - a_cr(1,ifulm)   * rhs(iful) )
    end do

  end do

  x(1:n) = rhs(1:n)

  return
end
subroutine r83_cr_sls ( n, a_cr, nb, b, x )

!*****************************************************************************80
!
!! R83_CR_SLS solves several linear systems factored by R83_CR_FA.
!
!  Discussion:
!
!    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.  If
!    the factors computed by R83_CR_FA are passed to R83_CR_SLS, then one or 
!    many linear systems involving the matrix A may be solved.
!
!    Note that R83_CR_FA does not perform pivoting, and so the solution 
!    produced by R83_CR_SLS may be less accurate than a solution produced 
!    by a standard Gauss algorithm.  However, such problems can be 
!    guaranteed not to occur if the matrix A is strictly diagonally 
!    dominant, that is, if the absolute value of the diagonal coefficient 
!    is greater than the sum of the absolute values of the two off diagonal 
!    coefficients, for each row of the matrix.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A_CR(3,0:2*N), factorization information 
!    computed by R83_CR_FA.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input, real ( kind = 8 ) B(N,NB), the right hand sides.
!
!    Output, real ( kind = 8 ) X(N,NB), the solutions of the linear systems.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a_cr(3,0:2*n)
  real ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) iful
  integer ( kind = 4 ) ifulm
  integer ( kind = 4 ) ihaf
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ipnt
  integer ( kind = 4 ) ipntp
  integer ( kind = 4 ) ndiv
  real ( kind = 8 ) rhs(0:2*n,nb)
  real ( kind = 8 ) x(n,nb)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83_CR_SLS - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive N = ', n
    stop
  end if

  if ( n == 1 ) then
    x(1,1:nb) = a_cr(2,1) * b(1,1:nb)
    return
  end if
!
!  Set up RHS.
!
  rhs(0,1:nb) = 0.0D+00
  rhs(1:n,1:nb) = b(1:n,1:nb)
  rhs(n+1:2*n,1:nb) = 0.0D+00

  il = n
  ndiv = 1
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    il = il / 2
    ndiv = ndiv * 2
    ihaf = ipntp

!dir$ ivdep
    do iful = ipnt + 2, ipntp, 2
      ihaf = ihaf + 1
      rhs(ihaf,1:nb) = rhs(iful,1:nb) &
        - a_cr(3,iful-1) * rhs(iful-1,1:nb) &
        - a_cr(1,iful)   * rhs(iful+1,1:nb)
    end do

  end do

  rhs(ihaf,1:nb) = rhs(ihaf,1:nb) * a_cr(2,ihaf)
  ipnt = ipntp

  do while ( 0 < ipnt )

    ipntp = ipnt
    ndiv = ndiv / 2
    il = n / ndiv
    ipnt = ipnt - il
    ihaf = ipntp

!dir$ ivdep
    do ifulm = ipnt + 1, ipntp, 2
      iful = ifulm + 1
      ihaf = ihaf + 1
      rhs(iful,1:nb) = rhs(ihaf,1:nb)
      rhs(ifulm,1:nb) = a_cr(2,ifulm) &
        * (                     rhs(ifulm,1:nb) &
            - a_cr(3,ifulm-1) * rhs(ifulm-1,1:nb) &
            - a_cr(1,ifulm)   * rhs(iful,1:nb) )
    end do

  end do

  x(1:n,1:nb) = rhs(1:n,1:nb)

  return
end
subroutine r83_gs_sl ( n, a, b, x, tol, it_max, job, it, diff )

!*****************************************************************************80
!
!! R83_GS_SL solves an R83 system using Gauss-Seidel iteration.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A(3,N), the R83 matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = 8 ) X(N), an approximate solution to 
!    the system.
!
!    Input, real ( kind = 8 ) TOL, a tolerance.  If the maximum change in
!    the solution is less than TOL, the iteration is terminated early.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations taken.
!
!    Output, real ( kind = 8 ) DIFF, the maximum change in the solution
!    on the last iteration.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) job
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_norm
  real ( kind = 8 ) x_old(n)
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_GS_SL - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero diagonal entry, index = ', i
      stop
    end if
  end do

  if ( job == 0 ) then

    do it_num = 1, it_max

      it = it_num

      x_old(1:n) = x(1:n)

      x(1) =   ( b(1)                   - a(3,1) * x(2)   ) / a(2,1)
      do i = 2, n - 1
        x(i) = ( b(i) - a(1,i) * x(i-1) - a(3,i) * x(i+1) ) / a(2,i)
      end do
      x(n) =   ( b(n) - a(1,n) * x(n-1)                   ) / a(2,n)

      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( x(1:n) - x_old(1:n) ) )

      if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
        exit
      end if

    end do

  else

    do it_num = 1, it_max

      it = it_num

      x_old(1:n) = x(1:n)

      x(1) =   ( b(1)                     - a(1,2) * x(2)     ) / a(2,1)
      do i = 2, n - 1
        x(i) = ( b(i) - a(3,i-1) * x(i-1) - a(1,i+1) * x(i+1) ) / a(2,i)
      end do
      x(n) =   ( b(n) - a(3,n-1) * x(n-1)                     ) / a(2,n)

      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( x(1:n) - x_old(1:n) ) )

      if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
        exit
      end if
   
    end do

  end if

  return
end
subroutine r83_indicator ( n, a )

!*****************************************************************************80
!
!! R83_INDICATOR sets up an R83 indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    Here are the values as stored in an indicator matrix:
!
!      00 12 23 34 45
!      11 22 33 44 55
!      21 32 43 54 00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Output, real ( kind = 8 ) A(3,N), the R83 indicator matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  a(1,1) = 0.0D+00
  do j = 2, n
    i = j - 1
    a(1,j) = real ( fac * i + j, kind = 8 )
  end do

  do j = 1, n
    i = j
    a(2,j) = real ( fac * i + j, kind = 8 )
  end do

  do j = 1, n-1
    i = j + 1
    a(3,j) = real ( fac * i + j, kind = 8 )
  end do
  a(3,n) = 0.0D+00

  return
end
subroutine r83_jac_sl ( n, a, b, x, tol, it_max, job, it, diff )

!*****************************************************************************80
!
!! R83_JAC_SL solves an R83 system using Jacobi iteration.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A(3,N), the R83 matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = 8 ) X(N), an approximate solution 
!    to the system.
!
!    Input, real ( kind = 8 ) TOL, a tolerance.  If the maximum change in
!    the solution is less than TOL, the iteration is terminated early.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations taken.
!
!    Output, real ( kind = 8 ) DIFF, the maximum change in the solution
!    on the last iteration.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) job
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_new(n)
  real ( kind = 8 ) x_norm
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_JAC_SL - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero diagonal entry, index = ', i
      stop
    end if
  end do

  if ( job == 0 ) then

    do it_num = 1, it_max

      it = it_num

      x_new(1) =   b(1)                   - a(3,1) * x(2)
      do i = 2, n - 1
        x_new(i) = b(i) - a(1,i) * x(i-1) - a(3,i) * x(i+1)
      end do
      x_new(n) =   b(n) - a(1,n) * x(n-1)
!
!  Divide by diagonal terms.
!
      x_new(1:n) = x_new(1:n) / a(2,1:n)
!
!  Measure norms of solution and change:
!
      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( x_new(1:n) - x(1:n) ) )
!
!  Update.
!
      x(1:n) = x_new(1:n)
!
!  Test for early termination.
!
      if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
        exit
      end if

    end do

  else

    do it_num = 1, it_max

      it = it_num

      x_new(1) =   b(1)                     - a(1,2) * x(2)
      do i = 2, n - 1
        x_new(i) = b(i) - a(3,i-1) * x(i-1) - a(1,i+1) * x(i+1)
      end do
      x_new(n) =   b(n) - a(3,n-1) * x(n-1)
!
!  Divide by diagonal terms.
!
      x_new(1:n) = x_new(1:n) / a(2,1:n)
!
!  Measure norms of solution and change:
!
      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( x_new(1:n) - x(1:n) ) )
!
!  Update:
!
      x(1:n) = x_new(1:n)
!
!  Test for early termination.
!
      if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
        exit
      end if

    end do

  end if

  return
end
subroutine r83_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! R83_MXV multiplies an R83 matrix times an R8VEC.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, real ( kind = 8 ) A(3,N), the R83 matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  b(1:n)   =            a(2,1:n)   * x(1:n)
  b(1:n-1) = b(1:n-1) + a(1,2:n)   * x(2:n)
  b(2:n)   = b(2:n)   + a(3,1:n-1) * x(1:n-1)

  return
end
subroutine r83_np_det ( n, a_lu, det )

!*****************************************************************************80
!
!! R83_NP_DET returns the determinant of an R83 system factored by R83_NP_FA.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A_LU(3,N), the LU factors computed by R83_NP_FA.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(3,n)
  real ( kind = 8 ) det

  det = product ( a_lu(2,1:n) )

  return
end
subroutine r83_np_fa ( n, a, info )

!*****************************************************************************80
!
!! R83_NP_FA factors an R83 matrix without pivoting.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    Because this routine does not use pivoting, it can fail even when
!    the matrix is not singular, and it is liable to make larger
!    errors.
!
!    R83_NP_FA and R83_NP_SL may be preferable to the corresponding
!    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
!    in one step, and does not save the factorization.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.  On output, factorization information.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  info = 0

  do i = 1, n-1

    if ( a(2,i) == 0.0D+00 ) then
      info = i
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_NP_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Store the multiplier in L.
!
    a(3,i) = a(3,i) / a(2,i)
!
!  Modify the diagonal entry in the next column.
!
    a(2,i+1) = a(2,i+1) - a(3,i) * a(1,i+1)

  end do

  if ( a(2,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83_NP_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop
  end if

  return
end
subroutine r83_np_fs ( n, a, b, x )

!*****************************************************************************80
!
!! R83_NP_FS factors and solves an R83 system.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_NP_FS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
      stop
    end if
  end do

  x(1:n) = b(1:n)

  do i = 2, n
    xmult = a(3,i-1) / a(2,i-1)
    a(2,i) = a(2,i) - xmult * a(1,i)
    x(i)   = x(i)   - xmult * x(i-1)
  end do

  x(n) = x(n) / a(2,n)
  do i = n-1, 1, -1
    x(i) = ( x(i) - a(1,i+1) * x(i+1) ) / a(2,i)
  end do

  return
end
subroutine r83_np_fss ( n, a, nb, b, x )

!*****************************************************************************80
!
!! R83_NP_FSS factors and solves multiple R83 systems.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!  Example:
!
!    Here is how a R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input, real ( kind = 8 ) B(N,NB), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N,NB), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n,nb)
  real ( kind = 8 ) xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_NP_FSS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
      stop
    end if
  end do

  x(1:n,1:nb) = b(1:n,1:nb)

  do i = 2, n
    xmult = a(3,i-1) / a(2,i-1)
    a(2,i) = a(2,i) - xmult * a(1,i)
    x(i,1:nb)   = x(i,1:nb)   - xmult * x(i-1,1:nb)
  end do

  x(n,1:nb) = x(n,1:nb) / a(2,n)
  do i = n-1, 1, -1
    x(i,1:nb) = ( x(i,1:nb) - a(1,i+1) * x(i+1,1:nb) ) / a(2,i)
  end do

  return
end
subroutine r83_np_ml ( n, a_lu, x, b, job )

!*****************************************************************************80
!
!! R83_NP_ML computes A * x or x * A, where A has been factored by R83_NP_FA.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A_LU(3,N), the LU factors from R83_FA.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A*x or A'*x.
!
!    Input, integer ( kind = 4 ) JOB, specifies the product to find.
!    0, compute A * x.
!    nonzero, compute A' * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(3,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  real ( kind = 8 ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Compute X := U * X
!
    do i = 1, n

      b(i) = a_lu(2,i) * b(i)

      if ( i < n ) then
        b(i) = b(i) + a_lu(1,i+1) * b(i+1)
      end if

    end do
!
!  Compute X: = L * X.
!
    do i = n, 2, -1
      b(i) = b(i) + a_lu(3,i-1) * b(i-1)
    end do

  else
!
!  Compute X: = L' * X.
!
    do i = 1, n-1
      b(i) = b(i) + a_lu(3,i) * b(i+1)
    end do
!
!  Compute X: = U' * X.
!
    do i = n, 2, -1
      b(i) = a_lu(2,i) * b(i)
      b(i) = b(i) + a_lu(1,i) * b(i-1)
    end do
    b(1) = a_lu(2,1) * b(1)

  end if

  return
end
subroutine r83_np_sl ( n, a_lu, b, job )

!*****************************************************************************80
!
!! R83_NP_SL solves an R83 system factored by R83_NP_FA.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A_LU(3,N), the LU factors from R83_NP_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution of the linear system.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(3,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job

  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    do i = 2, n
      b(i) = b(i) - a_lu(3,i-1) * b(i-1)
    end do
!
!  Solve U * X = Y.
!
    do i = n, 1, -1
      b(i) = b(i) / a_lu(2,i)
      if ( 1 < i ) then
        b(i-1) = b(i-1) - a_lu(1,i) * b(i)
      end if
    end do

  else
!
!  Solve U' * Y = B
!
    do i = 1, n
      b(i) = b(i) / a_lu(2,i)
      if ( i < n ) then
        b(i+1) = b(i+1) - a_lu(1,i+1) * b(i)
      end if
    end do
!
!  Solve L' * X = Y.
!
    do i = n-1, 1, -1
      b(i) = b(i) - a_lu(3,i) * b(i+1)
    end do

  end if

  return
end
subroutine r83_print ( n, a, title )

!*****************************************************************************80
!
!! R83_PRINT prints an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 September 2003
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
!    Input, real ( kind = 8 ) A(3,N), the R83 matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  character ( len = * ) title

  call r83_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r83_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R83_PRINT_SOME prints some of an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2003
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
!    Input, real ( kind = 8 ) A(3,N), the R83 matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column, to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - 1 )
    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + 1 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( 1 < i-j .or. 1 < j-i ) then
          ctemp(j2) = '              '
        else if ( j == i+1 ) then
          if ( r8_is_int ( a(1,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(1,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(1,j)
          end if
        else if ( j == i ) then
          if ( r8_is_int ( a(2,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(2,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(2,j)
          end if
        else if ( j == i-1 ) then
          if ( r8_is_int ( a(3,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(3,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(3,j)
          end if
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r83_random ( n, seed, a )

!*****************************************************************************80
!
!! R83_RANDOM randomizes an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) A(3,N), the R83 matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) seed

  a(1,1) = 0.0D+00
  call r8vec_uniform_01 ( n - 1, seed, a(1,2:n) )

  call r8vec_uniform_01 ( n,     seed, a(2,1:n) )

  call r8vec_uniform_01 ( n - 1, seed, a(3,1:n-1) )
  a(3,n) = 0.0D+00

  return
end
subroutine r83_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R83_TO_R8GE copies an R83 matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A(3,N), the R83 matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n

      if ( j == i-1 ) then
        b(i,j) = a(3,j)
      else if ( j == i ) then
        b(i,j) = a(2,j)
      else if ( j == i+1 ) then
        b(i,j) = a(1,j)
      else
        b(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine r83_vxm ( n, a, x, b )

!*****************************************************************************80
!
!! R83_VXM multiplies an R8VEC by an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, real ( kind = 8 ) A(3,N), the R83 matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) B(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  b(1:n)   =            a(2,1:n)   * x(1:n)
  b(1:n-1) = b(1:n-1) + a(3,1:n-1) * x(2:n) 
  b(2:n)   = b(2:n)   + a(1,2:n)   * x(1:n-1)

  return
end
subroutine r83np_fs ( n, a, b, x )

!*****************************************************************************80
!
!! R83NP_FS factors and solves an R83NP system.
!
!  Discussion:
!
!    The R83NP storage format is used for a tridiagonal matrix.
!    The subdiagonal   is in entries (1,2:N), 
!    the diagonal      is in entries (2,1:N), 
!    the superdiagonal is in entries (3,1:N-1). 
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!    The "R83NP" format used for this routine is different from the R83 format.
!    Here, we insist that the nonzero entries
!    for a given row now appear in the corresponding column of the
!    packed array.
!
!  Example:
!
!    Here is how an R83NP matrix of order 5 would be stored:
!
!       *  A21 A32 A43 A54
!      A11 A22 A33 A44 A55
!      A12 A23 A34 A45  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83NP_FS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
      stop
    end if
  end do

  x(1:n) = b(1:n)

  do i = 2, n
    a(2,i) = a(2,i) - a(3,i-1) * a(1,i) / a(2,i-1)
    x(i)   = x(i)   - x(i-1)   * a(1,i) / a(2,i-1)
  end do

  x(n) = x(n) / a(2,n)
  do i = n-1, 1, -1
    x(i) = ( x(i) - a(3,i) * x(i+1) ) / a(2,i)
  end do

  return
end
subroutine r83p_det ( n, a_lu, work4, det )

!*****************************************************************************80
!
!! R83P_DET computes the determinant of a matrix factored by R83P_FA.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = 8 ) A_LU(3,N), LU factors from R83P_FA.
!
!    Input, real ( kind = 8 ) WORK4, factorization information from R83P_FA.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(3,n)
  real ( kind = 8 ) det
  real ( kind = 8 ) work4

  det = product ( a_lu(2,1:n-1) ) * work4

  return
end
subroutine r83p_fa ( n, a, info, work2, work3, work4 )

!*****************************************************************************80
!
!! R83P_FA factors an R83P matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!    Once the matrix has been factored by R83P_FA, R83P_SL may be called
!    to solve linear systems involving the matrix.
!
!    The logical matrix has a form which is suggested by this diagram:
!
!      D1 U1          L1
!      L2 D2 U2
!         L3 D3 U3
!            L4 D4 U4
!               L5 D5 U5
!      U6          L6 D6
!
!    The algorithm treats the matrix as a border banded matrix:
!
!      ( A1  A2 )
!      ( A3  A4 )
!
!    where:
!
!      D1 U1          | L1
!      L2 D2 U2       |  0
!         L3 D3  U3    |  0
!            L4 D4 U4 |  0
!               L5 D5 | U5
!      ---------------+---
!      U6  0  0  0 L6 | D6
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Method:
!
!    The algorithm rewrites the system as:
!
!         X1 + inverse(A1) A2 X2 = inverse(A1) B1
!
!      A3 X1 +             A4 X2 = B2
!
!    The first equation can be "solved" for X1 in terms of X2:
!
!         X1 = - inverse(A1) A2 X2 + inverse(A1) B1
!
!    allowing us to rewrite the second equation for X2 explicitly:
!
!      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 3.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the periodic tridiagonal matrix.  
!    On output, the arrays have been modified to hold information
!    defining the border-banded factorization of submatrices A1
!    and A3.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
!    Output, real ( kind = 8 ) WORK2(N-1), WORK3(N-1), WORK4, 
!    factorization information.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  real ( kind = 8 ) work2(n-1)
  real ( kind = 8 ) work3(n-1)
  real ( kind = 8 ) work4
!
!  Compute inverse(A1):
!
  call r83_np_fa ( n - 1, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83P_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  R83_NP_FA returned INFO = ', info
    write ( *, '(a)' ) '  Factoring failed for column INFO.'
    write ( *, '(a)' ) '  The tridiagonal matrix A1 is singular.'
    write ( *, '(a)' ) '  This algorithm cannot continue!'
    stop
  end if
!
!  WORK2 := inverse(A1) * A2.
!
  work2(1) = a(3,n)
  work2(2:n-2) = 0.0D+00
  work2(n-1) = a(1,n)

  job = 0
  call r83_np_sl ( n - 1, a, work2, job )
!
!  WORK3 := inverse ( A1' ) * A3'.
!
  work3(1) = a(1,1)
  work3(2:n-2) = 0.0D+00
  work3(n-1) = a(3,n-1)

  job = 1
  call r83_np_sl ( n - 1, a, work3, job )
!
!  A4 := ( A4 - A3 * inverse(A1) * A2 )
!
  work4 = a(2,n) - a(1,1) * work2(1) - a(3,n-1) * work2(n-1)

  if ( work4 == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83P_FA - Fatal error!'
    write ( *, '(a)' ) '  The factored A4 submatrix is zero.'
    write ( *, '(a)' ) '  This algorithm cannot continue!'
    stop
  end if

  return
end
subroutine r83p_indicator ( n, a )

!*****************************************************************************80
!
!! R83P_INDICATOR sets up an R83P indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!    Here are the values as stored in an indicator matrix:
!
!      51 12 23 34 45
!      11 22 33 44 55
!      21 32 43 54 15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Output, real ( kind = 8 ) A(3,N), the R83P indicator matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  i = n
  j = 1
  a(1,j) = real ( fac * i + j, kind = 8 )
  do j = 2, n
    i = j - 1
    a(1,j) = real ( fac * i + j, kind = 8 )
  end do

  do j = 1, n
    i = j
    a(2,j) = real ( fac * i + j, kind = 8 )
  end do

  do j = 1, n-1
    i = j + 1
    a(3,j) = real ( fac * i + j, kind = 8 )
  end do
  i = 1
  j = n
  a(3,j) = real ( fac * i + j, kind = 8 )

  return
end
subroutine r83p_ml ( n, a_lu, x, b, job )

!*****************************************************************************80
!
!! R83P_ML computes A * x or x * A, where A has been factored by R83P_FA.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = 8 ) A_LU(3,N), the LU factors ffrom R83P_FA.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by the matrix.
!
!    Output, real ( kind = 8 ) B(N), the result of the multiplication.
!
!    Input, integer ( kind = 4 ) JOB, indicates what product should be computed.
!    0, compute A * x.
!    nonzero, compute A' * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(3,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) job
  real ( kind = 8 ) x(n)
!
!  Multiply A(1:N-1,1:N-1) and X(1:N-1).
!
  call r83_np_ml ( n - 1, a_lu, x, b, job )
!
!  Add terms from the border.
!
  if ( job == 0 ) then
    b(1) = b(1) + a_lu(3,n) * x(n)
    b(n-1) = b(n-1) + a_lu(1,n) * x(n)
    b(n) = a_lu(1,1) * x(1) + a_lu(3,n-1) * x(n-1) + a_lu(2,n) * x(n)
  else
    b(1) = b(1) + a_lu(1,1) * x(n)
    b(n-1) = b(n-1) + a_lu(3,n-1) * x(n)
    b(n) = a_lu(3,n) * x(1) + a_lu(1,n) * x(n-1) + a_lu(2,n) * x(n)
  end if

  return
end
subroutine r83p_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! R83P_MXV multiplies an R83P matrix by an R8VEC.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = 8 ) A(3,N), the R83P matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  b(1) =   a(3,n)   * x(n)   + a(2,1) * x(1) + a(1,2)   * x(2)

  do i = 2, n-1
    b(i) = a(3,i-1) * x(i-1) + a(2,i) * x(i) + a(1,i+1) * x(i+1)
  end do

  b(n) =   a(3,n-1) * x(n-1) + a(2,n) * x(n) + a(1,1)   * x(1)

  return
end
subroutine r83p_print ( n, a, title )

!*****************************************************************************80
!
!! R83P_PRINT prints an R83P matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
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
!    Input, real ( kind = 8 ) A(3,N), the R83P matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  character ( len = * ) title

  call r83p_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r83p_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R83P_PRINT_SOME prints some of an R83P matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2004
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
!    Input, real ( kind = 8 ) A(3,N), the R83P matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column, to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the column range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    if ( 1 < i2lo .or. j2hi < n ) then
      i2lo = max ( i2lo, j2lo - 1 )
    end if

    i2hi = min ( ihi, n )

    if ( i2hi < n .or. 1 < j2lo ) then
      i2hi = min ( i2hi, j2hi + 1 )
    end if

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i == n .and. j == 1 ) then
          if ( r8_is_int ( a(1,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(1,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(1,j)
          end if
        else if ( i == 1 .and. j == n ) then
          if ( r8_is_int ( a(3,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(3,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(3,j)
          end if
        else if ( 1 < i - j .or. 1 < j - i ) then
          ctemp(j2) = '              '
        else if ( j == i + 1 ) then
          if ( r8_is_int ( a(1,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(1,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(1,j)
          end if
        else if ( j == i ) then
          if ( r8_is_int ( a(2,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(2,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(2,j)
          end if
        else if ( j == i - 1 ) then
          if ( r8_is_int ( a(3,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(3,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(3,j)
          end if
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r83p_random ( n, seed, a )

!*****************************************************************************80
!
!! R83P_RANDOM randomizes an R83P matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 3.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) A(3,N), the R83P matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) seed

  call r8mat_uniform_01 ( 3, n, seed, a )

  return
end
subroutine r83p_sl ( n, a_lu, b, x, job, work2, work3, work4 )

!*****************************************************************************80
!
!! R83P_SL solves an R83P system.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!    The linear system must have been factored by R83P_FA.
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = 8 ) A_LU(3,N), the LU factors from R83P_FA.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution to the linear system.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
!    Input, real ( kind = 8 ) WORK2(N-1), WORK3(N-1), WORK4, 
!    factor data from R83P_FA.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(3,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) job
  real ( kind = 8 ) work2(n-1)
  real ( kind = 8 ) work3(n-1)
  real ( kind = 8 ) work4
  real ( kind = 8 ) x(n)

  x(1:n) = b(1:n)

  if ( job == 0 ) then
!
!  Solve A1 * X1 = B1.
!
    call r83_np_sl ( n - 1, a_lu, x, job )
!
!  X2 = B2 - A3 * X1
!
    x(n) = x(n) - a_lu(1,1) * x(1) - a_lu(3,n-1) * x(n-1)
!
!  Solve A4 * X2 = X2
!
    x(n) = x(n) / work4
!
!  X1 := X1 - inverse ( A1 ) * A2 * X2.
!
    x(1:n-1) = x(1:n-1) - work2(1:n-1) * x(n)

  else
!
!  Solve A1' * X1 = B1.
!
    call r83_np_sl ( n - 1, a_lu, x, job )
!
!  X2 := X2 - A2' * B1
!
    x(n) = x(n) - a_lu(3,n) * x(1) - a_lu(1,n) * x(n-1)
!
!  Solve A4 * X2 = X2.
!
    x(n) = x(n) / work4
!
!  X1 := X1 - transpose ( inverse ( A1 ) * A3 ) * X2.
!
    x(1:n-1) = x(1:n-1) - work3(1:n-1) * x(n)

  end if

  return
end
subroutine r83p_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R83P_TO_R8GE copies an R83P matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = 8 ) A(3,N), the R83P matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        b(i,j) = a(2,j)
      else if ( j == i-1 ) then
        b(i,j) = a(3,j)
      else if ( j == i+1 ) then
        b(i,j) = a(1,j)
      else if ( i == 1 .and. j == n ) then
        b(i,j) = a(3,j)
      else if ( i == n .and. j == 1 ) then
        b(i,j) = a(1,j)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r83p_vxm ( n, a, x, b )

!*****************************************************************************80
!
!! R83P_VXM multiplies an R8VEC by an R83P matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = 8 ) A(3,N), the R83P matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product X * A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  b(1) =   a(1,1)   * x(n)   + a(2,1) * x(1) + a(3,1)   * x(2)

  do i = 2, n-1
    b(i) = a(1,i)   * x(i-1) + a(2,i) * x(i) + a(3,i)   * x(i+1)
  end do

  b(n) =   a(1,n)   * x(n-1) + a(2,n) * x(n) + a(3,n)   * x(1)

  return
end
subroutine r85_indicator ( n, a )

!*****************************************************************************80
!
!! R85_INDICATOR sets up an R85 indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!    Here are the values as stored in an indicator matrix:
!
!      00 00 13 24 35 46
!      00 12 23 34 45 56
!      11 22 33 44 55 66
!      21 32 43 54 65 00
!      31 42 53 64 00 00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 3.
!
!    Output, real ( kind = 8 ) A(5,N), the indicator matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(5,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  a(1,1) = 0.0D+00
  a(1,2) = 0.0D+00
  do j = 3, n
    i = j - 2
    a(1,j) = real ( fac * i + j, kind = 8 )
  end do

  a(2,1) = 0.0D+00
  do j = 2, n
    i = j - 1
    a(2,j) = real ( fac * i + j, kind = 8 )
  end do

  do j = 1, n
    i = j
    a(3,j) = real ( fac * i + j, kind = 8 )
  end do

  do j = 1, n-1
    i = j + 1
    a(4,j) = real ( fac * i + j, kind = 8 )
  end do
  a(4,n) = 0.0D+00

  do j = 1, n-2
    i = j + 2
    a(5,j) = real ( fac * i + j, kind = 8 )
  end do
  a(5,n-1) = 0.0D+00
  a(5,n) = 0.0D+00

  return
end
subroutine r85_np_fs ( n, a, b, x )

!*****************************************************************************80
!
!! R85_NP_FS factors and solves an R85 linear system.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!    The factorization algorithm requires that each diagonal entry be nonzero.
!
!    No pivoting is performed, and therefore the algorithm may fail
!    in simple cases where the matrix is not singular.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 September 2003
!
!  Author:
!
!    Original FORTRAN77 version by Ward Cheney, David Kincaid.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Ward Cheney, David Kincaid,
!    Numerical Mathematics and Computing,
!    Brooks-Cole Publishing, 2004,
!    ISBN: 0534201121.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, real ( kind = 8 ) A(5,N),
!    On input, the pentadiagonal matrix.
!    On output, the data has been overwritten by factorization information.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B has been overwritten by factorization information.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmult

  do i = 1, n
    if ( a(3,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R85_NP_FS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(3,', i, ') = 0.'
      stop
    end if
  end do

  do i = 2, n-1

    xmult = a(2,i) / a(3,i-1)
    a(3,i) = a(3,i) - xmult * a(4,i-1)
    a(4,i) = a(4,i) - xmult * a(5,i-1)

    b(i) = b(i) - xmult * b(i-1)

    xmult = a(1,i+1) / a(3,i-1)
    a(2,i+1) = a(2,i+1) - xmult * a(4,i-1)
    a(3,i+1) = a(3,i+1) - xmult * a(5,i-1)

    b(i+1) = b(i+1) - xmult * b(i-1)

  end do

  xmult = a(2,n) / a(3,n-1)
  a(3,n) = a(3,n) - xmult * a(4,n-1)

  x(n) = ( b(n) - xmult * b(n-1) ) / a(3,n)
  x(n-1) = ( b(n-1) - a(4,n-1) * x(n) ) / a(3,n-1)

  do i = n-2, 1, -1
    x(i) = ( b(i) - a(4,i) * x(i+1) - a(5,i) * x(i+2) ) / a(3,i)
  end do

  return
end
subroutine r85_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! R85_MXV multiplies an R85 matrix by an R8VEC.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, real ( kind = 8 ) A(5,N), the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  b(1:n)   =            a(3,1:n)   * x(1:n)

  b(3:n)   = b(3:n)   + a(1,3:n)   * x(1:n-2)
  b(2:n)   = b(2:n)   + a(2,2:n)   * x(1:n-1)
  b(1:n-1) = b(1:n-1) + a(4,1:n-1) * x(2:n)
  b(1:n-2) = b(1:n-2) + a(5,1:n-2) * x(3:n)

  return
end
subroutine r85_print ( n, a, title )

!*****************************************************************************80
!
!! R85_PRINT prints an R85 matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 September 2003
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
!    Input, real ( kind = 8 ) A(5,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(5,n)
  character ( len = * ) title

  call r85_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r85_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R85_PRINT_SOME prints some of an R85 matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2004
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
!    Input, real ( kind = 8 ) A(5,N), the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column, to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(5,n)
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - 2 )

    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + 2 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( 2 < i-j .or. 2 < j-i ) then
          ctemp(j2) = '              '
        else if ( j == i+2 ) then
          if ( r8_is_int ( a(1,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(1,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(1,j)
          end if
        else if ( j == i+1 ) then
          if ( r8_is_int ( a(2,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(2,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(2,j)
          end if
        else if ( j == i ) then
          if ( r8_is_int ( a(3,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(3,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(3,j)
          end if
        else if ( j == i-1 ) then
          if ( r8_is_int ( a(4,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(4,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(4,j)
          end if
        else if ( j == i-2 ) then
          if ( r8_is_int ( a(5,j) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a(5,j)
          else
            write ( ctemp(j2), '(g14.6)' ) a(5,j)
          end if
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r85_random ( n, seed, a )

!*****************************************************************************80
!
!! R85_RANDOM randomizes an R85 matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) A(5,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  a(1,1) = 0.0D+00
  a(1,2) = 0.0D+00
  do j = 3, n
    a(1,j) = r8_uniform_01 ( seed )
  end do

  a(2,1) = 0.0D+00
  do j = 2, n
    a(2,j) = r8_uniform_01 ( seed )
  end do

  do j = 1, n
    a(3,j) = r8_uniform_01 ( seed )
  end do

  do j = 1, n-1
    a(4,j) = r8_uniform_01 ( seed )
  end do

  do j = 1, n-2
    a(5,j) = r8_uniform_01 ( seed )
  end do

  return
end
subroutine r85_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R85_TO_R8GE copies an R85 matrix into an R8GE matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = 8 ) A(5,N), the R85 matrix.
!
!    Output, real ( kind = 8 ) A(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n

      if ( j == i-2 ) then
        b(i,j) = a(1,i)
      else if ( j == i-1 ) then
        b(i,j) = a(2,i)
      else if ( i == j ) then
        b(i,j) = a(3,i)
      else if ( j == i+1 ) then
        b(i,j) = a(4,i)
      else if ( j == i+2 ) then
        b(i,j) = a(5,i)
      else
        b(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine r85_vxm ( n, a, x, b )

!*****************************************************************************80
!
!! R85_VXM multiplies an R8VEC by an R85 matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, real ( kind = 8 ) A(5,N), the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) B(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(n)

  b(1:n)   =            a(3,1:n)   * x(1:n)
  b(2:n)   = b(2:n)   + a(4,1:n-1) * x(1:n-1)
  b(3:n)   = b(3:n)   + a(5,1:n-2) * x(1:n-2)
  b(1:n-1) = b(1:n-1) + a(2,2:n)   * x(2:n)
  b(1:n-2) = b(1:n-2) + a(1,3:n)   * x(3:n)

  return
end
subroutine r8bb_add ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8BB_ADD adds a value to an entry in an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
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
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real ( kind = 8 ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2),
!    the R8BB matrix.
!
!    Input, integer ( kind = 4 ) I, J, the row and column of the entry to be
!    incremented.  Some combinations of I and J are illegal.
!
!    Input, real ( kind = 8 ) VALUE, the value to be added to the 
!    (I,J)-th entry.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j
  real ( kind = 8 ) value

  if ( value == 0.0D+00 ) then
    return
  end if

  if ( i <= 0 .or. n1 + n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_ADD - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop
  end if

  if ( j <= 0 .or. n1 + n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_ADD - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ', j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
!  Normally, we would check the condition MU < (J-I), but the storage
!  format requires extra entries be set aside in case of pivoting, which
!  means that the condition becomes MU+ML < (J-I).
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( (mu+ml) < (j-i) .or. ml < (i-j) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8BB_ADD - Warning!'
      write ( *, '(a,i8,a,i8,a)' ) '  Unable to add to entry (', i, ',', j, ').'
    else
      ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
    end if
!
!  The A2 block of the matrix.
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = (2*ml+mu+1)*n1+(j-n1-1)*n1 + i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2 + (i-n1)
  end if

  a(ij) = a(ij) + value

  return
end
subroutine r8bb_fa ( n1, n2, ml, mu, a, pivot, info )

!*****************************************************************************80
!
!! R8BB_FA factors an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!    Once the matrix has been factored by R8BB_FA, R8BB_SL may be called
!    to solve linear systems involving the matrix.
!
!    R8BB_FA uses LINPACK routines to carry out the factorization.
!
!
!    The linear system must be border banded, of the form:
!
!      ( A1 A2 ) (X1) = (B1)
!      ( A3 A4 ) (X2)   (B2)
!
!    where A1 is a (usually big) banded square matrix, A2 and A3 are
!    column and row strips which may be nonzero, and A4 is a dense
!    square matrix.
!
!    The algorithm rewrites the system as:
!
!         X1 + inv(A1) A2 X2 = inv(A1) B1
!
!      A3 X1 +         A4 X2 = B2
!
!    and then rewrites the second equation as
!
!      ( A4 - A3 inv(A1) A2 ) X2 = B2 - A3 inv(A1) B1
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
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input/output, real ( kind = 8 ) A( (2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2 ).
!    On input, the border-banded matrix to be factored.
!    On output, information describing a partial factorization
!    of the original coefficient matrix.  This information is required
!    by R8BB_SL in order to solve linear systems associated with that
!    matrix.
!
!    Output, integer ( kind = 4 ) PIVOT(N1+N2), contains pivoting information.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n1+n2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nband

  nband = (2*ml+mu+1) * n1
!
!  Factor the A1 band matrix, overwriting A1 by its factors.
!
  if ( 0 < n1 ) then

    call r8gb_fa ( n1, ml, mu, a, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8BB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  R8GB_FA returned INFO = ', info
      write ( *, '(a)' ) '  Factoring failed for column INFO.'
      write ( *, '(a)' ) '  The band matrix A1 is singular.'
      write ( *, '(a)' ) '  This algorithm cannot continue!'
      stop
    end if

  end if

  if ( 0 < n1 .and. 0 < n2 ) then
!
!  Solve A1 * x = -A2 for x, and overwrite A2 by the results.
!
    do i = nband+1, nband+n1*n2
      a(i) = -a(i)
    end do

    job = 0
    do j = 1, n2
      call r8gb_sl ( n1, ml, mu, a, pivot, a(nband+(j-1)*n1+1), job )
    end do
!
!  A4 := A4 + A3 * A2.
!
    do i = 1, n2
      do j = 1, n1
        ij = nband + n1*n2 + (j-1)*n2 + i
        do k = 1, n2
          ik = nband + 2*n1*n2 + (k-1)*n2 + i
          jk = nband + (k-1) * n1 + j
          a(ik) = a(ik) + a(ij) * a(jk)
        end do
      end do
    end do

  end if
!
!  Factor A4.
!
  if ( 0 < n2 ) then

    call r8ge_fa ( n2, a(nband+2*n1*n2+1), pivot(n1+1), info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8BB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  R8GE_FA returned INFO = ',info
      write ( *, '(a)' ) '  This indicates singularity in column INFO.'
      write ( *, '(a,i8)' ) '  of the A4 submatrix, which is column ', n1+info
      write ( *, '(a)' ) '  of the full matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  It is possible that the full matrix is '
      write ( *, '(a)' ) '  nonsingular, but the algorithm R8BB_FA may'
      write ( *, '(a)' ) '  not be used for this matrix.'
      stop
    end if
  end if

  return
end
subroutine r8bb_get ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8BB_GET returns an entry of an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense
!    blocks. N1 and N2 must be nonnegative, and at least one must be positive.
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
!    Input, integer ( kind = 4 ) I, J, the row and column of the entry to 
!    be retrieved.
!
!    Output, real ( kind = 8 ) VALUE, the value of the (I,J) entry.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j
  real ( kind = 8 ) value

  if ( i <= 0 .or. n1+n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_GET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop
  end if

  if ( j <= 0 .or. n1+n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_GET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ', j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
!  Normally, we would check the condition MU < (J-I), but the storage
!  format requires extra entries be set aside in case of pivoting, which
!  means that the condition becomes MU+ML < (J-I).
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( mu+ml < (j-i) .or. ml < (i-j) ) then
      value = 0.0D+00
      return
    else
      ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
    end if
!
!  The A2 block of the matrix.
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
  end if

  value = a(ij)

  return
end
subroutine r8bb_indicator ( n1, n2, ml, mu, a )

!*****************************************************************************80
!
!! R8BB_INDICATOR sets up an R8BB indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!  Example:
!
!    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
!
!       00
!       00  00
!       00  00  00 --- ---
!      A11 A12 A13  00 ---  A16 A17
!      A21 A22 A23 A24  00  A26 A27
!      --- A32 A33 A34 A35  A36 A37
!      --- --- A43 A44 A45  A46 A47
!      --- --- --- A54 A55  A56 A57
!                       00
!
!      A61 A62 A63 A64 A65  A66 A67
!      A71 A72 A73 A74 A75  A76 A77
!
!    The matrix is actually stored as a vector, and we will simply suggest
!    the structure and values of the indicator matrix as:
!
!      00 00 00 00 00
!      00 00 13 24 35     16 17     61 62 63 64 65     66 67
!      00 12 23 34 45  +  26 27  +  71 72 73 74 75  +  76 77
!      11 22 33 44 55     36 37     
!      21 32 43 54 00     46 47     
!                         56 57     
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Output, real ( kind = 8 ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  fac = 10 ** ( i4_log_10 ( n1 + n2 ) + 1 )
!
!  Set the banded matrix A1.
!
  do j = 1, n1
    do row = 1, 2 * ml + mu + 1
      i = row + j - ml - mu - 1
      if ( ml < row .and. 1 <= i .and. i <= n1 ) then
        a(row+(j-1)*(2*ml+mu+1)) = real ( fac * i + j, kind = 8 )
      else
        a(row+(j-1)*(2*ml+mu+1)) = 0.0D+00
      end if
    end do
  end do
!
!  Set the N1 by N2 rectangular strip A2.
!
  base = ( 2 * ml + mu + 1 ) * n1

  do i = 1, n1
    do j = n1 + 1, n1 + n2
      a(base + i + (j-n1-1)*n1 ) = real ( fac * i + j, kind = 8 )
    end do
  end do
!
!  Set the N2 by N1 rectangular strip A3.
!
  base = ( 2 * ml + mu + 1 ) * n1 + n1 * n2

  do i = n1 + 1, n1 + n2
    do j = 1, n1    
      a(base + i-n1 + (j-1)*n2 ) = real ( fac * i + j, kind = 8 )
    end do
  end do
!
!  Set the N2 by N2 square A4.
!
  base = ( 2 * ml + mu + 1 ) * n1 + n1 * n2 + n2 * n1

  do i = n1 + 1, n1 + n2
    do j = n1 + 1, n1 + n2
      a(base + i-n1 + (j-n1-1)*n2 ) = real ( fac * i + j, kind = 8 )
    end do
  end do

  return
end
subroutine r8bb_mxv ( n1, n2, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8BB_MXV multiplies an R8BB matrix by an R8VEC.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense
!    blocks  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
!    Input, real ( kind = 8 ) X(N1+N2), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N1+N2), the result of multiplying A by X.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) b(n1+n2)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n1+n2)
!
!  Initialize B.
!
  b(1:n1+n2) = 0.0D+00
!
!  Multiply by A1.
!
  do j = 1, n1

    ilo = max ( 1, j - mu - ml )
    ihi = min ( n1, j + ml )
    ij = (j-1) * (2*ml+mu+1) - j + ml + mu + 1

    b(ilo:ihi) = b(ilo:ihi) + a(ij+ilo:ij+ihi) * x(j)

  end do
!
!  Multiply by A2.
!
  do j = n1+1, n1+n2
    ij = (2*ml+mu+1)*n1+(j-n1-1)*n1

    b(1:n1) = b(1:n1) + a(ij+1:ij+n1) * x(j)

  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1+n2
    ij = (2*ml+mu+1)*n1+n1*n2+(j-1)*n2-n1

    b(n1+1:n1+n2) = b(n1+1:n1+n2) + a(ij+n1+1:ij+n1+n2) * x(j)

  end do

  return
end
subroutine r8bb_print ( n1, n2, ml, mu, a, title )

!*****************************************************************************80
!
!! R8BB_PRINT prints an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  character ( len = * ) title

  call r8bb_print_some ( n1, n2, ml, mu, a, 1, 1, n1+n2, n1+n2, title )

  return
end
subroutine r8bb_print_some ( n1, n2, ml, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8BB_PRINT_SOME prints some of an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n1+n2 )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n1+n2 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0D+00

        if ( i <= n1 .and. j <= n1 ) then
          if ( (j-i) <= mu+ml .and. (i-j) <= ml ) then
            ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
            aij = a(ij)
          end if
        else if ( i <= n1 .and. n1 < j ) then
          ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
          aij = a(ij)
        else if ( n1 < i ) then
          ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
          aij = a(ij)
        end if

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8bb_random ( n1, n2, ml, mu, seed, a )

!*****************************************************************************80
!
!! R8BB_RANDOM randomizes an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  integer ( kind = 4 ) row
  integer ( kind = 4 ) seed
!
!  Randomize the banded matrix A1.
!  We still believe that the "junk" entries should be set to 0.
!
  do j = 1, n1
    do row = 1, 2*ml+mu+1
      i = row + j - ml - mu - 1
      if ( ml < row .and. 1 <= i .and. i <= n1 ) then
        r = r8_uniform_01 ( seed )
      else
        r = 0.0D+00
      end if
      a(row+(j-1)*(2*ml+mu+1)) = r
    end do
  end do
!
!  Randomize the rectangular strips A2+A3+A4.
!
  ilo = (2*ml+mu+1) * n1 + 1

  call r8vec_uniform_01 ( n1*n2+n2*n1+n2*n2, seed, a(ilo:) )

  return
end
subroutine r8bb_set ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8BB_SET sets an entry of an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real ( kind = 8 ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2),
!    the R8BB matrix.
!
!    Input, integer ( kind = 4 ) I, J, the row and column of the entry to 
!    be set.
!
!    Input, real ( kind = 8 ) VALUE, the value to be assigned to the
!    (I,J) entry.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j
  real ( kind = 8 ) value

  if ( i <= 0 .or. n1+n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop
  end if

  if ( j <= 0 .or. n1+n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ', j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
!  Normally, we would check the condition MU < (J-I), but the storage
!  format requires extra entries be set aside in case of pivoting, which
!  means that the condition becomes MU+ML < (J-I).
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( mu+ml < (j-i) .or. ml < (i-j) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8BB_SET - Warning!'
      write ( *, '(a,i8,a,i8,a)' ) '  Unable to set entry (', i, ',', j, ').'
      return
    else
      ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
    end if
!
!  The A2 block of the matrix.
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
  end if

  a(ij) = value

  return
end
subroutine r8bb_sl ( n1, n2, ml, mu, a_lu, pivot, b )

!*****************************************************************************80
!
!! R8BB_SL solves an R8BB system factored by R8BB_FA.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!    The linear system A * x = b is decomposable into the block system:
!
!      ( A1 A2 ) * (X1) = (B1)
!      ( A3 A4 )   (X2)   (B2)
!
!    All the arguments except B are input quantities only, which are
!    not changed by the routine.  They should have exactly the same values
!    they had on exit from R8BB_FA.
!
!    If more than one right hand side is to be solved, with the same matrix,
!    R8BB_SL should be called repeatedly.  However, R8BB_FA only needs to be
!    called once to create the factorization.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A_LU( (2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2), 
!    the LU factors from R8BB_FA.
!
!    Input, integer ( kind = 4 ) PIVOT(N1+N2), the pivoting information 
!    from R8BB_FA.
!
!    Input/output, real ( kind = 8 ) B(N1+N2).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a_lu((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) b(n1+n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) pivot(n1+n2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) nband

  nband = (2*ml+mu+1)*n1
!
!  Set B1 := inverse(A1) * B1.
!
  if ( 0 < n1 ) then
    job = 0
    call r8gb_sl ( n1, ml, mu, a_lu, pivot, b, job )
  end if
!
!  Modify the right hand side of the second linear subsystem.
!  Set B2 := B2 - A3*B1.
!
  do i = 1, n2
    do j = 1, n1
      ij = nband + n1*n2 + (j-1)*n2 + i
      b(n1+i) = b(n1+i) - a_lu(ij) * b(j)
    end do
  end do
!
!  Set B2 := inverse(A4) * B2.
!
  if ( 0 < n2 ) then
    job = 0
    call r8ge_sl ( n2, a_lu(nband+2*n1*n2+1), pivot(n1+1), b(n1+1), job )
  end if
!
!  Modify the first subsolution.
!  Set B1 := B1 + A2*B2.
!
  do i = 1, n1
    do j = 1, n2
      ij = nband + (j-1)*n1 + i
      b(i) = b(i) + a_lu(ij) * b(n1+j)
    end do
  end do

  return
end
subroutine r8bb_to_r8ge ( n1, n2, ml, mu, a, b )

!*****************************************************************************80
!
!! R8BB_TO_R8GE copies an R8BB matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
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
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
!    Output, real ( kind = 8 ) B(N1+N2,N1+N2), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) b(n1+n2,n1+n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j

  do i = 1, n1
    do j = 1, n1

      if ( mu+ml < (j-i) .or. ml < (i-j) ) then
        b(i,j) = 0.0D+00
      else
        ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
        b(i,j) = a(ij)
      end if

    end do
  end do

  do i = 1, n1
    do j = n1+1, n2
      ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
      b(i,j) = a(ij)
    end do
  end do

  do i = n1+1, n2
    do j = 1, n1+n2
      ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
      b(i,j) = a(ij)
    end do
  end do

  return
end
subroutine r8bb_vxm ( n1, n2, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8BB_VXM multiplies an R8VEC by an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A((2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8BB matrix.
!
!    Input, real ( kind = 8 ) X(N1+N2), the vector to multiply A.
!
!    Output, real ( kind = 8 ) B(N1+N2), the product X times A.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) b(n1+n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n1+n2)
!
!  Initialize B.
!
  b(1:n1+n2) = 0.0D+00
!
!  Multiply by A1.
!
  do j = 1, n1
    ilo = max ( 1, j - mu - ml )
    ihi = min ( n1, j + ml )
    ij = (j-1) * (2*ml+mu+1) - j + ml + mu + 1
    do i = ilo, ihi
      b(j) = b(j) + x(i) * a(ij+i)
    end do
  end do
!
!  Multiply by A2.
!
  do j = n1+1, n1+n2
    ij = (2*ml+mu+1)*n1+(j-n1-1)*n1
    do i = 1, n1
      b(j) = b(j) + x(i) * a(ij+i)
    end do
  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1+n2
    ij = (2*ml+mu+1)*n1+n1*n2+(j-1)*n2-n1
    do i = n1+1, n1+n2
      b(j) = b(j) + x(i) * a(ij+i)
    end do
  end do

  return
end
subroutine r8blt_det ( n, ml, a, det )

!*****************************************************************************80
!
!! R8BLT_DET computes the determinant of an R8BLT matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth.
!
!    Input, real ( kind = 8 ) A(ML+1,N), the R8BLT matrix.
!
!    Output, real ( kind = 8 ) DET, the determinant of A.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) det

  det = product ( a(1,1:n) )

  return
end
subroutine r8blt_indicator ( n, ml, a )

!*****************************************************************************80
!
!! R8BLT_INDICATOR sets up an R8BLT indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!    The indicator matrix is stored as:
!
!      11  22  33  44  55
!      21  32  43  54   0
!      31  42  53   0   0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth.
!
!    Output, real ( kind = 8 ) A(ML+1,N), the R8BLT matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, n
    do j = max ( 1, i - ml ), i
      a(i-j+1,j) = real ( fac * i + j, kind = 8 )
    end do
  end do

  do i = n+1, n+ml
    do j = i-ml, n
      a(i-j+1,j) = 0.0D+00
    end do
  end do

  return
end
subroutine r8blt_mxv ( n, ml, a, x, b )

!*****************************************************************************80
!
!! R8BLT_MXV multiplies an R8BLT matrix by an R8VEC.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth.
!
!    Input, real ( kind = 8 ) A(ML+1,N), the R8BLT matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  real ( kind = 8 ) x(n)

  do i = 1, n
    b(i) = 0.0D+00
    jlo = max ( 1, i - ml )
    jhi = i
    do j = jlo, jhi
      b(i) = b(i) + a(i-j+1,j) * x(j)
    end do
  end do

  return
end
subroutine r8blt_print ( n, ml, a, title )

!*****************************************************************************80
!
!! R8BLT_PRINT prints an R8BLT matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth.
!
!    Input, real ( kind = 8 ) A(ML+1,N), the R8BLT matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  character ( len = * ) title

  call r8blt_print_some ( n, ml, a, 1, 1, n, n, title )

  return
end
subroutine r8blt_print_some ( n, ml, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8BLT_PRINT_SOME prints some of an R8BLT matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth.
!
!    Input, real ( kind = 8 ) A(ML+1,N), the R8BLT matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo )
    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + ml )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j <= i .and. i <= j + ml ) then
          aij = a(i-j+1,j)
        else
          aij = 0.0D+00
        end if

        if ( ml < i-j .or. 0 < j-i ) then
          ctemp(j2) = '              '
        else if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8blt_random ( n, ml, seed, a )

!*****************************************************************************80
!
!! R8BLT_RANDOM randomizes an R8BLT matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) A(ML+1,N), the R8BLT matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  do i = 1, n
    do j = max ( 1, i - ml ), i
      a(i-j+1,j) = r8_uniform_01 ( seed )
    end do
  end do
!
!  The junk entries can be thought of as corresponding to
!  elements of a phantom portion of the matrix.
!
  do i = n+1, n+ml
    do j = i-ml, n
      a(i-j+1,j) = 0.0D+00
    end do
  end do

  return
end
subroutine r8blt_sl ( n, ml, a, b, job )

!*****************************************************************************80
!
!! R8BLT_SL solves an R8BLT system.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!    No factorization of the lower triangular matrix is required.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth.
!
!    Input, real ( kind = 8 ) A(ML+1,N), the R8BLT matrix.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB, is 0 to solve the untransposed system,
!    nonzero to solve the transposed system.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job

  if ( job == 0 ) then

    do j = 1, n
      b(j) = b(j) / a(1,j)
      ihi = min ( j + ml, n )
      do i = j+1, ihi
        b(i) = b(i) - a(i-j+1,j) * b(j)
      end do
    end do

  else

    do j = n, 1, -1
      b(j) = b(j) / a(1,j)
      ilo = max ( j - ml, 1 )
      do i = ilo, j-1
        b(i) = b(i) - a(j-i+1,i) * b(j)
      end do
    end do

  end if

  return
end
subroutine r8blt_to_r8ge ( n, ml, a, b )

!*****************************************************************************80
!
!! R8BLT_TO_R8GE copies an R8BLT matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth of A.
!    ML must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(ML+1,N), the R8BLT matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n
      if ( j <= i .and. i <= j + ml ) then
        b(i,j) = a(i-j+1,j)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8blt_vxm ( n, ml, a, x, b )

!*****************************************************************************80
!
!! R8BLT_VXM multiplies an R8VEC by an R8BLT matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth.
!
!    Input, real ( kind = 8 ) A(ML+1,N), the R8BLT matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product X*A.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  do i = 1, n
    jlo = max ( 1, i - ml )
    jhi = i
    do j = jlo, jhi
      b(j) = b(j) + x(i) * a(i-j+1,j)
    end do
  end do

  return
end
subroutine r8bto_indicator ( m, l, a )

!*****************************************************************************80
!
!! R8BTO_INDICATOR sets up an R8BTO indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Example:
!
!    M = 2, L = 3
!
!    1 2 | 3 4 | 5 6
!    5 5 | 6 6 | 7 7
!    ----+-----+-----
!    7 8 | 1 2 | 3 4
!    8 8 | 5 5 | 6 6
!    ----+-----+-----
!    9 0 | 7 8 | 1 2
!    9 9 | 8 8 | 5 5
!
!    X = (/ 1, 2, 3, 4, 5, 6 /)
!
!    B = (/ 91, 134, 73, 125, 97, 129 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Output, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k

  fac = 10 ** ( i4_log_10 ( m * l ) + 1 )
!
!  Blocks 1 to L form the first row.
!
  j = 0

  do k = 1, l

    do j2 = 1, m
      j = j + 1
      do i = 1, m
        a(i,j2,k) = real ( fac * i + j, kind = 8 )
      end do
    end do
  end do
!
!  Blocks L+1 through 2*L-1 form the remainder of the first column.
!
  i = m

  do k = l+1, 2*l-1

    do i2 = 1, m
      i = i + 1
      do j = 1, m
        a(i2,j,k) = real ( fac * i + j, kind = 8 )
      end do
    end do

  end do

  return
end
subroutine r8bto_mxv ( m, l, a, x, b )

!*****************************************************************************80
!
!! R8BTO_MXV multiplies an R8BTO matrix by an R8VEC.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Example:
!
!    M = 2, L = 3
!
!    1 2 | 3 4 | 5 6
!    5 5 | 6 6 | 7 7
!    ----+-----+-----
!    7 8 | 1 2 | 3 4
!    8 8 | 5 5 | 6 6
!    ----+-----+-----
!    9 0 | 7 8 | 1 2
!    9 9 | 8 8 | 5 5
!
!    X = (/ 1, 2, 3, 4, 5, 6 /)
!
!    B = (/ 91, 134, 73, 125, 79, 138 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Input, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
!    Input, real ( kind = 8 ) X(M*L), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(M*L), the product A * X.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  real ( kind = 8 ) b(m,l)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,l)
!
!  Construct the right hand side by blocks.
!
  do i = 1, l

    b(1:m,i) = 0.0D+00

    do j = 1, i-1
      b(1:m,i) = b(1:m,i) + matmul ( a(1:m,1:m,l+i-j), x(1:m,j) )
    end do

    do j = i, l
      b(1:m,i) = b(1:m,i) + matmul ( a(1:m,1:m,j+1-i), x(1:m,j) )
    end do

  end do

  return
end
subroutine r8bto_print ( m, l, a, title )

!*****************************************************************************80
!
!! R8BTO_PRINT prints an R8BTO matrix.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Input, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  character ( len = * ) title

  call r8bto_print_some ( m, l, a, 1, 1, m*l, m*l, title )

  return
end
subroutine r8bto_print_some ( m, l, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8BTO_PRINT_SOME prints some of an R8BTO matrix.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Input, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3hi
  integer ( kind = 4 ) i3lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) j3hi
  integer ( kind = 4 ) j3lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) n
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  n = m * l
!
!  Print the columns of the matrix, in strips of 5.
!
  do j3lo = jlo, jhi, incx

    j3hi = j3lo + incx - 1
    j3hi = min ( j3hi, n )
    j3hi = min ( j3hi, jhi )

    inc = j3hi + 1 - j3lo

    write ( *, '(a)' ) ' '

    do j = j3lo, j3hi
      j3 = j + 1 - j3lo
      write ( ctemp(j3), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j3), j3 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i3lo = max ( ilo, 1 )
    i3hi = min ( ihi, n )

    do i = i3lo, i3hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j3 = 1, inc

        j = j3lo - 1 + j3
!
!  i = M * ( i1 - 1 ) + i2
!  j = M * ( j1 - 1 ) + j2
!
        i1 = ( i - 1 ) / m + 1
        i2 = i - m * ( i1 - 1 )
        j1 = ( j - 1 ) / m + 1
        j2 = j - m * ( j1 - 1 )

        if ( i1 <= j1 ) then
          aij = a(i2,j2,j1+1-i1)
        else
          aij = a(i2,j2,l+i1-j1)
        end if

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j3), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j3), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j3), j3 = 1, inc )

    end do

  end do

  return
end
subroutine r8bto_random ( m, l, seed, a )

!*****************************************************************************80
!
!! R8BTO_RANDOM randomizes an R8BTO matrix.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A1(M,M,2*L-1), the R8BTO matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) l

  real ( kind = 8 ) a(m,m,2*l-1)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  do i = 1, m
    do j = 1, m
      do k = 1, 2 * l - 1 
        a(i,j,k) = r8_uniform_01 ( seed )
      end do
    end do
  end do

  return
end
subroutine r8bto_sl ( m, l, a, b, x )

!*****************************************************************************80
!
!! R8BTO_SL solves an R8BTO system.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Example:
!
!    M = 2, L = 3
!
!    1 2 | 3 4 | 5 6
!    5 5 | 6 6 | 7 7
!    ----+-----+-----
!    7 8 | 1 2 | 3 4
!    8 8 | 5 5 | 6 6
!    ----+-----+-----
!    9 0 | 7 8 | 1 2
!    9 9 | 8 8 | 5 5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column.
!
!    Input, real ( kind = 8 ) A(M*M,2*L-1), the R8BTO matrix.
!
!    Input, real ( kind = 8 ) B(M*L), the right hand side vector.
!
!    Output, real ( kind = 8 ) X(M*L), the solution vector.  X and B
!    may share storage.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m*m,2*l-1)
  real ( kind = 8 ) b(m,l)
  real ( kind = 8 ) c1(m,m,l-1)
  real ( kind = 8 ) c2(m,m,l-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) n
  integer ( kind = 4 ) pivot(m)
  real ( kind = 8 ) r1(m,m)
  real ( kind = 8 ) r2(m,m)
  real ( kind = 8 ) r3(m,m)
  real ( kind = 8 ) r5(m,m)
  real ( kind = 8 ) r6(m,m)
  real ( kind = 8 ) x(m,l)
!
!  Solve the system with the principal minor of order M.
!
  i3 = 1
  do j = 1, m
    do i = 1, m
      c1(i,j,1) = a(i3,1)
      r1(i,j) = a(i3,1)
      i3 = i3 + 1
    end do
  end do

  r3(1:m,1:m) = r1(1:m,1:m)
  x(1:m,1) = b(1:m,1)

  call r8ge_fa ( m, r3, pivot, info )

  job = 0
  call r8ge_sl ( m, r3, pivot, x(1,1), job )

  if ( l == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system
!  with the block Toeplitz matrix for N = 2 through L.
!
  do n = 2, l
!
!  Compute multiples of the first and last block columns of
!  the inverse of the principal minor of order M*N.
!
    i3 = 1
    do j = 1, m
      do i = 1, m
        r5(i,j) = a(i3,l+n-1)
        r6(i,j) = a(i3,n)
        i3 = i3 + 1
      end do
    end do

    if ( 2 < n ) then

      c1(1:m,1:m,n-1) = r2(1:m,1:m)

      do i1 = 1, n-2
        i2 = n - i1
        do j = 1, m
          i3 = 1
!
!  My apologies, but we have an I1+L, followed in the next line by I1+1.
!
          do i = 1, m
            call daxpy ( m, c1(i,j,i2), a(i3,i1+l), 1, r5(1,j), 1 )
            call daxpy ( m, c2(i,j,i1), a(i3,i1+1), 1, r6(1,j), 1 )
            i3 = i3 + m
          end do
        end do
      end do

    end if

    do j = 1, m
      r2(1:m,j) = -r5(1:m,j)
      job = 0
      call r8ge_sl ( m, r3, pivot, r2(1,j), job )
    end do

    r3(1:m,1:m) = r6(1:m,1:m)
    r6(1:m,1:m) = -c1(1:m,1:m,1)

    do j = 1, m
      do i = 1, m
        c1(1:m,j,1) = c1(1:m,j,1) + r2(i,j) * r3(1:m,i)
      end do
    end do

    call r8ge_fa ( m, r6, pivot, info )

    do j = 1, m
      call r8ge_sl ( m, r6, pivot, r3(1,j), job )
      do i = 1, m
        r1(1:m,j) = r1(1:m,j) + r3(i,j) * r5(1:m,i)
      end do
    end do

    if ( 2 < n ) then

      r6(1:m,1:m) = c2(1:m,1:m,1)

      do i1 = 2, n-1

        if ( i1 /= n-1 ) then
          r5(1:m,1:m) = c2(1:m,1:m,i1)
        end if

        do j = 1, m
          c2(1:m,j,i1) = r6(1:m,j)
          do i = 1, m
            call daxpy ( m, r3(i,j), c1(1,i,i1), 1, c2(1,j,i1), 1 )
          end do
        end do

        do j = 1, m
          do i = 1, m
            call daxpy ( m, r2(i,j), r6(1,i), 1, c1(1,j,i1), 1 )
          end do
        end do

        r6(1:m,1:m) = r5(1:m,1:m)

      end do

    end if

    c2(1:m,1:m,1) = r3(1:m,1:m)
!
!  Compute the solution of the system with the principal minor of order M*N.
!
    r3(1:m,1:m) = r1(1:m,1:m)
    x(1:m,n) = b(1:m,n)

    do i1 = 1, n-1
      i2 = n - i1
      i3 = 1
      do i = 1, m
        call daxpy ( m, -x(i,i2), a(i3,i1+l), 1, x(1,n), 1 )
        i3 = i3 + m
      end do
    end do

    call r8ge_fa ( m, r3, pivot, info )

    call r8ge_sl ( m, r3, pivot, x(1,n), job )

    do i1 = 1, n-1
      do i = 1, m
        x(1:m,i1) = x(1:m,i1) + x(i,n) * c2(1:m,i,i1)
!       call daxpy ( m, x(i,n), c2(1,i,i1), 1, x(1,i1), 1 )
      end do
    end do

  end do

  return
end
subroutine r8bto_to_r8ge ( m, l, a, n, b )

!*****************************************************************************80
!
!! R8BTO_TO_R8GE copies an R8BTO matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the R8BTO matrix.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or column
!    of the R8BTO matrix.
!
!    Input, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
!    Output, integer ( kind = 4 ) N, the order of the matrix, N = M*L.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  real ( kind = 8 ) b(m*l,m*l)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n

  n = m * l

  do i = 1, n

    i1 = ( i - 1 ) / m + 1
    i2 = i - m * ( i1 - 1 )

    do j = 1, n

      j1 = ( j - 1 ) / m + 1
      j2 = j - m * ( j1 - 1 )

      if ( i1 <= j1 ) then
        b(i,j) = a(i2,j2,j1+1-i1)
      else
        b(i,j) = a(i2,j2,l+i1-j1)
      end if

    end do

  end do

  return
end
subroutine r8bto_vxm ( m, l, a, x, b )

!*****************************************************************************80
!
!! R8BTO_VXM multiplies an R8VEC by an R8BTO matrix.
!
!  Discussion:
!
!    The R8BTO storage format is for a block Toeplitz matrix. The matrix
!    can be regarded as an L by L array of blocks, each of size M by M.
!    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
!    that is, along its diagonal, the blocks repeat.
!
!    Storage for the matrix consists of the L blocks of the first row,
!    followed by the L-1 blocks of the first column (skipping the first row).
!    These items are stored in the natural way in an (M,M,2*L-1) array.
!
!  Example:
!
!    M = 2, L = 3
!
!    1 2 | 3 4 | 5 6
!    5 5 | 6 6 | 7 7
!    ----+-----+-----
!    7 8 | 1 2 | 3 4
!    8 8 | 5 5 | 6 6
!    ----+-----+-----
!    9 0 | 7 8 | 1 2
!    9 9 | 8 8 | 5 5
!
!    X = (/ 1, 2, 3, 4, 5, 6 /)
!
!    B = (/ 163, 122, 121, 130, 87, 96 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the blocks of the matrix A.
!
!    Input, integer ( kind = 4 ) L, the number of blocks in a row or 
!    column of A.
!
!    Input, real ( kind = 8 ) A(M,M,2*L-1), the R8BTO matrix.
!
!    Input, real ( kind = 8 ) X(M*L), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(M*L), the product X * A.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m,2*l-1)
  real ( kind = 8 ) b(m,l)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m,l)
!
!  Construct the right hand side by blocks.
!
  do i = 1, l

    b(1:m,i) = 0.0D+00

    do j = 1, i
      b(1:m,i) = b(1:m,i) + matmul ( transpose ( a(1:m,1:m,i+1-j) ), x(1:m,j) )
    end do

    do j = i+1, l
      b(1:m,i) = b(1:m,i) + matmul ( transpose ( a(1:m,1:m,l+j-i) ), x(1:m,j) )
    end do

  end do

  return
end
subroutine r8but_det ( n, mu, a, det )

!*****************************************************************************80
!
!! R8BUT_DET computes the determinant of an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8BUT matrix.
!
!    Output, real ( kind = 8 ) DET, the determinant of A.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) det

  det = product ( a(mu+1,1:n) )

  return
end
subroutine r8but_indicator ( n, mu, a )

!*****************************************************************************80
!
!! R8BUT_INDICATOR sets up an R8BUT indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!    The indicator matrix is stored as:
!
!       0   0  13  24  35
!       0  12  23  34  45
!      11  22  33  44  55
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth.
!
!    Output, real ( kind = 8 ) A(MU+1,N), the R8BUT matrix.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, n
    do j = i, min ( n, i + mu )
      a(i-j+mu+1,j) = real ( fac * i + j, kind = 8 )
    end do
  end do

  do i = 1, mu
    do j = 1, mu+1-i
      a(i,j) = 0.0D+00
    end do
  end do

  return
end
subroutine r8but_mxv ( n, mu, a, x, b )

!*****************************************************************************80
!
!! R8BUT_MXV multiplies an R8BUT matrix by an R8VEC.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8BUT matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  do i = 1, n
    do j = i, min ( n, i + mu )
      b(i) = b(i) + a(i-j+mu+1,j) * x(j)
    end do
  end do

  return
end
subroutine r8but_print ( n, mu, a, title )

!*****************************************************************************80
!
!! R8BUT_PRINT prints an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8BUT matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  character ( len = * ) title

  call r8but_print_some ( n, mu, a, 1, 1, n, n, title )

  return
end
subroutine r8but_print_some ( n, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8BUT_PRINT_SOME prints some of an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8BUT matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo )
    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + mu )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j .and. j <= i + mu ) then
          aij = a(i-j+mu+1,j)
          if ( r8_is_int ( aij ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) aij
          else
            write ( ctemp(j2), '(g14.6)' ) aij
          end if
        else
          ctemp(j2) = '              '
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8but_random ( n, mu, seed, a )

!*****************************************************************************80
!
!! R8BUT_RANDOM randomizes an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(MU+1,N), the R8BUT matrix.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  do i = 1, mu + 1

    do j = 1, mu + 1 - i
      a(i,j) = 0.0D+00
    end do

    do j = max ( 1, mu + 2 - i ), n
      a(i,j) = r8_uniform_01 ( seed )
    end do

  end do

  return
end
subroutine r8but_sl ( n, mu, a, b, job )

!*****************************************************************************80
!
!! R8BUT_SL solves an R8BUT system.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8BUT matrix.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB, is 0 to solve the untransposed system,
!    nonzero to solve the transposed system.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) job

  if ( job == 0 ) then

    do j = n, 1, -1
      b(j) = b(j) / a(j-j+mu+1,j)
      jlo = max ( 1, j - mu )
      b(jlo:j-1) = b(jlo:j-1) - a(jlo-j+mu+1:j-1-j+mu+1,j) * b(j)
    end do

  else

    do j = 1, n
      b(j) = b(j) / a(j-j+mu+1,j)
      ihi = min ( n, j + mu )
      do i = j + 1, ihi
        b(i) = b(i) - a(j-i+mu+1,i) * b(j)
      end do
    end do

  end if

  return
end
subroutine r8but_to_r8ge ( n, mu, a, b )

!*****************************************************************************80
!
!! R8BUT_TO_R8GE copies an R8BUT matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth of A.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8BUT matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n
      if ( i <= j .and. j <= i+mu ) then
        b(i,j) = a(i-j+mu+1,j)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8but_vxm ( n, mu, a, x, b )

!*****************************************************************************80
!
!! R8BUT_VXM multiplies an R8VECr by an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8BUT matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product X*A.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilo
  real ( kind = 8 ) x(n)

  do i = 1, n
    ilo = max ( 1, i - mu )
    b(i) = sum ( x(ilo:i) * a(ilo-i+mu+1:mu+1,i) )
  end do

  return
end
subroutine r8cb_det ( n, ml, mu, a_lu, det )

!*****************************************************************************80
!
!! R8CB_DET computes the determinant of an R8CB matrix factored by R8CB_NP_FA.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 March 2004
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
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A_LU(ML+MU+1,N), the LU factors from R8CB_NP_FA.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(ml+mu+1,n)
  real ( kind = 8 ) det

  det = product ( a_lu(mu+1,1:n) )

  return
end
subroutine r8cb_indicator ( m, n, ml, mu, a )

!*****************************************************************************80
!
!! R8CB_INDICATOR sets up an R8CB indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically ML+MU+1 by N.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Output, real ( kind = 8 ) A(ML+MU+1,N), the R8CB matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+mu+1,n)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) value

  fac = 10 ** ( i4_log_10 ( n ) + 1 )
  k = 0

  do j = 1, n
    do diag = 1, ml + mu + 1

      i = diag + j - mu - 1

      if ( 1 <= i .and. i <= m .and. i - ml <= j .and. j <= i + mu ) then
        value = real ( fac * i + j, kind = 8 )
      else
        k = k + 1
        value = - real ( k, kind = 8 )
      end if

      a(diag,j) = value

    end do
  end do


  return
end
subroutine r8cb_ml ( n, ml, mu, a_lu, x, b, job )

!*****************************************************************************80
!
!! R8CB_ML computes A * x or A' * X, using R8CB_NP_FA factors.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!    It is assumed that R8CB_NP_FA has overwritten the original matrix
!    information by LU factors.  R8CB_ML is able to reconstruct the
!    original matrix from the LU factor data.
!
!    R8CB_ML allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 October 2003
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
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A_LU(ML+MU+1,N), the LU factors from R8CB_NP_FA.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(N), the result of the multiplication.
!
!    Input, integer ( kind = 4 ) JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * x.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) job
  real ( kind = 8 ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      ilo = max ( 1, j - mu )
      do i = ilo, j - 1
        b(i) = b(i) + a_lu(i-j+mu+1,j) * b(j)
      end do
      b(j) = a_lu(j-j+mu+1,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n-1, 1, -1
      ihi = min ( n, j + ml )
      b(j+1:ihi) = b(j+1:ihi) - a_lu(mu+2:ihi-j+mu+1,j) * b(j)
    end do

  else
!
!  Y = ( PL )' * X.
!
    do j = 1, n-1

      ihi = min ( n, j + ml )
      do i = j+1, ihi
        b(j) = b(j) - b(i) * a_lu(i-j+mu+1,j)
      end do

    end do
!
!  B = U' * Y = ( PL * U )' * X = A' * X.
!
    do i = n, 1, -1
      jhi = min ( n, i + mu )
      do j = i+1, jhi
        b(j) = b(j) + b(i) * a_lu(i-j+mu+1,j)
      end do
      b(i) = b(i) * a_lu(i-i+mu+1,i)
    end do

  end if

  return
end
subroutine r8cb_mxv ( n, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8CB_MXV multiplies an R8CB matrix by an R8VEC.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 1998
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
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(ML+MU+1,N), the R8CB matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  real ( kind = 8 ) x(n)

  do i = 1, n
    b(i) = 0.0D+00
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      b(i) = b(i) + a(i-j+mu+1,j) * x(j)
    end do
  end do

  return
end
subroutine r8cb_np_fa ( n, ml, mu, a, info )

!*****************************************************************************80
!
!! R8CB_NP_FA factors an R8CB matrix by Gaussian elimination.
!
!  Discussion:
!
!    The R8CB storage format is appropriate for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!    R8CB_NP_FA is a version of the LINPACK routine R8GBFA, modifed to use
!    no pivoting, and to be applied to the R8CB compressed band matrix storage
!    format.  It will fail if the matrix is singular, or if any zero
!    pivot is encountered.
!
!    If R8CB_NP_FA successfully factors the matrix, R8CB_NP_SL may be called
!    to solve linear systems involving the matrix.
!
!    The matrix is stored in a compact version of LINPACK general
!    band storage, which does not include the fill-in entires.
!    The following program segment will store the entries of a banded
!    matrix in the compact format used by this routine:
!
!      m = mu+1
!      do j = 1, n
!        i1 = max ( 1, j-mu )
!        i2 = min ( n, j+ml )
!        do i = i1, i2
!          k = i-j+m
!          a(k,j) = afull(i,j)
!        end do
!      end do
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
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input/output, real ( kind = 8 ) A(ML+MU+1,N), the compact band matrix.
!    On input, the coefficient matrix of the linear system.
!    On output, the LU factors of the matrix.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+mu+1,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
!
!  The value of M is MU + 1 rather than ML + MU + 1.
!
  m = mu + 1
  info = 0
  ju = 0

  do k = 1, n-1
!
!  If our pivot entry A(MU+1,K) is zero, then we must give up.
!
    if ( a(m,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CB_NP_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  LM counts the number of nonzero elements that lie below the current
!  diagonal entry, A(K,K).
!
!  Multiply the LM entries below the diagonal by -1/A(K,K), turning
!  them into the appropriate "multiplier" terms in the L matrix.
!
    lm = min ( ml, n-k )
    a(m+1:m+lm,k) = - a(m+1:m+lm,k) / a(m,k)
!
!  MM points to the row in which the next entry of the K-th row is, A(K,J).
!  We then add L(I,K)*A(K,J) to A(I,J) for rows I = K+1 to K+LM.
!
    ju = max ( ju, mu + k )
    ju = min ( ju, n )
    mm = m

    do j = k+1, ju
      mm = mm - 1
      a(mm+1:mm+lm,j) = a(mm+1:mm+lm,j) + a(mm,j) * a(m+1:m+lm,k)
    end do

  end do

  if ( a(m,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CB_NP_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop
  end if

  return
end
subroutine r8cb_np_sl ( n, ml, mu, a_lu, b, job )

!*****************************************************************************80
!
!! R8CB_NP_SL solves an R8CB system factored by R8CB_NP_FA.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!    R8CB_NP_SL can also solve the related system A' * x = b.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 1998
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
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A_LU(ML+MU+1,N), the LU factors from R8CB_NP_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B contains the right hand side of the linear system, B.
!    On output, B contains the solution of the linear system, X.
!
!    Input, integer ( kind = 4 ) JOB.
!    If JOB is zero, the routine will solve A * x = b.
!    If JOB is nonzero, the routine will solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
!
!  The value of M is ML + 1, rather than MU + ML + 1.
!
  m = mu + 1
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve PL * Y = B.
!
    if ( 0 < ml ) then
      do k = 1, n-1
        lm = min ( ml, n-k )
        b(k+1:k+lm) = b(k+1:k+lm) + b(k) * a_lu(m+1:m+lm,k)
      end do
    end if
!
!  Solve U * X = Y.
!
    do k = n, 1, -1

      b(k) = b(k) / a_lu(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm

      b(lb:lb+lm-1) = b(lb:lb+lm-1) - b(k) * a_lu(la:la+lm-1,k)

    end do
!
!  Solve A' * X = B.
!
  else
!
!  Solve U' * Y = B.
!
    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm

      b(k) = ( b(k) - sum ( a_lu(la:la+lm-1,k) * b(lb:lb+lm-1) ) ) &
        / a_lu(m,k)

    end do
!
!  Solve ( PL )' * X = Y.
!
    if ( 0 < ml ) then

      do k = n-1, 1, -1
        lm = min ( ml, n-k )
        b(k) = b(k) + sum ( a_lu(m+1:m+lm,k) * b(k+1:k+lm) )
      end do

    end if

  end if

  return
end
subroutine r8cb_print ( m, n, ml, mu, a, title )

!*****************************************************************************80
!
!! R8CB_PRINT prints an R8CB matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = 8 ) A(ML+MU+1,N), the R8CB matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+mu+1,n)
  integer ( kind = 4 ) m
  character ( len = * ) title

  call r8cb_print_some ( m, n, ml, mu, a, 1, 1, m, n, title )

  return
end
subroutine r8cb_print_some ( m, n, ml, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8CB_PRINT_SOME prints some of an R8CB matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = 8 ) A(ML+MU+1,N), the R8CB matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+mu+1,n)
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) m
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - mu )
    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi + ml )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( ml < i - j .or. mu < j - i ) then
          ctemp(j2) = '              '
        else if ( r8_is_int ( a(i-j+mu+1,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i-j+mu+1,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+mu+1,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8cb_random ( n, ml, mu, seed, a )

!*****************************************************************************80
!
!! R8CB_RANDOM randomizes an R8CB matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2003
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
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) A(ML+MU+1,N), the R8CB matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+mu+1,n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
!
!  Set the entries that correspond to matrix elements.
!
  do j = 1, n

    ilo = max ( 1, j - mu )
    ihi = min ( n, j + ml )

    do i = j - mu, 0
      a(i-j+mu+1,j) = 0.0D+00
    end do

    do i = ilo, ihi
      a(i-j+mu+1,j) = r8_uniform_01 ( seed )
    end do

    do i = n + 1, j + ml
      a(i-j+mu+1,j) = 0.0D+00 
    end do

  end do

  return
end
subroutine r8cb_to_r8vec ( m, n, ml, mu, a, x )

!*****************************************************************************80
!
!! R8CB_TO_R8VEC copies an R8CB matrix to an R8VEC.
!
!  Discussion:
!
!    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
!    a data item carries its dimensionality implicitly, and so cannot be
!    regarded sometimes as a vector and sometimes as an array.
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the array.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!
!    Input, real ( kind = 8 ) A(ML+MU+1,N), the array to be copied.
!
!    Output, real ( kind = 8 ) X((ML+MU+1)*N), the vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+mu+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  real ( kind = 8 ) x((ml+mu+1)*n)

  do j = 1, n

    ihi = min ( mu, mu + 1 - j )
    do i = 1, ihi
      x(i+(j-1)*(ml+mu+1)) = 0.0D+00
    end do

    ilo = max ( ihi + 1, 1 )
    ihi = min ( ml+mu+1, mu+1+m-j )
    do i = ilo, ihi
      x(i+(j-1)*(ml+mu+1)) = a(i,j)
    end do

    ilo = ihi + 1
    ihi = ml+mu+1
    do i = ilo, ihi
      x(i+(j-1)*(ml+mu+1)) = 0.0D+00
    end do

  end do

  return
end
subroutine r8cb_to_r8ge ( n, ml, mu, a, b )

!*****************************************************************************80
!
!! R8CB_TO_R8GE copies an R8CB matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths of A.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(ML+MU+1,N), the R8CB matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+mu+1,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n
      if ( j-mu <= i .and. i <= j+ml ) then
        b(i,j) = a(mu+1+i-j,j)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8cb_vxm ( n, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8CB_VXM multiplies an R8VECr by an R8CB matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 October 1998
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
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(ML+MU+1,N), the R8CB matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product X*A.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  do i = 1, n
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      b(j) = b(j) + x(i) * a(i-j+mu+1,j)
    end do
  end do

  return
end
subroutine r8cbb_add ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8CBB_ADD adds a value to an entry of an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real ( kind = 8 ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8CBB matrix.
!
!    Input, integer ( kind = 4 ) I, J, the indices of the entry to be 
!    incremented.
!
!    Input, real ( kind = 8 ) VALUE, the value to be added to the (I,J) entry.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j
  real ( kind = 8 ) value

  if ( value == 0.0D+00 ) then
    return
  end if
!
!  Check for I or J out of bounds.
!
  if ( i <= 0 .or. n1+n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_ADD - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop
  end if

  if ( j <= 0 .or. n1+n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_ADD - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ',j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( mu < (j-i) .or. ml < (i-j) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CBB_ADD - Warning!'
      write ( *, '(a,i8,a,i8,a)' ) '  Unable to add to entry (', i, ',', j, ').'
      return
    else
      ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
    end if
!
!  The A2 block of the matrix:
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = (ml+mu+1)*n1+(j-n1-1)*n1 + i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2 + (i-n1)
  end if

  a(ij) = a(ij) + value

  return
end
subroutine r8cbb_error ( n1, n2, ml, mu, ierror )

!*****************************************************************************80
!
!! R8CBB_ERROR checks the dimensions of an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1 - 1.
!
!    Output, integer ( kind = 4 ) IERROR, reports whether errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if ML is illegal;
!    IERROR = IERROR + 2 if MU is illegal;
!    IERROR = IERROR + 4 if N1 is illegal;
!    IERROR = IERROR + 8 if N2 is illegal;
!    IERROR = IERROR + 16 if neither N1 nor N2 is positive.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  ierror = 0

  if ( ml < 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal ML = ', ml
    write ( *, '(a)' ) '  but ML must be greater than or equal to 0.'
  else if ( max ( n1 - 1, 0 ) < ml ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal ML = ', ml
    write ( *, '(a,i8)' ) '  but ML must be <= Max ( N1 - 1, 0 ) = ', &
      max ( n1 - 1, 0 )
  end if

  if ( mu < 0  ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal MU = ', mu
    write ( *, '(a)' ) '  but MU must be greater than or equal to 0.'
  else if ( max ( n1 - 1, 0 ) < ml ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal MU = ', mu
    write ( *, '(a,i8)' ) '  but MU must be <= Max ( N1 - 1, 0 ) = ', &
      max ( n1 - 1, 0 )
  end if

  if ( n1 < 0 ) then
    ierror = ierror + 4
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal N1 = ', n1
  end if

  if ( n2 < 0 ) then
    ierror = ierror + 8
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal N2 = ', n2
  end if

  if ( n1 + n2 <= 0 ) then
    ierror = ierror + 16
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal N1+N2 = ', n1+n2
  end if

  return
end
subroutine r8cbb_fa ( n1, n2, ml, mu, a, info )

!*****************************************************************************80
!
!! R8CBB_FA factors an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!    Once the matrix has been factored by SCCB_FA, SCCB_SL may be called
!    to solve linear systems involving the matrix.
!
!    SCCB_FA uses special non-pivoting versions of LINPACK routines to
!    carry out the factorization.  The special version of the banded
!    LINPACK solver also results in a space saving, since no entries
!    need be set aside for fill in due to pivoting.
!
!    The linear system must be border banded, of the form:
!
!      ( A1 A2 ) (X1) = (B1)
!      ( A3 A4 ) (X2)   (B2)
!
!    where A1 is a (usually big) banded square matrix, A2 and A3 are
!    column and row strips which may be nonzero, and A4 is a dense
!    square matrix.
!
!    The algorithm rewrites the system as:
!
!         X1 + inverse(A1) A2 X2 = inverse(A1) B1
!
!      A3 X1 +             A4 X2 = B2
!
!    and then rewrites the second equation as
!
!      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
!
!    The algorithm will certainly fail if the matrix A1 is singular,
!    or requires pivoting.  The algorithm will also fail if the A4 matrix,
!    as modified during the process, is singular, or requires pivoting.
!    All these possibilities are in addition to the failure that will
!    if the total matrix A is singular.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real ( kind = 8 ) A( (ML+MU+1)*N1 + 2*N1*N2 + N2*N2).
!    On input, A contains the compact border-banded coefficient matrix.
!    On output, A contains information describing a partial factorization
!    of the original coefficient matrix.  
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nband

  nband = (ml+mu+1)*n1
!
!  Factor the A1 band matrix, overwriting A1 by its factors.
!
  if ( 0 < n1 ) then

    call r8cb_np_fa ( n1, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CBB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  R8CB_NP_FA returned INFO = ', info
      write ( *, '(a)' ) '  Factoring failed for column INFO.'
      write ( *, '(a)' ) '  The band matrix A1 is singular.'
      write ( *, '(a)' ) '  This algorithm cannot continue!'
      stop
    end if

  end if

  if ( 0 < n1 .and. 0 < n2 ) then
!
!  Set A2 := -inverse(A1) * A2.
!
    a(nband+1:nband+n1*n2) = -a(nband+1:nband+n1*n2)
    job = 0

    do j = 1, n2
      call r8cb_np_sl ( n1, ml, mu, a, a(nband+(j-1)*n1+1), job )
    end do
!
!  Set A4 := A4 + A3*A2
!
    do i = 1, n2
      do j = 1, n1
        ij = nband + n1*n2 + (j-1)*n2 + i
        do k = 1, n2
          ik = nband + 2*n1*n2 + (k-1)*n2 + i
          jk = nband + (k-1)*n1 + j
          a(ik) = a(ik) + a(ij) * a(jk)
        end do
      end do
    end do

  end if
!
!  Factor A4.
!
  if ( 0 < n2 ) then

    call r8ge_np_fa ( n2, a(nband+2*n1*n2+1), info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CBB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  R8GE_NP_FA returned INFO = ',info
      write ( *, '(a)' ) '  This indicates singularity in column INFO'
      info = n1 + info
      write ( *, '(a,i8)' ) '  of the A4 submatrix, which is column ',info
      write ( *, '(a)' ) '  of the full matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  It is possible that the full matrix is '
      write ( *, '(a)' ) '  nonsingular, but the algorithm R8CBB_FA may'
      write ( *, '(a)' ) '  not be used for this matrix.'
      stop
    end if

  end if

  return
end
subroutine r8cbb_get ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8CBB_GET returns the value of an entry of an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8CBB matrix.
!
!    Input, integer ( kind = 4 ) I, J, the row and column of the entry to 
!    retrieve.
!
!    Output, real ( kind = 8 ) VALUE, the value of the (I,J) entry.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j
  real ( kind = 8 ) value
!
!  Check for I or J out of bounds.
!
  if ( i <= 0 .or. n1+n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_GET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop
  end if

  if ( j <= 0 .or. n1+n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_GET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ', j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( mu < (j-i) .or. ml < (i-j) ) then
      value = 0.0D+00
      return
    else
      ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
    end if
!
!  The A2 block of the matrix:
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = (ml+mu+1)*n1+(j-n1-1)*n1+i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
  end if

  value = a(ij)

  return
end
subroutine r8cbb_indicator ( n1, n2, ml, mu, a )

!*****************************************************************************80
!
!! R8CBB_INDICATOR sets up an R8CBB indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Output, real ( kind = 8 ) A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  fac = 10 ** ( i4_log_10 ( n1 + n2 ) + 1 )
!
!  Set the banded matrix A1.
!
  do j = 1, n1
    do row = 1, ml + mu + 1
      i = row + j - mu - 1
      if ( 1 <= i .and. i <= n1 ) then
        a(row+(j-1)*(ml+mu+1)) = real ( fac * i + j, kind = 8 )
      else
        a(row+(j-1)*(ml+mu+1)) = 0.0D+00
      end if
    end do
  end do
!
!  Set the N1 by N2 rectangular strip A2.
!
  base = ( ml + mu + 1 ) * n1

  do i = 1, n1
    do j = n1 + 1, n1 + n2
      a(base + i + (j-n1-1)*n1 ) = real ( fac * i + j, kind = 8 )
    end do
  end do
!
!  Set the N2 by N1 rectangular strip A3.
!
  base = ( ml + mu + 1 ) * n1 + n1 * n2

  do i = n1 + 1, n1 + n2
    do j = 1, n1    
      a(base + i-n1 + (j-1)*n2 ) = real ( fac * i + j, kind = 8 )
    end do
  end do
!
!  Set the N2 by N2 square A4.
!
  base = ( ml + mu + 1 ) * n1 + n1 * n2 + n2 * n1

  do i = n1 + 1, n1 + n2
    do j = n1 + 1, n1 + n2
      a(base + i-n1 + (j-n1-1)*n2 ) = real ( fac * i + j, kind = 8 )
    end do
  end do

  return
end
subroutine r8cbb_mxv ( n1, n2, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8CBB_MXV multiplies an R8CBB matrix by an R8VEC.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, real ( kind = 8 ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8CBB matrix.
!
!    Input, real ( kind = 8 ) X(N1+N2), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N1+N2), the result of multiplying A by X.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) b(n1+n2)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n1+n2)
!
!  Set B to zero.
!
  b(1:n1+n2) = 0.0D+00
!
!  Multiply by A1.
!
  do j = 1, n1
    ilo = max ( 1, j-mu )
    ihi = min ( n1, j+ml )
    ij = (j-1)*(ml+mu+1)-j+mu+1
    b(ilo:ihi) = b(ilo:ihi) + a(ij+ilo:ij+ihi) * x(j)
  end do
!
!  Multiply by A2.
!
  do j = n1+1, n1+n2
    ij = (ml+mu+1)*n1+(j-n1-1)*n1
    b(1:n1) = b(1:n1) + a(ij+1:ij+n1) * x(j)
  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1+n2
    ij = (ml+mu+1)*n1+n1*n2+(j-1)*n2-n1
    b(n1+1:n1+n2) = b(n1+1:n1+n2) + a(ij+n1+1:ij+n1+n2) * x(j)
  end do

  return
end
subroutine r8cbb_print ( n1, n2, ml, mu, a, title )

!*****************************************************************************80
!
!! R8CBB_PRINT prints an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  character ( len = * ) title

  call r8cbb_print_some ( n1, n2, ml, mu, a, 1, 1, n1+n2, n1+n2, title )

  return
end
subroutine r8cbb_print_some ( n1, n2, ml, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8CBB_PRINT_SOME prints some of an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n1+n2 )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n1+n2 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0D+00

        if ( i <= n1 .and. j <= n1 ) then
          if ( j - i <= mu .and. i - j <= ml ) then
            ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
            aij = a(ij)
          end if
        else if ( i <= n1 .and. n1 < j ) then
          ij = (ml+mu+1)*n1+(j-n1-1)*n1+i
          aij = a(ij)
        else if ( n1 < i ) then
          ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
          aij = a(ij)
        end if

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8cbb_random ( n1, n2, ml, mu, seed, a )

!*****************************************************************************80
!
!! R8CBB_RANDOM randomizes an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2), 
!    the R8CBB matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  integer ( kind = 4 ) row
  integer ( kind = 4 ) seed
!
!  Randomize the banded matrix A1.
!  We still believe that the "junk" entries should be set to 0.
!
  do j = 1, n1
    do row = 1, ml+mu+1
      i = row + j - mu - 1
      if ( 1 <= i .and. i <= n1 ) then
        r = r8_uniform_01 ( seed )
      else
        r = 0.0D+00
      end if
      a(row+(j-1)*(ml+mu+1)) = r
    end do
  end do
!
!  Randomize the rectangular strips A2+A3+A4.
!
  ilo = (ml+mu+1) * n1 + 1

  call r8vec_uniform_01 ( n1*n2+n2*n1+n2*n2, seed, a(ilo:) )

  return
end
subroutine r8cbb_set ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8CBB_SET sets the value of an entry in an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real ( kind = 8 ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8CBB matrix.
!
!    Input, integer ( kind = 4 ) I, J, the row and column of the entry to set.
!
!    Input, real ( kind = 8 ) VALUE, the value to be assigned to the
!    (I,J) entry.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j
  real ( kind = 8 ) value
!
!  Check for I or J out of bounds.
!
  if ( i <= 0 .or. n1+n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop
  end if

  if ( j <= 0 .or. n1+n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ', j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( mu < (j-i) .or. ml < (i-j) ) then
      if ( value /= 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CBB_SET - Warning!'
        write ( *, '(a,i8,a,i8,a)' ) '  Unable to set entry (', i, ',', j, ').'
      end if
      return
    else
      ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
    end if
!
!  The A2 block of the matrix:
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = (ml+mu+1)*n1 + (j-n1-1)*n1 + i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = (ml+mu+1)*n1 + n2*n1 + (j-1)*n2 + (i-n1)
  end if

  a(ij) = value

  return
end
subroutine r8cbb_sl ( n1, n2, ml, mu, a_lu, b )

!*****************************************************************************80
!
!! R8CBB_SL solves an R8CBB system factored by R8CBB_FA.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!
!    The linear system A * x = b is decomposable into the block system:
!
!      ( A1 A2 ) * (X1) = (B1)
!      ( A3 A4 )   (X2)   (B2)
!
!    where A1 is a (usually big) banded square matrix, A2 and A3 are
!    column and row strips which may be nonzero, and A4 is a dense
!    square matrix.
!
!    All the arguments except B are input quantities only, which are
!    not changed by the routine.  They should have exactly the same values
!    they had on exit from R8CBB_FA.
!
!    If more than one right hand side is to be solved, with the same
!    matrix, R8CBB_SL should be called repeatedly.  However, R8CBB_FA only
!    needs to be called once to create the factorization.
!
!    See the documentation of R8CBB_FA for details on the matrix storage.
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
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A_LU( (ML+MU+1)*N1 + 2*N1*N2 + N2*N2).
!    the LU factors from R8CBB_FA.
!
!    Input/output, real ( kind = 8 ) B(N1+N2).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a_lu((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) b(n1+n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) nband

  nband = (ml+mu+1)*n1
!
!  Set B1 := inverse(A1) * B1.
!
  if ( 0 < n1 ) then
    job = 0
    call r8cb_np_sl ( n1, ml, mu, a_lu, b, job )
  end if
!
!  Modify the right hand side of the second linear subsystem.
!  Replace B2 by B2-A3*B1.
!
  do j = 1, n1
    do i = 1, n2
      ij = nband + n1*n2 + (j-1)*n2 + i
      b(n1+i) = b(n1+i) - a_lu(ij) * b(j)
    end do
  end do
!
!  Solve A4*B2 = B2.
!
  if ( 0 < n2 ) then
    job = 0
    call r8ge_np_sl ( n2, a_lu(nband+2*n1*n2+1), b(n1+1), job )
  end if
!
!  Modify the first subsolution.
!  Set B1 = B1+A2*B2.
!
  do i = 1, n1
    do j = 1, n2
      ij = nband + (j-1)*n1 + i
      b(i) = b(i) + a_lu(ij) * b(n1+j)
    end do
  end do

  return
end
subroutine r8cbb_to_r8ge ( n1, n2, ml, mu, a, b )

!*****************************************************************************80
!
!! R8CBB_TO_R8GE copies an R8CBB matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = 8 ) A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
!
!    Output, real ( kind = 8 ) B(N1+N2,N1+N2), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) b(n1+n2,n1+n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j

  do i = 1, n1
    do j = 1, n1

      if ( mu+ml < (j-i) .or. ml < (i-j) ) then
        b(i,j) = 0.0D+00
      else
        ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
        b(i,j) = a(ij)
      end if

    end do
  end do

  do i = 1, n1
    do j = n1+1, n2
      ij = (ml+mu+1)*n1+(j-n1-1)*n1+i
      b(i,j) = a(ij)
    end do
  end do

  do i = n1+1, n2
    do j = 1, n1+n2
      ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
      b(i,j) = a(ij)
    end do
  end do

  return
end
subroutine r8cbb_vxm ( n1, n2, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8CBB_VXM multiplies an R8VEC by an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, integer ( kind = 4 ) N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, real ( kind = 8 ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8CBB matrix.
!
!    Input, real ( kind = 8 ) X(N1+N2), the vector to multiply the matrix.
!
!    Output, real ( kind = 8 ) B(N1+N2), the product X * A.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = 8 ) b(n1+n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n1+n2)
!
!  Set B to zero.
!
  b(1:n1+n2) = 0.0D+00
!
!  Multiply by A1.
!
  do j = 1, n1
    ilo = max ( 1, j-mu )
    ihi = min ( n1, j+ml )
    ij = (j-1)*(ml+mu+1)-j+mu+1
    b(j) = b(j) + sum ( x(ilo:ihi) * a(ij+ilo:ij+ihi) )
  end do
!
!  Multiply by A2.
!
  do j = n1+1, n1+n2
    ij = (ml+mu+1)*n1+(j-n1-1)*n1
    b(j) = b(j) + sum ( x(1:n1) * a(ij+1:ij+n1) )
  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1+n2
    ij = (ml+mu+1)*n1+n1*n2+(j-1)*n2-n1
    do i = n1+1, n1+n2
      b(j) = b(j) + x(i) * a(ij+i)
    end do
  end do

  return
end
subroutine r8cc_get ( m, n, nz_num, col, row, a, i, j, aij )

!*****************************************************************************80
!
!! R8CC_GET gets a value of an R8CC matrix.
!
!  Discussion:
!
!    It is legal to request entries of the matrix for which no storage
!    was set aside.  In that case, a zero value will be returned.
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) COL(N+1), indicate where each column's data
!    begins.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), the row indices.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero entries.
!
!    Input, integer ( kind = 4 ) I, J, the indices of the value to retrieve.
!
!    Output, real ( kind = 8 ) AIJ, the value of A(I,J).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) col(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row(nz_num)
!
!  Seek sparse index K corresponding to full index (I,J).
!
  call r8cc_ijk ( m, n, nz_num, col, row, i, j, k )
!
!  If no K was found, then be merciful, and simply return 0.
!
  if ( k == -1 ) then
    aij = 0.0D+00
  else
    aij = a(k)
  end if

  return
end
subroutine r8cc_ijk ( m, n, nz_num, col, row, i, j, k )

!*****************************************************************************80
!
!! R8CC_IJK seeks the sparse index K of (I,J), the full index of an R8CC matrix.
!
!  Discussion:
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) COL(N+1), indicate where each column's 
!    data begins.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), the row indices.
!
!    Input, integer ( kind = 4 ) I, J, the indices of the value to retrieve.
!
!    Output, integer ( kind = 4 ) K, the index of the sparse matrix in which 
!    entry (I,J) is stored, or -1 if no such entry exists.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) col(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row(nz_num)
!
!  Determine the part of ROW containing row indices of entries
!  in column J.
!
  k1 = col(j)
  k2 = col(j+1)-1
!
!  Seek the location K for which ROW(K) = I.
!  
  call i4vec_search_binary_a ( k2+1-k1, row(k1:k2), i, k )

  if ( k /= -1 ) then
    k = k + k1 - 1
  end if

  return
end
subroutine r8cc_inc ( m, n, nz_num, col, row, a, i, j, aij )

!*****************************************************************************80
!
!! R8CC_INC increments a value of an R8CC matrix.
!
!  Discussion:
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) COL(N+1), indicate where each column's 
!    data begins.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), the row indices.
!
!    Input/output, real ( kind = 8 ) A(NZ_NUM), the nonzero entries.
!    On output, entry (I,J) has been incremented.
!
!    Input, integer ( kind = 4 ) I, J, the indices of the value to retrieve.
!
!    Input, real ( kind = 8 ) AIJ, the value to be added to A(I,J).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) col(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row(nz_num)
!
!  Seek sparse index K corresponding to full index (I,J).
!
  call r8cc_ijk ( m, n, nz_num, col, row, i, j, k )
!
!  If no K was found, we fail.
!
  if ( k == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CC_INC - Fatal error!'
    write ( *, '(a)' ) '  R8CC_IJK could not find the entry.'
    write ( *, '(a,i8)' ) '  Row I = ', i
    write ( *, '(a,i8)' ) '  Col J = ', j
    stop
  end if

  a(k) = a(k) + aij

  return
end
subroutine r8cc_indicator ( m, n, nz_num, col, row, a )

!*****************************************************************************80
!
!! R8CC_INDICATOR sets up an R8CC indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer ( kind = 4 ) COL(N+1), points to the first element of 
!    each column.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), contains the row indices 
!    of the elements.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the R8CC matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(n+1)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row(nz_num)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do j = 1, n
    do k = col(j), col(j+1) - 1
      i = row(k)
      a(k) = real ( fac * i + j, kind = 8 )
    end do
  end do

  return
end
subroutine r8cc_kij ( m, n, nz_num, col, row, k, i, j )

!*****************************************************************************80
!
!! R8CC_KIJ seeks the full index (I,J) of K, the sparse index of an R8CC matrix.
!
!  Discussion:
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) COL(N+1), indicate where each column's data 
!    begins.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), the row indices.
!
!    Input, integer ( kind = 4 ) K, the sparse index of an entry of the matrix.
!    1 <= K <= NZ_NUM.
!
!    Output, integer ( kind = 4 ) I, J, the full indices corresponding to the 
!    sparse index K.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) col(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row(nz_num)

  i = -1
  j = -1

  if ( k < 1 .or. nz_num < k ) then
    return
  end if
!
!  The row index is easy.
!
  i = row(k)
!
!  Determine the column by bracketing in COL.
!
  do jj = 1, n
    k1 = col(jj)
    k2 = col(jj+1)-1
    if ( k1 <= k .and. k <= k2 ) then
      j = jj
      exit
    end if
  end do

  if ( j == -1 ) then
    return
  end if

  return
end
subroutine r8cc_mxv ( m, n, nz_num, col, row, a, x, b )

!*****************************************************************************80
!
!! R8CC_MXV multiplies an R8CC matrix by an R8VEC.
!
!  Discussion:
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer ( kind = 4 ) COL(N+1), points to the first element of 
!    each column.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), contains the row indices of 
!    the elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the R8CC matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(M), the product A*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) col(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) x(n)

  b(1:m) = 0.0D+00

  do j = 1, n
    do k = col(j), col(j+1) - 1
      i = row(k)
      b(i) = b(i) + a(k) * x(j)
    end do
  end do

  return
end
subroutine r8cc_print ( m, n, nz_num, col, row, a, title )

!*****************************************************************************80
!
!! R8CC_PRINT prints an R8CC matrix.
!
!  Discussion:
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer ( kind = 4 ) COL(N+1), points to the first element of 
!    each column.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), contains the row indices of 
!    the elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the R8CC matrix.
!
!    Input, character ( len = * ) TITLE, a title.
! 
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(n+1)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) title

  call r8cc_print_some ( m, n, nz_num, col, row, a, 1, 1, n, n, title )

  return
end
subroutine r8cc_print_some ( m, n, nz_num, col, row, a, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! R8CC_PRINT_SOME prints some of an R8CC matrix.
!
!  Discussion:
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer ( kind = 4 ) COL(N+1), points to the first element of 
!    each column.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), contains the row indices of 
!    the elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the R8CC matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) col(n+1)
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
!  1) Assume everything is zero.
!
      aij = 0.0D+00
      do j2 = 1, inc
        write ( ctemp(j2), '(f8.0,6x)' ) aij
      end do
!
!  2) Now consider each column J in J2LO to J2HI, 
!     and look at every nonzero, and check if it occurs in row I.
!
      do j = j2lo, j2hi
        do k = col(j), col(j+1)-1
          if ( row(k) == i ) then
            j2 = j - j2lo + 1
            if ( r8_is_int ( a(k) ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) a(k)
            else
              write ( ctemp(j2), '(g14.6)' ) a(k)
            end if
          end if
        end do
      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8cc_random ( m, n, nz_num, col, row, a, seed )

!*****************************************************************************80
!
!! R8CC_RANDOM randomizes an R8CC matrix.
!
!  Discussion:
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer ( kind = 4 ) COL(N+1), points to the first element of each 
!    column.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), contains the row indices of the 
!    elements.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the R8CC matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(n+1)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row(nz_num)
  integer ( kind = 4 ) seed

  do j = 1, n
    do k = col(j), col(j+1) - 1
      a(k) = r8_uniform_01 ( seed )
    end do
  end do

  return
end
subroutine r8cc_read ( col_file, row_file, a_file, m, n, nz_num, col, row, a )

!*****************************************************************************80
!
!! R8CC_READ reads an R8CC matrix from three files.
!
!  Discussion:
!
!    This routine needs the values of M, N, and NZ_NUM, which can be 
!    determined by a call to R8CC_READ_SIZE.
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, character ( len = * ) COL_FILE, ROW_FILE, A_FILE, the names of the 
!    files containing the column pointers, row indices, and matrix entries.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the 
!    matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Output, integer ( kind = 4 ) COL(N+1), the column pointers.
!
!    Output, integer ( kind = 4 ) ROW(NZ_NUM), the row indices.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  character ( len = * ) a_file
  integer ( kind = 4 ) col(n+1)
  character ( len = * ) col_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) row_file

  call get_unit ( input_unit )
!
!  Read the column information.
!
  open ( unit = input_unit, file = col_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( col_file ) // '".'
    stop
  end if

  do k = 1, n + 1

    read ( input_unit, *, iostat = ios ) col(k)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
      write ( *, '(a,i8,a)' ) '  I/O error while reading record ', k, &
        ' of "' // trim ( col_file ) // '".'
      stop
    end if

  end do

  close ( unit = input_unit )
!
!  Read the row information.
!
  open ( unit = input_unit, file = row_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( row_file ) // '".'
    stop
  end if

  do k = 1, nz_num

    read ( input_unit, *, iostat = ios ) row(k)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
      write ( *, '(a,i8,a)' ) '  I/O error while reading record ', k, &
        ' of "' // trim ( row_file ) // '".'
      stop
    end if

  end do

  close ( unit = input_unit )
!
!  Read the value information.
!
  open ( unit = input_unit, file = a_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( a_file ) // '".'
    stop
  end if

  do k = 1, nz_num

    read ( input_unit, *, iostat = ios ) a(k)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
      write ( *, '(a,i8,a)' ) '  I/O error while reading record ', k, &
        ' of "' // trim ( a_file ) // '".'
      stop
    end if

  end do

  close ( unit = input_unit )

  return
end
subroutine r8cc_read_size ( col_file, row_file, m, n, nz_num, base )

!*****************************************************************************80
!
!! R8CC_READ_SIZE reads the sizes of an R8CC sparse matrix from a file.
!
!  Discussion:
!
!    The value of M is "guessed" to be the largest value that occurs in
!    the ROW file.  However, if a row index of 0 is encountered, then
!    the value of M is incremented by 1.
!
!    The value of N is the number of records in the COL file minus 1.
!
!    The value of NZ_NUM is simply the number of records in the ROW file.
!
!    The value of BASE is 0 or 1, depending on whether the program
!    "guesses" that the row and column indices are 0-based or 1-based.
!    Although the first entry of the COL array might be used as evidence,
!    this program makes its determination based on whether it encounters
!    a 0 index in the ROW file.
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, character ( len = * ) COL_FILE, ROW_FILE, the names of the 
!    column and row files that describe the structure of the matrix.
!
!    Output, integer ( kind = 4 ) M, N, the inferred number of rows and columns 
!    in the sparse matrix.
!
!    Output, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries in the
!    sparse matrix.
!
!    Output, integer ( kind = 4 ) BASE, is 0 if the row indexing is believed
!    to be 0-based, and 1 if the row-index is believed to be
!    1-based.  In uncertain cases, BASE = 1 is the default.
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) col
  character ( len = * ) col_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ) row
  character ( len = * ) row_file
!
!  Default values.
!
  m = -1
  n = -1
  nz_num = -1
  base = -1
!
!  Check the COL file first.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = col_file, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CC_READ_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open  the input file "' &
      // trim ( col_file ) // '".'
    stop
  end if

  n = -1

  do

    read ( input_unit, *, iostat = ios ) col

    if ( ios /= 0 ) then
      exit
    end if

    n = n + 1

  end do

  close ( unit = input_unit )
!
!  Check the ROW file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = row_file, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CC_READ_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( row_file ) // '".'
    stop
  end if

  base = 1
  m = 0
  nz_num = 0
  
  do

    read ( input_unit, *, iostat = ios ) row

    if ( ios /= 0 ) then
      exit
    end if

    nz_num = nz_num + 1
    m = max ( m, row )
    if ( row == 0 ) then
      base = 0
    end if

  end do

  close ( unit = input_unit )

  return
end
subroutine r8cc_set ( m, n, nz_num, col, row, a, i, j, aij )

!*****************************************************************************80
!
!! R8CC_SET sets a value of an R8CC matrix.
!
!  Discussion:
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) COL(N+1), indicate where each column's data 
!    begins.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), the row indices.
!
!    Input/output, real ( kind = 8 ) A(NZ_NUM), the nonzero entries.
!    On output, the entry of A corresponding to (I,J) has been reset.
!
!    Input, integer ( kind = 4 ) I, J, the indices of the value to retrieve.
!
!    Input, real ( kind = 8 ) AIJ, the new value of A(I,J).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) col(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row(nz_num)
!
!  Seek sparse index K corresponding to full index (I,J).
!
  call r8cc_ijk ( m, n, nz_num, col, row, i, j, k )
!
!  If no K was found, we fail.
!
  if ( k == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CC_SET - Fatal error!'
    write ( *, '(a)' ) '  R8CC_IJK could not find the entry.'
    write ( *, '(a,i8)' ) '  Row I = ', i
    write ( *, '(a,i8)' ) '  Col J = ', j
    stop
  end if

  a(k) = aij

  return
end
subroutine r8cc_to_r8ge ( m, n, nz_num, col, row, a, b )

!*****************************************************************************80
!
!! R8CC_TO_R8GE converts an R8CC matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer ( kind = 4 ) COL(N+1), points to the first element of 
!    each column.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), contains the row indices 
!    of the elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the R8CC matrix.
!
!    Output, real ( kind = 8 ) B(M,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(m,n)
  integer ( kind = 4 ) col(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)

  b(1:m,1:n) = 0.0D+00

  do j = 1, n
    do k = col(j), col(j+1)-1
      b(row(k),j) = a(k)
    end do
  end do

  return
end
subroutine r8cc_vxm ( m, n, nz_num, col, row, a, x, b )

!*****************************************************************************80
!
!! R8CC_VXM multiplies an R8VEC times an R8CC matrix.
!
!  Discussion:
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer ( kind = 4 ) COL(N+1), points to the first element 
!    of each column.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), contains the row indices 
!    of the elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the R8CC matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(N), the product A'*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) col(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) x(m)

  b(1:n) = 0.0D+00

  do j = 1, n
    do k = col(j), col(j+1) - 1
      i = row(k)
      b(j) = b(j) + a(k) * x(i)
    end do
  end do

  return
end
subroutine r8cc_write ( col_file, row_file, a_file, m, n, nz_num, col, row, a )

!*****************************************************************************80
!
!! R8CC_WRITE writes an R8CC matrix to three files.
!
!  Discussion:
!
!    The R8CC format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The R8CC format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, character ( len = * ) COL_FILE, ROW_FILE, A_FILE, the names of the 
!    files containing the column pointers, row entries, and matrix entries.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements 
!    in the matrix.
!
!    Input, integer ( kind = 4 ) COL(N+1), the column pointers.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), the row indices.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  character ( len = * ) a_file
  integer ( kind = 4 ) col(n+1)
  character ( len = * ) col_file
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) row_file

  call get_unit ( output_unit )
!
!  Write the column information.
!
  open ( unit = output_unit, file = col_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CC_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' &
      // trim ( col_file ) // '".'
    stop
  end if

  do k = 1, n + 1

    write ( output_unit, '(i8)' ) col(k)

  end do

  close ( unit = output_unit )
!
!  Write the row information.
!
  open ( unit = output_unit, file = row_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CC_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' &
      // trim ( row_file ) // '".'
    stop
  end if

  do k = 1, nz_num

    write ( output_unit, '(i8)' ) row(k)

  end do

  close ( unit = output_unit )
!
!  Write the value information.
!
  open ( unit = output_unit, file = a_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CC_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' &
      // trim ( a_file ) // '".'
    stop
  end if

  do k = 1, nz_num

    write ( output_unit, '(g14.6)' ) a(k)

  end do

  close ( unit = output_unit )

  return
end
subroutine r8ci_eval ( n, a, lambda )

!*****************************************************************************80
!
!! R8CI_EVAL returns the eigenvalues of an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis,
!    Circulant Matrices,
!    Second Edition,
!    Chelsea, 1994,
!    ISBN: 0828403384,
!    LC: QA188.D37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8CI matrix.
!
!    Output, complex ( kind = 8 ) LAMBDA(N), the eigenvalues.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  complex ( kind = 8 ) lambda(n)
  complex ( kind = 8 ) w(n)

  call c8vec_unity ( n, w )

  lambda(1:n) = cmplx ( a(n), 0.0D+00, kind = 8 )
  do i = n-1, 1, -1
    lambda(1:n) = lambda(1:n) * w(1:n) + cmplx ( a(i), 0.0D+00, kind = 8 )
  end do

  call c8vec_sort_a2 ( n, lambda )

  return
end
subroutine r8ci_indicator ( n, a )

!*****************************************************************************80
!
!! R8CI_INDICATOR sets up an R8CI indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 January 2004
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
!    Output, real ( kind = 8 ) A(N), the R8CI matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  i = 1

  do j = 1, n
    a(j) = real ( fac * i + j, kind = 8 )
  end do

  return
end
subroutine r8ci_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! R8CI_MXV multiplies an R8CI matrix by an R8VEC.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8CI matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    b(i) = sum ( a(n+2-i:n) * x(1:i-1) ) &
         + sum ( a(1:n+1-i) * x(i:n) )
  end do

  return
end
subroutine r8ci_print ( n, a, title )

!*****************************************************************************80
!
!! R8CI_PRINT prints an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
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
!    Input, real ( kind = 8 ) A(N), the R8CI matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  character ( len = * ) title

  call r8ci_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r8ci_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8CI_PRINT_SOME prints some of an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
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
!    Input, real ( kind = 8 ) A(N), the R8CI matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(j+1-i)
        else
          aij = a(n+j+1-i)
        end if

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8ci_random ( n, seed, a )

!*****************************************************************************80
!
!! R8CI_RANDOM randomizes an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2003
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
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(N), the R8CI matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  call r8vec_uniform_01 ( n, seed, a )

  return
end
subroutine r8ci_sl ( n, a, b, x, job )

!*****************************************************************************80
!
!! R8CI_SL solves an R8CI system.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
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
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8CI matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) r5
  real ( kind = 8 ) r6
  real ( kind = 8 ) work(2*n-2)
  real ( kind = 8 ) x(n)

  if ( job == 0 ) then
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = 0.0D+00
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(n+2-nsub)
      r6 = a(nsub)

      if ( 2 < nsub ) then

        work(nsub-1) = r2

        do i = 1, nsub-2
          r5 = r5 + a(n+1-i) * work(nsub-i)
          r6 = r6 + a(i+1) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( 2 < nsub ) then

        r6 = work(n)
        work(n-1+nsub-1) = 0.0D+00
        do i = 2, nsub-1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = 0.0D+00
      do i = 1, nsub-1
        r5 = r5 + a(n+1-i) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      x(1:nsub-1) = x(1:nsub-1) + work(n:n+nsub-2) * r6
      x(nsub) = r6

    end do

  else
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = 0.0D+00
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(nsub)
      r6 = a(n+2-nsub)

      if ( 2 < nsub ) then

        work(nsub-1) = r2

        do i = 1, nsub-2
          r5 = r5 + a(i+1) * work(nsub-i)
          r6 = r6 + a(n+1-i) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( 2 < nsub ) then

        r6 = work(n)
        work(n-1+nsub-1) = 0.0D+00
        do i = 2, nsub-1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = 0.0D+00
      do i = 1, nsub-1
        r5 = r5 + a(i+1) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      do i = 1, nsub-1
        x(i) = x(i) + work(n-1+i) * r6
      end do

      x(nsub) = r6

    end do

  end if

  return
end
subroutine r8ci_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R8CI_TO_R8GE copies an R8CI matrix into an R8GE matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8CI matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i

  do i = 1, n
    b(i,1:i-1) = a(n+2-i:n+2*1-2*i)
    b(i,i:n) = a(1:n+1-i)
  end do

  return
end
subroutine r8ci_vxm ( n, a, x, b )

!*****************************************************************************80
!
!! R8CI_VXM multiplies an R8VEC by an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8CI matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A' * X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    b(i) = sum ( a(i:1:-1) * x(1:i) ) &
         + sum ( a(n:i+1:-1) * x(i+1:n) )
  end do

  return
end
subroutine r8col_swap ( m, n, a, i, j )

!*****************************************************************************80
!
!! R8COL_SWAP swaps columns I and J of a real array of column data.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!    Output:
!
!      A = (
!        1.  4.  3.  2.
!        5.  8.  7.  6.
!        9. 12. 11. 10. )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  if ( 1 <= i .and. i <= n .and. 1 <= j .and. j <= n ) then

    do k = 1, m
      call r8_swap ( a(k,i), a(k,j) )
    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    write ( *, '(a,i8)' ) '  I =    ', i
    write ( *, '(a,i8)' ) '  J =    ', j
    write ( *, '(a,i8)' ) '  NCOL = ', n
    stop

  end if

  return
end
subroutine r8gb_det ( n, ml, mu, a_lu, pivot, det )

!*****************************************************************************80
!
!! R8GB_DET computes the determinant of a matrix factored by R8GB_FA or R8GB_TRF.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
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
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A_LU(2*ML+MU+1,N), the LU factors from 
!    R8GB_FA or R8GB_TRF.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector, as computed 
!    by R8GB_FA or R8GB_TRF.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(2*ml+mu+1,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pivot(n)

  det = product ( a_lu(ml+mu+1,1:n) )

  do i = 1, n
    if ( pivot(i) /= i ) then
      det = -det
    end if
  end do

  return
end
subroutine r8gb_fa ( n, ml, mu, a, pivot, info )

!*****************************************************************************80
!
!! R8GB_FA performs a LINPACK-style PLU factorization of an R8GB matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    This routine is based on the LINPACK routine SGBFA.
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    The following program segment will set up the input.
!
!      m = ml + mu + 1
!      do j = 1, n
!        i1 = max ( 1, j-mu )
!        i2 = min ( n, j+ml )
!        do i = i1, i2
!          k = i - j + m
!          a(k,j) = afull(i,j)
!        end do
!      end do
!
!    This uses rows ML+1 through 2*ML+MU+1 of the array A.
!    In addition, the first ML rows in the array are used for
!    elements generated during the triangularization.
!
!    The ML+MU by ML+MU upper left triangle and the
!    ML by ML lower right triangle are not referenced.
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
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt
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
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input/output, real ( kind = 8 ) A(2*ML+MU+1,N), on input, 
!    the matrix in band storage, on output, information about 
!    the LU factorization.
!
!    Output, integer ( kind = 4 ) PIVOT(N), the pivot vector.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) jz
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  real ( kind = 8 ) t

  m = ml + mu + 1
  info = 0
!
!  Zero out the initial fill-in columns.
!
  j0 = mu + 2
  j1 = min ( n, m ) - 1

  do jz = j0, j1
    i0 = m + 1 - jz
    a(i0:ml,jz) = 0.0D+00
  end do

  jz = j1
  ju = 0

  do k = 1, n-1
!
!  Zero out the next fill-in column.
!
    jz = jz + 1
    if ( jz <= n ) then
      a(1:ml,jz) = 0.0D+00
    end if
!
!  Find L = pivot index.
!
    lm = min ( ml, n-k )

    l = m
    do j = m+1, m+lm
      if ( abs ( a(l,k) ) < abs ( a(j,k) ) ) then
        l = j
      end if
    end do

    pivot(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Interchange if necessary.
!
    t      = a(l,k)
    a(l,k) = a(m,k)
    a(m,k) = t
!
!  Compute multipliers.
!
    a(m+1:m+lm,k) = - a(m+1:m+lm,k) / a(m,k)
!
!  Row elimination with column indexing.
!
    ju = max ( ju, mu+pivot(k) )
    ju = min ( ju, n )
    mm = m

    do j = k+1, ju

      l = l - 1
      mm = mm - 1

      if ( l /= mm ) then
        t       = a(l,j)
        a(l,j)  = a(mm,j)
        a(mm,j) = t
      end if

      a(mm+1:mm+lm,j) = a(mm+1:mm+lm,j) + a(mm,j) * a(m+1:m+lm,k)

    end do

  end do

  pivot(n) = n

  if ( a(m,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GB_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop
  end if

  return
end
subroutine r8gb_indicator ( m, n, ml, mu, a )

!*****************************************************************************80
!
!! R8GB_INDICATOR sets up an R8GB indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    Note that the R8GB storage format includes extra room for
!    fillin entries that occur during Gauss elimination.  These entries
!    are not normally seen or used by the user.  This routine will
!    set those values to zero.
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Output, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) value

  fac = 10 ** ( i4_log_10 ( n ) + 1 )
  k = 0

  do j = 1, n
    do diag = 1, 2 * ml + mu + 1

      i = diag + j - ml - mu - 1

      if ( 1 <= i .and. i <= m .and. i - ml <= j .and. j <= i + mu ) then
        value = real ( fac * i + j, kind = 8 )
      else if ( 1 <= i .and. i <= m .and. &
        i - ml <= j .and. j <= i + mu + ml ) then
        value = 0.0D+00
      else
        k = k + 1
        value = - real ( k, kind = 8 )
      end if

      a(diag,j) = value

    end do
  end do

  return
end
subroutine r8gb_ml ( n, ml, mu, a_lu, pivot, x, b, job )

!*****************************************************************************80
!
!! R8GB_ML computes A * x or A' * X, using R8GB_FA factors.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    It is assumed that R8GB_FA has overwritten the original matrix
!    information by LU factors.  R8GB_ML is able to reconstruct the
!    original matrix from the LU factor data.
!
!    R8GB_ML allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 1998
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
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A_LU(2*ML+MU+1,N), the LU factors from R8GB_FA.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector computed by R8GB_FA.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(N), the result of the multiplication.
!
!    Input, integer ( kind = 4 ) JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * X.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(2*ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      ilo = max ( 1, j - ml - mu )
      do i = ilo, j - 1
        b(i) = b(i) + a_lu(i-j+ml+mu+1,j) * b(j)
      end do
      b(j) = a_lu(j-j+ml+mu+1,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n-1, 1, -1

      ihi = min ( n, j + ml )
      do i = j+1, ihi
        b(i) = b(i) - a_lu(i-j+ml+mu+1,j) * b(j)
      end do

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

    end do

  else
!
!  Y = ( PL )' * X.
!
    do j = 1, n-1

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

      jhi = min ( n, j + ml )
      do i = j + 1, jhi
        b(j) = b(j) - b(i) * a_lu(i-j+ml+mu+1,j)
      end do

    end do
!
!  B = U' * Y = ( PL * U )' * X = A' * X.
!
    do i = n, 1, -1

      jhi = min ( n, i + ml + mu )
      do j = i + 1, jhi
        b(j) = b(j) + b(i) * a_lu(i-j+ml+mu+1,j)
      end do
      b(i) = b(i) * a_lu(i-i+ml+mu+1,i)
    end do

  end if

  return
end
subroutine r8gb_mu ( n, ml, mu, a_lu, pivot, x, b, job )

!*****************************************************************************80
!
!! R8GB_MU computes A * x or A' * X, using R8GB_TRF factors.
!
!  Warning:
!
!    This routine needs to be updated to allow for rectangular matrices.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    It is assumed that R8GB_TRF has overwritten the original matrix
!    information by LU factors.  R8GB_MU is able to reconstruct the
!    original matrix from the LU factor data.
!
!    R8GB_MU allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A_LU(2*ML+MU+1,N), the LU factors from R8GB_TRF.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector computed by R8GB_TRF.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(N), the result of the multiplication.
!
!    Input, integer ( kind = 4 ) JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * X.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(2*ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      ilo = max ( 1, j - ml - mu )
      do i = ilo, j - 1
        b(i) = b(i) + a_lu(i-j+ml+mu+1,j) * b(j)
      end do
      b(j) = a_lu(j-j+ml+mu+1,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n-1, 1, -1

      ihi = min ( n, j + ml )
      do i = j+1, ihi
        b(i) = b(i) + a_lu(i-j+ml+mu+1,j) * b(j)
      end do

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

    end do

  else
!
!  Y = ( PL )' * X.
!
    do j = 1, n-1

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

      jhi = min ( n, j + ml )
      do i = j+1, jhi
        b(j) = b(j) + b(i) * a_lu(i-j+ml+mu+1,j)
      end do

    end do
!
!  B = U' * Y = ( PL * U )' * X = A' * X.
!
    do i = n, 1, -1

      jhi = min ( n, i + ml + mu )
      do j = i + 1, jhi
        b(j) = b(j) + b(i) * a_lu(i-j+ml+mu+1,j)
      end do
      b(i) = b(i) * a_lu(i-i+ml+mu+1,i)
    end do

  end if

  return
end
subroutine r8gb_mxv ( m, n, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8GB_MXV multiplies an R8GB matrix by an R8VEC.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    LINPACK and LAPACK storage of general band matrices requires
!    an extra ML upper diagonals for possible fill in entries during
!    Gauss elimination.  This routine does not access any entries
!    in the fill in diagonals, because it assumes that the matrix
!    has NOT had Gauss elimination applied to it.  If the matrix
!    has been Gauss eliminated, then the routine R8GB_MU must be
!    used instead.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 1999
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
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(M), the product A * x.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  real ( kind = 8 ) x(n)

  do i = 1, m
    b(i) = 0.0D+00
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      b(i) = b(i) + a(i-j+ml+mu+1,j) * x(j)
    end do
  end do

  return
end
subroutine r8gb_nz_num ( m, n, ml, mu, a, nz_num )

!*****************************************************************************80
!
!! R8GB_NZ_NUM counts the nonzeroes in an R8GB matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will examine
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Output, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries in A.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) nz_num

  nz_num = 0

  do i = 1, m
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu + ml )
    do j = jlo, jhi
      if ( a(i-j+ml+mu+1,j) /= 0.0D+00 ) then
        nz_num = nz_num + 1
      end if
    end do
  end do

  return
end
subroutine r8gb_print ( m, n, ml, mu, a, title )

!*****************************************************************************80
!
!! R8GB_PRINT prints an R8GB matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1..
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  integer ( kind = 4 ) m
  character ( len = * ) title

  call r8gb_print_some ( m, n, ml, mu, a, 1, 1, m, n, title )

  return
end
subroutine r8gb_print_some ( m, n, ml, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8GB_PRINT_SOME prints some of an R8GB matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1..
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) m
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - mu - ml )
    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi + ml )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i < j - ml - mu  .or. j + ml < i ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+ml+mu+1,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8gb_random ( m, n, ml, mu, seed, a )

!*****************************************************************************80
!
!! R8GB_RANDOM randomizes an R8GB matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine assumes it is setting
!    up an unfactored matrix, so it only uses the first MU upper bands,
!    and does not place nonzero values in the fillin bands.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row
  integer ( kind = 4 ) seed

  do j = 1, n
    do row = 1, 2*ml+mu+1
      i = row + j - ml - mu - 1
      if ( ml < row .and. 1 <= i .and. i <= m ) then
        a(row,j) = r8_uniform_01 ( seed )
      else
        a(row,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8gb_sl ( n, ml, mu, a_lu, pivot, b, job )

!*****************************************************************************80
!
!! R8GB_SL solves a system factored by R8GB_FA.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A_LU(2*ML+MU+1,N), the LU factors from R8GB_FA.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from R8GB_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution.
!
!    Input, integer ( kind = 4 ) JOB.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(2*ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  real ( kind = 8 ) t

  m = mu + ml + 1
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    if ( 1 <= ml ) then

      do k = 1, n-1

        lm = min ( ml, n-k )
        l = pivot(k)

        if ( l /= k ) then
          t    = b(l)
          b(l) = b(k)
          b(k) = t
        end if

        b(k+1:k+lm) = b(k+1:k+lm) + b(k) * a_lu(m+1:m+lm,k)

      end do
    end if
!
!  Solve U * X = Y.
!
    do k = n, 1, -1

      b(k) = b(k) / a_lu(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm

      b(lb:lb+lm-1) = b(lb:lb+lm-1) - b(k) * a_lu(la:la+lm-1,k)

    end do
!
!  Solve A' * X = B.
!
  else
!
!  Solve U' * Y = B.
!
    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      b(k) = ( b(k) - sum ( a_lu(la:la+lm-1,k) * b(lb:lb+lm-1) ) ) &
        / a_lu(m,k)
    end do
!
!  Solve L' * X = Y.
!
    if ( 1 <= ml ) then

      do k = n-1, 1, -1

        lm = min ( ml, n-k )
        b(k) = b(k) + sum ( a_lu(m+1:m+lm,k) * b(k+1:k+lm) )
        l = pivot(k)

        if ( l /= k ) then
          t    = b(l)
          b(l) = b(k)
          b(k) = t
        end if

      end do

    end if

  end if

  return
end
subroutine r8gb_to_r8s3 ( m, n, ml, mu, a, nz_num, isym, row, col, b )

!*****************************************************************************80
!
!! R8GB_TO_R8S3 copies an R8GB matrix to an R8S3 matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will copy nonzero
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!    This routine reorders the entries of A so that the first N entries
!    are exactly the diagonal entries of the matrix, in order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrices.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrices.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths of A1.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries in A.
!    This number can be obtained by calling R8GB_NZ_NUM.
!
!    Output, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric,
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!    For this routine, ISYM is always output 0.
!
!    Output, integer ( kind = 4 ) ROW(NZ_NUM), the row indices.
!
!    Output, integer ( kind = 4 ) COL(NZ_NUM), the column indices.
!
!    Output, real ( kind = 8 ) B(NZ_NUM), the R8S3 matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(2*ml+mu+1,n)
  real ( kind = 8 ) b(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) row(nz_num)

  isym = 0
  nz = 0

  do i = 1, m
    do j = 1, n
      if ( i - ml <= j .and. j <= i + mu + ml ) then
        if ( a(ml+mu+1+i-j,j) /= 0.0D+00 ) then

          if ( nz_num <= nz ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8GB_TO_R8S3 - Fatal error!'
            write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
            write ( *, '(a)' ) '  But the matrix has more nonzeros than that!'
            stop
          end if

          nz = nz + 1
          row(nz) = i
          col(nz) = j
          b(nz) = a(ml+mu+1+i-j,j)

        end if
      end if
    end do
  end do

  if ( nz < nz_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GB_TO_R8S3 - Warning!'
    write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
    write ( *, '(a,i8)' ) '  But the number of nonzeros is ', nz
  end if

  return
end
subroutine r8gb_to_r8sp ( m, n, ml, mu, a, nz_num, row, col, b )

!*****************************************************************************80
!
!! R8GB_TO_R8SP copies an R8GB matrix to an R8SP matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will copy nonzero
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrices.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrices.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths of A1.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries in A.
!    This number can be obtained by calling R8GB_NZ_NUM.
!
!    Output, integer ( kind = 4 ) ROW(NZ_NUM), the row indices.
!
!    Output, integer ( kind = 4 ) COL(NZ_NUM), the column indices.
!
!    Output, real ( kind = 8 ) B(NZ_NUM), the R8SP matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(2*ml+mu+1,n)
  real ( kind = 8 ) b(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) row(nz_num)

  nz = 0

  do i = 1, m

    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu + ml )

    do j = jlo, jhi

      if ( a(ml+mu+1+i-j,j) == 0.0D+00 ) then
        cycle
      end if

      if ( nz_num <= nz ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8GB_TO_R8SP - Fatal error!'
        write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
        write ( *, '(a)' ) '  But the matrix has more nonzeros than that!'
        stop
      end if

      nz = nz + 1
      row(nz) = i
      col(nz) = j
      b(nz) = a(ml+mu+1+i-j,j)

    end do
  end do

  if ( nz < nz_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GB_TO_R8SP - Warning!'
    write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
    write ( *, '(a,i8)' ) '  But the number of nonzeros is ', nz
  end if

  return
end
subroutine r8gb_to_r8vec ( m, n, ml, mu, a, x )

!*****************************************************************************80
!
!! R8GB_TO_R8VEC copies an R8GB matrix to an R8VEC.
!
!  Discussion:
!
!    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
!    a data item carries its dimensionality implicitly, and so cannot be
!    regarded sometimes as a vector and sometimes as an array.
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the array.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the array to be copied.
!
!    Output, real ( kind = 8 ) X((2*ML+MU+1)*N), the vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  real ( kind = 8 ) x((2*ml+mu+1)*n)

  do j = 1, n

    ihi = min ( ml + mu, ml + mu + 1 - j )
    do i = 1, ihi
      x(i+(j-1)*(2*ml+mu+1)) = 0.0D+00
    end do

    ilo = max ( ihi + 1, 1 )
    ihi = min ( 2*ml+mu+1, ml+mu+m+1-j )
    do i = ilo, ihi
      x(i+(j-1)*(2*ml+mu+1)) = a(i,j)
    end do

    ilo = ihi + 1
    ihi = 2*ml+mu+1
    do i = ilo, ihi
      x(i+(j-1)*(2*ml+mu+1)) = 0.0D+00
    end do

  end do

  return
end
subroutine r8gb_to_r8ge ( m, n, ml, mu, a, b )

!*****************************************************************************80
!
!! R8GB_TO_R8GE copies an R8GB matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will copy nonzero
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrices.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrices.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths of A1.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Output, real ( kind = 8 ) B(M,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  real ( kind = 8 ) b(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m
    do j = 1, n
      if ( i - ml <= j .and. j <= i + mu + ml ) then
        b(i,j) = a(ml+mu+1+i-j,j)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8gb_trf ( m, n, ml, mu, a, pivot, info )

!*****************************************************************************80
!
!! R8GB_TRF performs a LAPACK-style PLU factorization of an R8GB matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    This is a simplified, standalone version of the LAPACK
!    routine R8GBTRF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    Original FORTRAN77 version by the LAPACK group.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.  
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.  
!    0 <= N.
!
!    Input, integer ( kind = 4 ) ML, the number of subdiagonals within the
!    band of A.  0 <= ML.
!
!    Input, integer ( kind = 4 ) MU, the number of superdiagonals within 
!    the band of A.  0 <= MU.
!
!    Input/output, real ( kind = 8 ) A(2*ML+MU+1,N).  On input, the matrix A
!    in band storage, and on output, information about the PLU factorization.
!
!    Output, integer ( kind = 4 ) PIVOT(min(M,N)), the pivot indices;
!    for 1 <= i <= min(M,N), row i of the matrix was interchanged with
!    row IPIV(i).
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    = 0: successful exit;
!    < 0: an input argument was illegal;
!    > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!         has been completed, but the factor U is exactly
!         singular, and division by zero will occur if it is used
!         to solve a system of equations.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km
  integer ( kind = 4 ) kv
  real ( kind = 8 ) piv
!
  info = 0
!
!  KV is the number of superdiagonals in the factor U, allowing for fill-in.
!
  kv = mu + ml
!
!  Set fill-in elements in columns MU+2 to KV to zero.
!
  do j = mu + 2, min ( kv, n )
    do i = kv - j + 2, ml
      a(i,j) = 0.0D+00
    end do
  end do
!
!  JU is the index of the last column affected by the current stage
!  of the factorization.
!
  ju = 1

  do j = 1, min ( m, n )
!
!  Set the fill-in elements in column J+KV to zero.
!
    if ( j + kv <= n ) then
      a(1:ml,j+kv) = 0.0D+00
    end if
!
!  Find the pivot and test for singularity.
!  KM is the number of subdiagonal elements in the current column.
!
    km = min ( ml, m-j )

    piv = abs ( a(kv+1,j) )
    jp = kv + 1

    do i = kv + 2, kv + km + 1
      if ( piv < abs ( a(i,j) ) ) then
        piv = abs ( a(i,j) )
        jp = i
      end if
    end do

    jp = jp - kv

    pivot(j) = jp + j - 1

    if ( a(kv+jp,j) /= 0.0D+00 ) then

      ju = max ( ju, min ( j+mu+jp-1, n ) )
!
!  Apply interchange to columns J to JU.
!
      if ( jp /= 1 ) then

        do i = 0, ju - j
          call r8_swap ( a(kv+jp-i,j+i), a(kv+1-i,j+i) )
        end do

      end if
!
!  Compute the multipliers.
!
      if ( 0 < km ) then

        a(kv+2:kv+km+1,j) = a(kv+2:kv+km+1,j) / a(kv+1,j)
!
!  Update the trailing submatrix within the band.
!
        if ( j < ju ) then

          do k = 1, ju-j

            if ( a(kv+1-k,j+k) /= 0.0D+00 ) then

              do i = 1, km
                a(kv+i+1-k,j+k) = a(kv+i+1-k,j+k) - a(kv+i+1,j) * a(kv+1-k,j+k)
              end do

            end if

          end do

        end if

      end if

    else
!
!  If pivot is zero, set INFO to the index of the pivot
!  unless a zero pivot has already been found.
!
      if ( info == 0 ) then
        info = j
      end if

    end if

  end do

  return
end
subroutine r8gb_trs ( n, ml, mu, nrhs, trans, a, pivot, b, info )

!*****************************************************************************80
!
!! R8GB_TRS solves an R8GB linear system factored by R8GB_TRF.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 1999
!
!  Author:
!
!    Original FORTRAN77 version by the LAPACK group.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, the number of subdiagonals within the 
!    band of A.  ML must be at least 0, and no greater than N - 1.
!
!    Input, integer ( kind = 4 ) MU, the number of superdiagonals within the 
!    band of A.  MU must be at least 0, and no greater than N - 1.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides and the 
!    number of columns of the matrix B.  NRHS must be positive.
!
!    Input, character TRANS, specifies the form of the system.
!    'N':  A * x = b  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the LU factorization of the 
!    band matrix A, computed by R8GB_TRF.  
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot indices; for 1 <= I <= N, 
!    row I of the matrix was interchanged with row PIVOT(I).
!
!    Input/output, real ( kind = 8 ) B(N,NRHS),
!    On entry, the right hand side vectors B.
!    On exit, the solution vectors, X.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    = 0:  successful exit
!    < 0: if INFO = -K, the K-th argument had an illegal value
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrhs

  real ( kind = 8 ) a(2*ml+mu+1,n)
  real ( kind = 8 ) b(n,nrhs)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kd
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lm
  real ( kind = 8 ) temp
  character trans
!
!  Test the input parameters.
!
  info = 0

  if ( trans /= 'N' .and. trans /= 'n' .and. &
       trans /= 'T' .and. trans /= 't' .and. &
       trans /= 'C' .and. trans /= 'c' ) then
    info = -1
  else if ( n <= 0 ) then
    info = -2
  else if ( ml < 0 ) then
    info = -3
  else if ( mu < 0 ) then
    info = -4
  else if ( nrhs <= 0 ) then
    info = -5
  end if

  if ( info /= 0 ) then
    return
  end if

  kd = mu + ml + 1
!
!  Solve A * x = b.
!
!  Solve L * x = b, overwriting b with x.
!
!  L is represented as a product of permutations and unit lower
!  triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!  where each transformation L(i) is a rank-one modification of
!  the identity matrix.
!
  if ( trans == 'N' .or. trans == 'n' ) then

    if ( 0 < ml ) then

      do j = 1, n - 1

        lm = min ( ml, n-j )
        l = pivot(j)

        do i = 1, nrhs
          call r8_swap ( b(l,i), b(j,i) )
        end do

        do k = 1, nrhs
          if ( b(j,k) /= 0.0D+00 ) then
            b(j+1:j+lm,k) = b(j+1:j+lm,k) - a(kd+1:kd+lm,j) * b(j,k)
          end if
        end do

      end do

    end if
!
!  Solve U * x = b, overwriting b with x.
!
    do i = 1, nrhs

      do j = n, 1, -1
        if ( b(j,i) /= 0.0D+00 ) then
          l = ml + mu + 1 - j
          b(j,i) = b(j,i) / a(ml+mu+1,j)
          do k = j - 1, max ( 1, j - ml - mu ), -1
            b(k,i) = b(k,i) - a(l+k,j) * b(j,i)
          end do
        end if
      end do

    end do

  else
!
!  Solve A' * x = b.
!
!  Solve U' * x = b, overwriting b with x.
!
    do i = 1, nrhs

      do j = 1, n
        temp = b(j,i)
        l = ml + mu + 1 - j
        do k = max ( 1, j - ml - mu ), j - 1
          temp = temp - a(l+k,j) * b(k,i)
        end do
        b(j,i) = temp / a(ml+mu+1,j)
      end do

    end do
!
!  Solve L' * x = b, overwriting b with x.
!
    if ( 0 < ml ) then

      do j = n - 1, 1, -1

        lm = min ( ml, n-j )

        do k = 1, nrhs
          b(j,k) = b(j,k) - sum ( b(j+1:j+lm,k) * a(kd+1:kd+lm,j) )
        end do

        l = pivot(j)

        do i = 1, nrhs
          call r8_swap ( b(l,i), b(j,i) )
        end do

      end do

    end if

  end if

  return
end
subroutine r8gb_vxm ( m, n, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8GB_VXM multiplies an R8VEC by an R8GB matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    LINPACK and LAPACK storage of general band matrices requires
!    an extra ML upper diagonals for possible fill in entries during
!    Gauss elimination.  This routine does not access any entries
!    in the fill in diagonals, because it assumes that the matrix
!    has NOT had Gauss elimination applied to it.  If the matrix
!    has been Gauss eliminated, then the routine R8GB_MU must be
!    used instead.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product X*A.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m)

  do j = 1, n
    b(j) = 0.0D+00
    ilo = max ( 1, j - mu )
    ihi = min ( m, j + ml )
    do i = ilo, ihi
      b(j) = b(j) + x(i) * a(i-j+ml+mu+1,j)
    end do
  end do

  return
end
subroutine r8gd_error ( n, ndiag, ierror )

!*****************************************************************************80
!
!! R8GD_ERROR checks the dimensions of an R8GD matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 1998
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Output, integer ( kind = 4 ) IERROR, reports whether any errors were 
!    detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if N is illegal;
!    IERROR = IERROR + 2 if NDIAG is illegal.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  ierror = 0

  if ( n < 1 ) then
    ierror = ierror + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'R8GD_ERROR - Illegal N = ', n
  end if

  if ( ndiag < 1 .or. 2 * n - 1 < ndiag ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'R8GD_ERROR - Illegal NDIAG = ', ndiag
  end if

  return
end
subroutine r8gd_indicator ( n, ndiag, offset, a )

!*****************************************************************************80
!
!! R8GD_INDICATOR sets up an R8GD indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Output, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) offset(ndiag)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        a(i,diag) = real ( fac * i + j, kind = 8 )
      else
        a(i,diag) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8gd_mxv ( n, ndiag, offset, a, x, b )

!*****************************************************************************80
!
!! R8GD_MXV multiplies an R8GD matrix by an R8VEC.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Example:
!
!    The "offset" value is printed near the first entry of each diagonal
!    in the original matrix, and above the columns in the new matrix.
!
!    Original matrix               New Matrix
!
!      0    1   2   3   4   5        -3  -2   0   1   3   5
!       \    \   \   \   \   \
!        11  12   0  14   0  16      --  --  11  12  14  16
!   -1 =  0  22  23   0  25   0      --  --  22  23  25  --
!   -2 = 31   0  33  34   0  36      --  31  33  34  36  --
!   -3 = 41  42   0  44  45   0      41  42  44  45  --  --
!   -4 =  0  52  53   0  55  56      52  53  55  56  --  --
!   -5 =  0   0  63  64  65  66      63  64  66  --  --  --
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal 
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) offset(ndiag)
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        b(i) = b(i) + a(i,diag) * x(j)
      end if
    end do
  end do

  return
end
subroutine r8gd_print ( n, ndiag, offset, a, title )

!*****************************************************************************80
!
!! R8GD_PRINT prints an R8GD matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the 
!    diagonal storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) offset(ndiag)
  character ( len = * ) title

  call r8gd_print_some ( n, ndiag, offset, a, 1, 1, n, n, title )

  return
end
subroutine r8gd_print_some ( n, ndiag, offset, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8GD_PRINT_SOME prints some of an R8GD matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal 
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) off
  integer ( kind = 4 ) offset(ndiag)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0D+00
        off = j - i
        do diag = 1, ndiag
          if ( off == offset(diag) ) then
            aij = a(i,diag)
          end if
        end do

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8gd_random ( n, ndiag, offset, seed, a )

!*****************************************************************************80
!
!! R8GD_RANDOM randomizes an R8GD matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal 
!    storage.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) offset(ndiag)
  integer ( kind = 4 ) seed

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        a(i,diag) = r8_uniform_01 ( seed )
      else
        a(i,diag) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8gd_to_r8ge ( n, ndiag, offset, a, b )

!*****************************************************************************80
!
!! R8GD_TO_R8GE copies an R8GD matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal 
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) offset(ndiag)

  b(1:n,1:n) = 0.0D+00

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        b(i,j) = a(i,diag)
      end if
    end do
  end do

  return
end
subroutine r8gd_vxm ( n, ndiag, offset, a, x, b )

!*****************************************************************************80
!
!! R8GD_VXM multiplies an R8VEC by an R8GD matrix.
!
!  Discussion:
!
!    The R8GD storage format is suitable for matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0.
!    Each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
!
!    Now, assuming that only a few of these diagonals contain nonzeros,
!    then for the I-th diagonal to be saved, we stored its offset in
!    OFFSET(I), and its entries in column I of the matrix.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal 
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8GD matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product X*A.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) offset(ndiag)
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        b(j) = b(j) + x(i) * a(i,diag)
      end if
    end do
  end do

  return
end
subroutine r8ge_co ( n, a, pivot, rcond, z )

!*****************************************************************************80
!
!! R8GE_CO factors an R8GE matrix and estimates its condition number.
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
!    EPSILON/RCOND.
!
!    If RCOND is so small that the logical expression
!      1.0D+00 + rcond == 1.0D+00
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate
!    underflows.
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
!    Output, real ( kind = 8 ) RCOND, an estimate of the reciprocal
!    condition number of A.
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
  real ( kind = 8 ) ek
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) pivot(n)
  real ( kind = 8 ) rcond
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
!  RCOND = 1 / ( norm(A) * (estimate of norm(inverse(A))) )
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

    if ( k+1 <= n ) then

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

  if ( anorm /= 0.0D+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0D+00
  end if

  return
end
subroutine r8ge_det ( n, a_lu, pivot, det )

!*****************************************************************************80
!
!! R8GE_DET computes the determinant of a matrix factored by R8GE_FA or R8GE_TRF.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 March 2003
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
!    Input, real ( kind = 8 ) A_LU(N,N), the LU factors from R8GE_FA or R8GE_TRF.
!
!    Input, integer ( kind = 4 ) PIVOT(N), as computed by R8GE_FA or R8GE_TRF.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pivot(n)

  det = 1.0D+00

  do i = 1, n
    det = det * a_lu(i,i)
    if ( pivot(i) /= i ) then
      det = - det
    end if
  end do

  return
end
subroutine r8ge_dilu ( m, n, a, d )

!*****************************************************************************80
!
!! R8GE_DILU produces the diagonal incomplete LU factor of an R8GE matrix.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
!    Output, real ( kind = 8 ) D(M), the DILU factor.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) d(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m
    if ( i <= n ) then
      d(i) = a(i,i)
    else
      d(i) = 0.0D+00
    end if
  end do

  do i = 1, min ( m, n )
    d(i) = 1.0D+00 / d(i)
    do j = i + 1, min ( m, n )
      d(j) = d(j) - a(j,i) * d(i) * a(i,j)
    end do
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
subroutine r8ge_fs ( n, a, b, info )

!*****************************************************************************80
!
!! R8GE_FS factors and solves an R8GE system.
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
!    R8GE_FS does not save the LU factors of the matrix, and hence cannot
!    be used to efficiently solve multiple linear systems, or even to
!    factor A at one time, and solve a single linear system at a later time.
!
!    R8GE_FS uses partial pivoting, but no pivot vector is required.
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
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B is the right hand side of the linear system.
!    On output, B is the solution of the linear system.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipiv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  real ( kind = 8 ) piv
  real ( kind = 8 ) row(n)
  real ( kind = 8 ) temp

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a(jcol,jcol) )
    ipiv = jcol
    do i = jcol + 1, n
      if ( piv < abs ( a(i,jcol) ) ) then
        piv = abs ( a(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0D+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_FS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a(jcol,1:n)
      a(jcol,1:n) = a(ipiv,1:n)
      a(ipiv,1:n) = row(1:n)

      temp = b(jcol)
      b(jcol) = b(ipiv)
      b(ipiv) = temp

    end if
!
!  Scale the pivot row.
!
    a(jcol,jcol+1:n) = a(jcol,jcol+1:n) / a(jcol,jcol)
    b(jcol) = b(jcol) / a(jcol,jcol)
    a(jcol,jcol) = 1.0D+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a(i,jcol) /= 0.0D+00 ) then
        temp = - a(i,jcol)
        a(i,jcol) = 0.0D+00
        a(i,jcol+1:n) = a(i,jcol+1:n) + temp * a(jcol,jcol+1:n)
        b(i) = b(i) + temp * b(jcol)
      end if
    end do

  end do
!
!  Back solve.
!
  do jcol = n, 2, -1
    b(1:jcol-1) = b(1:jcol-1) - a(1:jcol-1,jcol) * b(jcol)
  end do

  return
end
subroutine r8ge_fss ( n, a, nb, b, info )

!*****************************************************************************80
!
!! R8GE_FSS factors and solves multiple R8GE systems.
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
!    This routine does not save the LU factors of the matrix, and hence cannot
!    be used to efficiently solve multiple linear systems, or even to
!    factor A at one time, and solve a single linear system at a later time.
!
!    This routine uses partial pivoting, but no pivot vector is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 June 2009
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
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input/output, real ( kind = 8 ) B(N,NB).
!    On input, the right hand sides of the linear system.
!    On output, the solutions of the linear systems.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipiv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  real ( kind = 8 ) piv
  real ( kind = 8 ) row(n)
  real ( kind = 8 ) t(nb)
  real ( kind = 8 ) temp

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a(jcol,jcol) )
    ipiv = jcol
    do i = jcol + 1, n
      if ( piv < abs ( a(i,jcol) ) ) then
        piv = abs ( a(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0D+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_FSS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a(jcol,1:n)
      a(jcol,1:n) = a(ipiv,1:n)
      a(ipiv,1:n) = row(1:n)

      t(1:nb)      = b(jcol,1:nb)
      b(jcol,1:nb) = b(ipiv,1:nb)
      b(ipiv,1:nb) = t(1:nb)

    end if
!
!  Scale the pivot row.
!
    a(jcol,jcol+1:n) = a(jcol,jcol+1:n) / a(jcol,jcol)
    b(jcol,1:nb) = b(jcol,1:nb) / a(jcol,jcol)
    a(jcol,jcol) = 1.0D+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a(i,jcol) /= 0.0D+00 ) then
        temp = - a(i,jcol)
        a(i,jcol) = 0.0D+00
        a(i,jcol+1:n) = a(i,jcol+1:n) + temp * a(jcol,jcol+1:n)
        b(i,1:nb) = b(i,1:nb) + temp * b(jcol,1:nb)
      end if
    end do

  end do
!
!  Back solve.
!
  do j = 1, nb
    do jcol = n, 2, -1
      b(1:jcol-1,j) = b(1:jcol-1,j) - a(1:jcol-1,jcol) * b(jcol,j)
    end do
  end do

  return
end
subroutine r8ge_identity ( n, a )

!*****************************************************************************80
!
!! R8GE_IDENTITY copies the identity matrix to an R8GE matrix.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) A(N,N), the N by N identity matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ), parameter :: one = 1.0D+00

  a(1:n,1:n) = 0.0D+00

  call r8mat_diag_set_scalar ( n, a, one )

  return
end
subroutine r8ge_ilu ( m, n, a, l, u )

!*****************************************************************************80
!
!! R8GE_ILU produces the incomplete LU factors of an R8GE matrix.
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
!    The incomplete LU factors of the M by N matrix A are:
!
!      L, an M by M unit lower triangular matrix,
!      U, an M by N upper triangular matrix
!
!    with the property that L and U are computed in the same way as
!    the usual LU factors, except that, whenever an off diagonal element
!    of the original matrix is zero, then the corresponding value of
!    U is forced to be zero.
!
!    This condition means that it is no longer the case that A = L*U.
!
!    On the other hand, L and U will have a simple sparsity structure
!    related to that of A.  The incomplete LU factorization is generally
!    used as a preconditioner in iterative schemes applied to sparse
!    matrices.  It is presented here merely for illustration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
!    Output, real ( kind = 8 ) L(M,M), the M by M unit lower triangular factor.
!
!    Output, real ( kind = 8 ) U(M,N), the M by N upper triangular factor.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) l(m,m)
  real ( kind = 8 ) u(m,n)
!
!  Initialize:
!
!    L := M by M Identity
!    U := A
!
  call r8ge_identity ( m, l )

  u(1:m,1:n) = a(1:m,1:n)

  do j = 1, min ( m - 1, n )
!
!  Zero out the entries in column J, from row J+1 to M.
!
    do i = j + 1, m

      if ( u(i,j) /= 0.0D+00 ) then

        l(i,j) = u(i,j) / u(j,j)
        u(i,j) = 0.0D+00

        do k = j+1, n
          if ( u(i,k) /= 0.0D+00 ) then
            u(i,k) = u(i,k) - l(i,j) * u(j,k)
          end if
        end do

      end if

    end do

  end do

  return
end
subroutine r8ge_indicator ( m, n, a )

!*****************************************************************************80
!
!! R8GE_INDICATOR sets up an R8GE indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, m
    do j = 1, n
      a(i,j) = real ( fac * i + j, kind = 8 )
    end do
  end do

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
subroutine r8ge_ml ( n, a_lu, pivot, x, b, job )

!*****************************************************************************80
!
!! R8GE_ML computes A * x or A' * x, using R8GE_FA factors.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!    It is assumed that R8GE_FA has overwritten the original matrix
!    information by LU factors.  R8GE_ML is able to reconstruct the
!    original matrix from the LU factor data.
!
!    R8GE_ML allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 1998
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
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector computed by R8GE_FA.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(N), the result of the multiplication.
!
!    Input, integer ( kind = 4 ) JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      b(1:j-1) = b(1:j-1) + a_lu(1:j-1,j) * b(j)
      b(j) = a_lu(j,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n-1, 1, -1

      b(j+1:n) = b(j+1:n) - a_lu(j+1:n,j) * b(j)

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

    end do

  else
!
!  Y = (PL)' * X:
!
    do j = 1, n - 1

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

      b(j) = b(j) - sum ( b(j+1:n) * a_lu(j+1:n,j) )

    end do
!
!  B = U' * Y = ( PL * U )' * X = A' * X.
!
    do i = n, 1, -1
      b(i+1:n) = b(i+1:n) + b(i) * a_lu(i,i+1:n)
      b(i) = b(i) * a_lu(i,i)
    end do

  end if

  return
end
subroutine r8ge_mu ( m, n, a_lu, trans, pivot, x, b )

!*****************************************************************************80
!
!! R8GE_MU computes A * x or A' * x, using R8GE_TRF factors.
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
!    It is assumed that R8GE_TRF has overwritten the original matrix
!    information by PLU factors.  R8GE_MU is able to reconstruct the
!    original matrix from the PLU factor data.
!
!    R8GE_MU allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix.
!
!    Input, real ( kind = 8 ) A_LU(M,N), the LU factors from R8GE_TRF.
!
!    Input, character TRANS, specifies the form of the system of equations:
!    'N':  A * x = b  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, integer ( kind = 4 ) PIVOT(*), the pivot vector computed by R8GE_TRF.
!
!    Input, real ( kind = 8 ) X(*), the vector to be multiplied.
!    For the untransposed case, X should have N entries.
!    For the transposed case, X should have M entries.
!
!    Output, real ( kind = 8 ) B(*), the result of the multiplication.
!    For the untransposed case, B should have M entries.
!    For the transposed case, B should have N entries.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(m,n)
  real ( kind = 8 ) b(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pivot(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn_max
  integer ( kind = 4 ) npiv
  real ( kind = 8 ) t
  character trans
  real ( kind = 8 ) x(*)
  real ( kind = 8 ), allocatable, dimension ( : ) :: y

  npiv = min ( m - 1, n )
  mn_max = max ( m, n )
  allocate ( y(1:mn_max) )

  if ( trans == 'n' .or. trans == 'N' ) then
!
!  Y[MN] = U[MNxN] * X[N].
!
    y(1:n) = 0.0D+00

    do j = 1, n

      do i = 1, min ( j, m )
        y(i) = y(i) + a_lu(i,j) * x(j)
      end do

    end do
!
!  Z[M] = L[MxMN] * Y[MN] = L[MxMN] * U[MNxN] * X[N].
!
    do i = 1, m

      if ( i <= n ) then
        b(i) = y(i)
      else
        b(i) = 0.0D+00
      end if

    end do

    do j = min ( m-1, n ), 1, -1
      b(j+1:m) = b(j+1:m) + a_lu(j+1:m,j) * y(j)
    end do
!
!  B = P * Z = P * L * Y = P * L * U * X = A * x.
!
    do j = npiv, 1, -1

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

    end do

  else if ( trans == 't' .or. trans == 'T' .or. &
            trans == 'c' .or. trans == 'C' ) then
!
!  Y = tranpose(P) * X:
!
    do i = 1, npiv

      k = pivot(i)

      if ( k /= i ) then
        t    = x(k)
        x(k) = x(i)
        x(i) = t
      end if

    end do

    do i = 1, n

      if ( i <= m ) then
        b(i) = x(i)
      else
        b(i) = 0.0D+00
      end if

    end do
!
!  Z = tranpose(L) * Y:
!
    do j = 1, min ( m - 1, n )
      b(j) = b(j) + sum ( x(j+1:m) * a_lu(j+1:m,j) )
    end do
!
!  B = U' * Z.
!
    do i = m, 1, -1
      b(i+1:n) = b(i+1:n) + b(i) * a_lu(i,i+1:n)
      if ( i <= n ) then
        b(i) = b(i) * a_lu(i,i)
      end if
    end do
!
!  Now restore X.
!
    do i = npiv, 1, -1

      k = pivot(i)

      if ( k /= i ) then
        t    = x(k)
        x(k) = x(i)
        x(i) = t
      end if

    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_MU - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of TRANS = ' // trans
    stop

  end if

  deallocate ( y )

  return
end
subroutine r8ge_mxm ( n1, n2, n3, a, b, c )

!*****************************************************************************80
!
!! R8GE_MXM multiplies two R8GE matrices.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the order of the matrices.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N1,N2), B(N2,N3), the R8GE factor matrices.
!
!    Output, real ( kind = 8 ) C(N1,N3), the R8GE product matrix.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  real ( kind = 8 ) a(n1,n2)
  real ( kind = 8 ) b(n2,n3)
  real ( kind = 8 ) c(n1,n3)

  c(1:n1,1:n3) = matmul ( a(1:n1,1:n2), b(1:n2,1:n3) )

  return
end
subroutine r8ge_mxv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8GE_MXV multiplies an R8GE matrix by an R8VEC.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(M), the product A * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) x(n)

  b(1:m) = matmul ( a(1:m,1:n), x(1:n) )

  return
end
subroutine r8ge_np_det ( n, a_lu, det )

!*****************************************************************************80
!
!! R8GE_NP_DET computes the determinant of a matrix factored by R8GE_NP_FA.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 March 2004
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
!    Input, real ( kind = 8 ) A_LU(N,N), the LU factors from R8GE_NP_FA.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i

  det = 1.0D+00

  do i = 1, n
    det = det * a_lu(i,i)
  end do

  return
end
subroutine r8ge_np_fa ( n, a, info )

!*****************************************************************************80
!
!! R8GE_NP_FA factors an R8GE matrix by nonpivoting Gaussian elimination.
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
!    R8GE_NP_FA is a version of the LINPACK routine SGEFA, but uses no
!    pivoting.  It will fail if the matrix is singular, or if any zero
!    pivot is encountered.
!
!    If R8GE_NP_FA successfully factors the matrix, R8GE_NP_SL may be called
!    to solve linear systems involving the matrix.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, A contains the matrix to be factored.
!    On output, A contains information about the factorization,
!    which must be passed unchanged to R8GE_NP_SL for solutions.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  info = 0

  do k = 1, n-1

    if ( a(k,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_NP_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if

    a(k+1:n,k) = - a(k+1:n,k) / a(k,k)
    do j = k+1, n
      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)
    end do

  end do

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop
  end if

  return
end
subroutine r8ge_np_inverse ( n, a )

!*****************************************************************************80
!
!! R8GE_NP_INVERSE computes the inverse of a matrix factored by R8GE_NP_FA.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 1999
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
!    On input, the factor information computed by R8GE_NP_FA.
!    On output, the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) temp
  real ( kind = 8 ) work(n)
!
!  Compute Inverse(U).
!
  do k = 1, n

    a(k,k) = 1.0D+00 / a(k,k)
    a(1:k-1,k) = - a(1:k-1,k) * a(k,k)

    do j = k + 1, n

      temp = a(k,j)
      a(k,j) = 0.0D+00
      a(1:k,j) = a(1:k,j) + temp * a(1:k,k)

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

  end do

  return
end
subroutine r8ge_np_ml ( n, a_lu, x, b, job )

!*****************************************************************************80
!
!! R8GE_NP_ML computes A * x or x * A, for a matrix factored by R8GE_NP_FA.
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
!    The matrix A is assumed to have been factored by R8GE_NP_FA.
!
!    R8GE_NP_ML allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 1998
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
!    Input, real ( kind = 8 ) A_LU(N,N), the LU factors from R8GE_NP_FA.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(N), the result of the multiplication.
!
!    Input, integer ( kind = 4 ) JOB, determines the multiplication to
!    be carried out:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real ( kind = 8 ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Compute U * X = Y:
!
    do i = 1, n
      b(i) = sum ( a_lu(i,i:n) * b(i:n) )
    end do
!
!  Compute L * Y = B:
!
    do j = n-1, 1, -1
      b(j+1:n) = b(j+1:n) - a_lu(j+1:n,j) * b(j)
    end do

  else
!
!  Compute L' * X = Y:
!
    do j = 1, n-1
      b(j) = b(j) - sum ( b(j+1:n) * a_lu(j+1:n,j) )
    end do
!
!  Compute U' * Y = B:
!
    do j = n, 1, -1
      b(j) = sum ( b(1:j) * a_lu(1:j,j) )
    end do

  end if

  return
end
subroutine r8ge_np_sl ( n, a_lu, b, job )

!*****************************************************************************80
!
!! R8GE_NP_SL solves a system factored by R8GE_NP_FA.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 1998
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
!    Input, real ( kind = 8 ) A_LU(N,N), the LU factors from R8GE_NP_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!
!    On input, B contains the right hand side vector B.
!    On output, B contains the solution X.
!
!    Input, integer ( kind = 4 ) JOB.
!    If JOB is zero, the routine will solve A * x = b.
!    If JOB is nonzero, the routine will solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
!
!  Solve A * x = b.
!
  if ( job == 0 ) then

    do k = 1, n-1
      b(k+1:n) = b(k+1:n) + a_lu(k+1:n,k) * b(k)
    end do

    do k = n, 1, -1
      b(k) = b(k) / a_lu(k,k)
      b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
    end do
!
!  Solve A' * X = B.
!
  else

    do k = 1, n
      b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
    end do

    do k = n-1, 1, -1
      b(k) = b(k) + sum ( b(k+1:n) * a_lu(k+1:n,k) )
    end do

  end if

  return
end
subroutine r8ge_np_trf ( m, n, a, info )

!*****************************************************************************80
!
!! R8GE_NP_TRF computes the LU factorization of an R8GE matrix.
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
!    R8GE_NP_TRF is a nonpivoting version of R8GE_TRF, and will fail if
!    a zero element is encountered along the diagonal.
!
!    The factorization has the form
!      A = L * U
!    where L is lower triangular with unit diagonal elements (lower
!    trapezoidal if N < M), and U is upper triangular (upper trapezoidal
!    if M < N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.  
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.  
!    0 <= N.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On entry, the M by N matrix to be factored.
!    On exit, the factors L and U from the factorization
!    A = L*U; the unit diagonal elements of L are not stored.
!
!    Output, integer ( kind = 4 ) INFO.
!    = 0: successful exit
!    < 0: if INFO = -K, the K-th argument had an illegal value
!    > 0: if INFO = K, U(K,K) is exactly zero. The factorization
!         has been completed, but the factor U is exactly
!         singular, and division by zero will occur if it is used
!         to solve a system of equations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
!
!  Test the input parameters.
!
  info = 0

  if ( m < 0 ) then
    info = - 1
    return
  else if ( n < 0 ) then
    info = - 2
    return
  end if

  if ( m == 0 .or. n == 0 ) then
    return
  end if

  do j = 1, min ( m, n )
!
!  Compute elements J+1:M of the J-th column.
!
    if ( a(j,j) /= 0.0D+00 ) then
      a(j+1:m,j) = a(j+1:m,j) / a(j,j)
    else if ( info == 0 ) then
      info = j
    end if
!
!  Update the trailing submatrix.
!
    if ( j < min ( m, n ) ) then

      do ii = j+1, m
        a(ii,j+1:n) = a(ii,j+1:n) - a(ii,j) * a(j,j+1:n)
      end do

    end if

  end do

  return
end
subroutine r8ge_np_trm ( m, n, a, x, b, job )

!*****************************************************************************80
!
!! R8GE_NP_TRM computes A * x or A' * x, for a matrix factored by R8GE_NP_TRF.
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
!    The matrix A is assumed to have been factored by R8GE_NP_TRF.
!
!    R8GE_NP_TRM allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the matrix.  M and N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix factors computed
!    by R8GE_NP_TRF.
!
!    Input, real ( kind = 8 ) X(*), the vector to be multiplied.
!    If JOB is 0, X must have dimension N.
!    If JOB is nonzero, X must have dimension M.
!
!    Output, real ( kind = 8 ) B(*), the result of the multiplication.
!    If JOB is 0, B must have dimension M.
!    If JOB is nonzero, B must have dimension N.
!
!    Input, integer ( kind = 4 ) JOB, determines the multiplication to
!    be carried out:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  real ( kind = 8 ) x(*)

  if ( job == 0 ) then

    b(1:m) = 0.0D+00
!
!  Compute U * X = Y:
!
    do i = 1, m
      b(i) = sum ( a(i,i:n) * x(i:n) )
    end do
!
!  Compute L * Y = B:
!
    do i = min ( m, n + 1 ), 2, -1
      b(i) = b(i) + sum ( a(i,1:i-1) * b(1:i-1) )
    end do
  else

    b(1:n) = 0.0D+00
!
!  Compute L' * X = Y:
!
    do i = 1, min ( m, n )
      b(i) = x(i) + sum ( a(i+1:m,i) * x(i+1:m) )
    end do
!
!  Compute U' * Y = B:
!
    do i = min ( m, n ), 1, -1
      b(i) = sum ( a(1:i,i) * b(1:i) )
    end do

  end if

  return
end
subroutine r8ge_np_trs ( n, nrhs, trans, a, b, info )

!*****************************************************************************80
!
!! R8GE_NP_TRS solves a system of linear equations factored by R8GE_NP_TRF.
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
!    R8GE_NP_TRS is a nonpivoting version of R8GE_TRS.
!
!    R8GE_TRS solves a system of linear equations
!      A * x = b  or  A' * X = B
!    with a general N by N matrix A using the LU factorization computed
!    by R8GE_NP_TRF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.  0 <= N.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.  
!    0 <= NRHS.
!
!    Input, character TRANS, pecifies the form of the system of equations:
!    'N':  A * x = b  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, real ( kind = 8 ) A(N,N), the factors L and U from the factorization
!    A = L*U as computed by R8GE_NP_TRF.
!
!    Input/output, real ( kind = 8 ) B(N,NRHS).
!    On entry, the right hand side matrix B.
!    On exit, the solution matrix X.
!
!    Output, integer ( kind = 4 ) INFO
!    = 0:  successful exit
!    < 0:  if INFO = -I, the I-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrhs

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,nrhs)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character trans

  info = 0

  if ( trans /= 'n' .and. trans /= 'N' .and. &
       trans /= 't' .and. trans /= 'T' .and. &
       trans /= 'c' .and. trans /= 'C' ) then
    info = - 1
    return
  else if ( n < 0 ) then
    info = - 2
    return
  else if ( nrhs < 0 ) then
    info = - 3
    return
  end if

  if ( n == 0 .or. nrhs == 0 ) then
    return
  end if

  if ( trans == 'n' .or. trans == 'N' ) then
!
!  Solve L * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = 1, n - 1
        b(j+1:n,k) = b(j+1:n,k) - a(j+1:n,j) * b(j,k)
      end do
    end do
!
!  Solve U * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = n, 1, -1
        b(j,k) = b(j,k) / a(j,j)
        b(1:j-1,k) = b(1:j-1,k) - a(1:j-1,j) * b(j,k)
      end do
    end do

  else
!
!  Solve U' * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = 1, n
        b(j,k) = b(j,k) / a(j,j)
        b(j+1:n,k) = b(j+1:n,k) - a(j,j+1:n) * b(j,k)
      end do
    end do
!
!  Solve L' * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = n, 2, -1
        b(1:j-1,k) = b(1:j-1,k) - a(j,1:j-1) * b(j,k)
      end do
    end do

  end if

  return
end
subroutine r8ge_plu ( m, n, a, p, l, u )

!*****************************************************************************80
!
!! R8GE_PLU produces the PLU factors of an R8GE matrix.
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
!    The PLU factors of the M by N matrix A are:
!
!      P, an M by M permutation matrix P,
!      L, an M by M unit lower triangular matrix,
!      U, an M by N upper triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
!    Output, real ( kind = 8 ) P(M,M), the M by M permutation factor.
!
!    Output, real ( kind = 8 ) L(M,M), the M by M unit lower triangular factor.
!
!    Output, real ( kind = 8 ) U(M,N), the M by N upper triangular factor.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) l(m,m)
  real ( kind = 8 ) p(m,m)
  integer ( kind = 4 ) pivot_row
  real ( kind = 8 ) pivot_value
  real ( kind = 8 ) u(m,n)
!
!  Initialize:
!
!    P: = M by M Identity
!    L: = M by M Identity
!    U: = A
!
  call r8ge_identity ( m, l )
  call r8ge_identity ( m, p )

  u(1:m,1:n) = a(1:m,1:n)
!
!  On step J, find the pivot row and the pivot value.
!
  do j = 1, min ( m-1, n )

    pivot_value = 0.0D+00
    pivot_row = -1

    do i = j, m

      if ( pivot_value < abs ( u(i,j) ) ) then
        pivot_value = abs ( u(i,j) )
        pivot_row = i
      end if

    end do
!
!  If the pivot row is nonzero, swap rows J and PIVOT_ROW.
!
    if ( pivot_row /= -1 ) then

      call r8row_swap ( m, n, u, j, pivot_row )

      call r8row_swap ( m, m, l, j, pivot_row )
      call r8col_swap ( m, m, l, j, pivot_row )

      call r8col_swap ( m, m, p, j, pivot_row )
!
!  Zero out the entries in column J, from row J+1 to M.
!
      do i = j + 1, m

        if ( u(i,j) /= 0.0D+00 ) then

          l(i,j) = u(i,j) / u(j,j)
          u(i,j) = 0.0D+00
          u(i,j+1:n) = u(i,j+1:n) - l(i,j) * u(j,j+1:n)

        end if

      end do

    end if

  end do

  return
end
subroutine r8ge_poly ( n, a, p )

!*****************************************************************************80
!
!! R8GE_POLY computes the characteristic polynomial of an R8GE matrix.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 1998
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
!    Input, real ( kind = 8 ) A(N,N), the R8GE matrix.
!
!    Output, real ( kind = 8 ) P(0:N), the coefficients of the characteristic
!    polynomial of A.  P(I) contains the coefficient of X**I.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) order
  real ( kind = 8 ) p(0:n)
  real ( kind = 8 ) trace
  real ( kind = 8 ) work1(n,n)
  real ( kind = 8 ) work2(n,n)
!
!  Initialize WORK1 to the identity matrix.
!
  call r8ge_identity ( n, work1 )

  p(n) = 1.0D+00

  do order = n-1, 0, -1
!
!  Work2 = A * WORK1.
!
    work2(1:n,1:n) = matmul ( a(1:n,1:n), work1(1:n,1:n) )
!
!  Take the trace.
!
    call r8mat_diag_get_vector ( n, work2, diag )

    trace = sum ( diag(1:n) )
!
!  P(ORDER) = - Trace ( WORK2 ) / ( N - ORDER )
!
    p(order) = - trace / real ( n - order, kind = 8 )
!
!  WORK1 := WORK2 + P(ORDER) * Identity.
!
    work1(1:n,1:n) = work2(1:n,1:n)

    call r8mat_diag_add_scalar ( n, work1, p(order) )

  end do

  return
end
subroutine r8ge_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8GE_PRINT prints an R8GE matrix.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8ge_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8ge_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8GE_PRINT_SOME prints some of an R8GE matrix.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( r8_is_int ( a(i,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8ge_random ( m, n, seed, a )

!*****************************************************************************80
!
!! R8GE_RANDOM randomizes an R8GE matrix.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) seed

  call r8mat_uniform_01 ( m, n, seed, a )

  return
end
subroutine r8ge_res ( m, n, a, x, b, r )

!*****************************************************************************80
!
!! R8GE_RES computes the residual vector for an R8GE system.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.  M and N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the original, UNFACTORED R8GE matrix.
!
!    Input, real ( kind = 8 ) X(N), the estimated solution.
!
!    Input, real ( kind = 8 ) B(M), the right hand side vector.
!
!    Output, real ( kind = 8 ) R(M), the residual vector, b - A * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) r(m)
  real ( kind = 8 ) x(n)

  r(1:m) = b(1:m) - matmul ( a(1:m,1:n), x(1:n) )

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
subroutine r8ge_sl_it ( n, a, a_lu, pivot, b, job, x, r )

!*****************************************************************************80
!
!! R8GE_SL_IT applies one step of iterative refinement following R8GE_SL.
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
!    It is assumed that:
!
!    * the original matrix A has been factored by R8GE_FA;
!    * the linear system A * x = b has been solved once by R8GE_SL.
!
!    (Actually, it is not necessary to solve the system once using R8GE_SL.
!    You may simply supply the initial estimated solution X = 0.)
!
!    Each time this routine is called, it will compute the residual in
!    the linear system, apply one step of iterative refinement, and
!    add the computed correction to the current solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 May 1999
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
!    Input, real ( kind = 8 ) A(N,N), the original, UNFACTORED R8GE matrix.
!
!    Input, real ( kind = 8 ) A_LU(N,N), the LU factors from R8GE_FA.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from R8GE_FA.
!
!    Input, real ( kind = 8 ) B(N), the right hand side vector.
!
!    Input, integer ( kind = 4 ) JOB, specifies the operation.
!    0, solve A*X=B.
!    nonzero, solve A'*X=B.
!
!    Input/output, real ( kind = 8 ) X(N), an estimate of the solution 
!    of A * x = b.  On output, the solution has been improved by one 
!    step of iterative refinement.
!
!    Output, real ( kind = 8 ) R(N), contains the correction terms added to X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) job
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) x(n)
!
!  Compute the residual vector.
!
  call r8ge_res ( n, n, a, x, b, r )
!
!  Solve A * dx = r
!
  call r8ge_sl ( n, a_lu, pivot, r, job )
!
!  Add dx to x.
!
  x(1:n) = x(1:n) + r(1:n)

  return
end
subroutine r8ge_to_r8gb ( m, n, ml, mu, a, b )

!*****************************************************************************80
!
!! R8GE_TO_R8GB copies an R8GE matrix to an R8GB matrix.
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
!    It usually doesn't make sense to try to store a general matrix
!    in a band matrix format.  You can always do it, but it will take
!    more space, unless the general matrix is actually banded.
!
!    The purpose of this routine is to allow a user to set up a
!    banded matrix in the easy-to-use general format, and have this
!    routine take care of the compression of the data into general
!    format.  All the user has to do is specify the bandwidths.
!
!    Note that this routine "believes" what the user says about the
!    bandwidth.  It will assume that all entries in the general matrix
!    outside of the bandwidth are zero.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.
!
!    The R8GB storage format is for an M by N banded matrix, with lower
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
!    extra superdiagonals, which may be required to store nonzero entries
!    generated during Gaussian elimination.
!
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will copy nonzero
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrices.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrices.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths of 
!    the matrix.  ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
!    Output, real ( kind = 8 ) B(2*ML+MU+1,N), the R8GB matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(2*ml+mu+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo

  b(1:2*ml+mu+1,1:n) = 0.0D+00

  do i = 1, m
    jlo = max ( i - ml, 1 )
    jhi = min ( i + mu, n )
    do j = jlo, jhi
      b(ml+mu+1+i-j,j) = a(i,j)
    end do
  end do

  return
end
subroutine r8ge_to_r8vec ( m, n, a, x )

!*****************************************************************************80
!
!! R8GE_TO_R8VEC copies an R8GE matrix to an R8VEC.
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
!    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
!    a data item carries its dimensionality implicitly, and so cannot be
!    regarded sometimes as a vector and sometimes as an array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the array.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be copied.
!
!    Output, real ( kind = 8 ) X(M*N), the vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(m*n)

  k = 0
  do j = 1, n
    do i = 1, m
      k = k + 1
      x(k) = a(i,j)
    end do
  end do

  return
end
subroutine r8ge_trf ( m, n, a, pivot, info )

!*****************************************************************************80
!
!! R8GE_TRF performs a LAPACK-style PLU factorization of an R8GE matrix.
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
!    R8GE_TRF is a standalone version of the LAPACK routine SGETRF.
!
!    The factorization uses partial pivoting with row interchanges,
!    and has the form
!      A = P * L * U
!    where P is a permutation matrix, L is lower triangular with unit
!    diagonal elements (lower trapezoidal if N < M), and U is upper
!    triangular (upper trapezoidal if M < N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 December 1998
!
!  Author:
!
!    Original FORTRAN77 version by the LAPACK group.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.  
!    0 <= M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.  
!    0 <= N.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On entry, the M by N matrix to be factored.
!    On exit, the factors L and U from the factorization
!    A = P*L*U; the unit diagonal elements of L are not stored.
!
!    Output, integer ( kind = 4 ) PIVOT(min(M,N)), the pivot indices;
!    for 1 <= I <= min(M,N), row i of the matrix was interchanged with
!    row PIVOT(I).
!
!    Output, integer ( kind = 4 ) INFO.
!    = 0: successful exit
!    < 0: if INFO = -K, the K-th argument had an illegal value
!    > 0: if INFO = K, U(K,K) is exactly zero. The factorization
!         has been completed, but the factor U is exactly
!         singular, and division by zero will occur if it is used
!         to solve a system of equations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jp
  real ( kind = 8 ) t
!
!  Test the input parameters.
!
  info = 0

  if ( m < 0 ) then
    info = - 1
    return
  else if ( n < 0 ) then
    info = - 2
    return
  end if

  if ( m == 0 .or. n == 0 ) then
    return
  end if

  do j = 1, min ( m, n )
!
!  Find the pivot.
!
    t = abs ( a(j,j) )
    jp = j
    do i = j + 1, m
      if ( t < abs ( a(i,j) ) ) then
        t = abs ( a(i,j) )
        jp = i
      end if
    end do

    pivot(j) = jp
!
!  Apply the interchange to columns 1:N.
!  Compute elements J+1:M of the J-th column.
!
    if ( a(jp,j) /= 0.0D+00 ) then

      if ( jp /= j ) then
        do jj = 1, n
          t        = a(j,jj)
          a(j,jj)  = a(jp,jj)
          a(jp,jj) = t
        end do
      end if

      if ( j < m ) then
        a(j+1:m,j) = a(j+1:m,j) / a(j,j)
      end if

    else if ( info == 0 ) then

      info = j

    end if
!
!  Update the trailing submatrix.
!
    if ( j < min ( m, n ) ) then

      do ii = j+1, m
        a(ii,j+1:n) = a(ii,j+1:n) - a(ii,j) * a(j,j+1:n)
      end do

    end if

  end do

  return
end
subroutine r8ge_trs ( n, nrhs, trans, a, pivot, b, info )

!*****************************************************************************80
!
!! R8GE_TRS solves a system of linear equations factored by R8GE_TRF.
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
!    R8GE_TRS is a standalone version of the LAPACK routine SGETRS.
!
!    R8GE_TRS solves a system of linear equations
!      A * x = b  or  A' * X = B
!    with a general N by N matrix A using the PLU factorization computed
!    by R8GE_TRF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 December 1998
!
!  Author:
!
!    Original FORTRAN77 version by the LAPACK group.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.  0 <= N.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.  
!    0 <= NRHS.
!
!    Input, character TRANS, specifies the form of the system of equations:
!    'N':  A * x = b  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, real ( kind = 8 ) A(N,N), the factors L and U from the factorization
!    A = P*L*U as computed by R8GE_TRF.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot indices from R8GE_TRF.
!
!    Input/output, real ( kind = 8 ) B(N,NRHS).
!    On entry, the right hand side matrix B.
!    On exit, the solution matrix X.
!
!    Output, integer ( kind = 4 ) INFO
!    = 0:  successful exit
!    < 0:  if INFO = -I, the I-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrhs

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,nrhs)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character              trans

  info = 0

  if ( trans /= 'n' .and. trans /= 'N' .and. &
       trans /= 't' .and. trans /= 'T' .and. &
       trans /= 'c' .and. trans /= 'C' ) then
    info = - 1
    return
  else if ( n < 0 ) then
    info = - 2
    return
  else if ( nrhs < 0 ) then
    info = - 3
    return
  end if

  if ( n == 0 .or. nrhs == 0 ) then
    return
  end if

  if ( trans == 'n' .or. trans == 'N' ) then
!
!  Apply row interchanges to the right hand sides.
!
    do i = 1, n
      if ( pivot(i) /= i ) then
        do k = 1, nrhs
          call r8_swap ( b(i,k), b(pivot(i),k) )
        end do
      end if
    end do
!
!  Solve L * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = 1, n - 1
        b(j+1:n,k) = b(j+1:n,k) - a(j+1:n,j) * b(j,k)
      end do
    end do
!
!  Solve U * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = n, 1, -1
        b(j,k) = b(j,k) / a(j,j)
        b(1:j-1,k) = b(1:j-1,k) - a(1:j-1,j) * b(j,k)
      end do
    end do

  else
!
!  Solve U' * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = 1, n
        b(j,k) = b(j,k) / a(j,j)
        b(j+1:n,k) = b(j+1:n,k) - a(j,j+1:n) * b(j,k)
      end do
    end do
!
!  Solve L' * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = n, 2, -1
        b(1:j-1,k) = b(1:j-1,k) - a(j,1:j-1) * b(j,k)
      end do
    end do
!
!  Apply row interchanges to the solution vectors.
!
    do i = n, 1, -1
      if ( pivot(i) /= i ) then
        do k = 1, nrhs
          call r8_swap ( b(i,k), b(pivot(i),k) )
        end do
      end if
    end do

  end if

  return
end
subroutine r8ge_vxm ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8GE_VXM multiplies an R8VEC by an R8GE matrix.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) x(m)

  b(1:n) = matmul ( transpose ( a(1:m,1:n) ), x(1:m) )

  return
end
subroutine r8lt_det ( n, a, det )

!*****************************************************************************80
!
!! R8LT_DET computes the determinant of an R8LT matrix.
!
!  Discussion:
!
!    The R8LT storage format is used for an M by N lower triangular matrix,
!    and sets aside storage even for the entries that must be zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
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
!    Input, real ( kind = 8 ) A(N,N), the R8LT matrix.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det
  real ( kind = 8 ) diag(n)

  call r8mat_diag_get_vector ( n, a, diag )

  det = product ( diag(1:n) )

  return
end
subroutine r8lt_indicator ( m, n, a )

!*****************************************************************************80
!
!! R8LT_INDICATOR sets up an R8LT indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8LT storage format is used for an M by N lower triangular matrix,
!    and sets aside storage even for the entries that must be zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    1 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.  M and N must be positive.
!
!    Output, real ( kind = 8 ) A(M,N), the R8LT matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, m
    do j = 1, min ( i, n )
      a(i,j) = real ( fac * i + j, kind = 8 )
    end do
    do j = i+1, n
      a(i,j) = 0.0D+00
    end do
  end do

  return
end
subroutine r8lt_inverse ( n, a )

!*****************************************************************************80
!
!! R8LT_INVERSE computes the inverse of an R8LT matrix.
!
!  Discussion:
!
!    The R8LT storage format is used for an M by N lower triangular matrix,
!    and sets aside storage even for the entries that must be zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Second edition,
!    Academic Press, 1978,
!    ISBN 0-12-519260-6
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!
!    On input, the lower triangular matrix to be inverted.
!    On output, the inverse of the lower triangular matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

!  Check.
!
  do i = 1, n
    if ( a(i,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8LT_INVERSE - Fatal error!'
      write ( *, '(a)' ) '  Zero diagonal element.'
      stop
    end if
  end do

  do j = 1, n

    do i = 1, n

      if ( i < j ) then

        a(i,j) = 0.0D+00

      else if ( i == j ) then

        a(i,j) = 1.0D+00 / a(i,j)

      else if ( j < i ) then

        a(i,j) = - sum ( a(i,j:i-1) * a(j:i-1,j) ) / a(i,i)

      end if

    end do
  end do

  return
end
subroutine r8lt_mxm ( n, a, b, c )

!*****************************************************************************80
!
!! R8LT_MXM multiplies two R8LT matrices.
!
!  Discussion:
!
!    The R8LT storage format is used for an M by N lower triangular matrix,
!    and sets aside storage even for the entries that must be zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N,N), B(N,N), the R8LT factor matrices.
!
!    Output, real ( kind = 8 ) C(N,N), the R8LT product matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  c(1:n,1:n) = 0.0D+00

  do i = 1, n
    do j = 1, i
      c(i,j) = sum ( a(i,j:i) * b(j:i,j) )
    end do
  end do

  return
end
subroutine r8lt_mxv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8LT_MXV multiplies an R8LT matrix by an R8VEC.
!
!  Discussion:
!
!    The R8LT storage format is used for an M by N lower triangular matrix,
!    and sets aside storage even for the entries that must be zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8LT matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(M), the product A * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jmax
  real ( kind = 8 ) x(n)

  do i = 1, m
    jmax = min ( i, n )
    b(i) = sum ( a(i,1:jmax) * x(1:jmax) )
  end do

  return
end
subroutine r8lt_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8LT_PRINT prints an R8LT matrix.
!
!  Discussion:
!
!    The R8LT storage format is used for an M by N lower triangular matrix,
!    and sets aside storage even for the entries that must be zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8LT matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8lt_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8lt_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8LT_PRINT_SOME prints some of an R8LT matrix.
!
!  Discussion:
!
!    The R8LT storage format is used for an M by N lower triangular matrix,
!    and sets aside storage even for the entries that must be zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8LT matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( ilo < jlo ) then
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i < j ) then
          ctemp(j2) = '              '
        else if ( r8_is_int ( a(i,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8lt_random ( m, n, seed, a )

!*****************************************************************************80
!
!! R8LT_RANDOM randomizes an R8LT matrix.
!
!  Discussion:
!
!    The R8LT storage format is used for an M by N lower triangular matrix,
!    and sets aside storage even for the entries that must be zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix.  M and N must be positive.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(M,N), the R8LT matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  do i = 1, m
    do j = 1, min ( i, n )
      a(i,j) = r8_uniform_01 ( seed )
    end do
    do j = i + 1, n
      a(i,j) = 0.0D+00
    end do
  end do

  return
end
subroutine r8lt_sl ( n, a, b, job )

!*****************************************************************************80
!
!! R8LT_SL solves an R8LT system.
!
!  Discussion:
!
!    The R8LT storage format is used for an M by N lower triangular matrix,
!    and sets aside storage even for the entries that must be zero.
!
!    No factorization of the lower triangular matrix is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the R8LT matrix.
!
!    Input/output, real ( kind = 8 ) B(N).
!
!    On input, the right hand side.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB, is 0 to solve the untransposed system,
!    nonzero to solve the transposed system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job

  if ( job == 0 ) then

    do j = 1, n
      b(j) = b(j) / a(j,j)
      b(j+1:n) = b(j+1:n) - a(j+1:n,j) * b(j)
    end do

  else

    do j = n, 1, -1
      b(j) = b(j) / a(j,j)
      b(1:j-1) = b(1:j-1) - a(j,1:j-1) * b(j)
    end do

  end if

  return
end
subroutine r8lt_vxm ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8LT_VXM multiplies an R8VEC by an R8LT matrix.
!
!  Discussion:
!
!    The R8LT storage format is used for an M by N lower triangular matrix,
!    and sets aside storage even for the entries that must be zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8LT matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(m)

  do i = 1, n
    b(i) = sum ( x(i:m) * a(i:m,i) )
  end do

  return
end
subroutine r8mat_diag_add_scalar ( n, a, s )

!*****************************************************************************80
!
!! R8MAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of a matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N), the N by N matrix to be modified.
!
!    Input, real ( kind = 8 ) S, the value to be added to the diagonal 
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) s

  do i = 1, n
    a(i,i) = a(i,i) + s
  end do

  return
end
subroutine r8mat_diag_get_vector ( n, a, v )

!*****************************************************************************80
!
!! R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of a matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) V(N), the diagonal entries of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) v(n)

  do i = 1, n
    v(i) = a(i,i)
  end do

  return
end
subroutine r8mat_diag_set_scalar ( n, a, s )

!*****************************************************************************80
!
!! R8MAT_DIAG_SET_SCALAR sets the diagonal of a matrix to a scalar value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N), the N by N matrix to be modified.
!
!    Input, real ( kind = 8 ) S, the value to be assigned to the 
!    diagonal of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) s

  do i = 1, n
    a(i,i) = s
  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints a real matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of a real matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)') j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( r8_is_int ( a(i,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8's.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8ncf_indicator ( m, n, nz_num, rowcol, a )

!*****************************************************************************80
!
!! R8NCF_INDICATOR sets up an R8NCF indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) ROWCOL(2,NZ_NUM), the coordinates of 
!    the nonzero entries.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the indicator matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rowcol(2,nz_num)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do k = 1, nz_num
    i = rowcol(1,k)
    j = rowcol(2,k)
    a(k) = real ( fac * i + j, kind = 8 )
  end do

  return
end
subroutine r8ncf_print ( m, n, nz_num, rowcol, a, title )

!*****************************************************************************80
!
!! R8NCF_PRINT prints an R8NCF matrix.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ROWCOL(2,NZ_NUM), the row and column indices
!    of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
! 
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rowcol(2,nz_num)
  character ( len = * ) title

  call r8ncf_print_some ( m, n, nz_num, rowcol, a, 1, 1, m, n, title )

  return
end
subroutine r8ncf_print_some ( m, n, nz_num, rowcol, a, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! R8NCF_PRINT_SOME prints some of an R8NCF matrix.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements 
!    in the matrix.
!
!    Input, integer ( kind = 4 ) ROWCOL(2,NZ_NUM), the row and column indices
!    of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  logical nonzero
  integer ( kind = 4 ) rowcol(2,nz_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      nonzero = .false.

      aij = 0.0D+00
      do j2 = 1, inc
        write ( ctemp(j2), '(f8.0,6x)' ) aij
      end do

      do k = 1, nz_num

        if ( &
          i == rowcol(1,k) .and. &
          j2lo <= rowcol(2,k) .and. &
          rowcol(2,k) <= j2hi ) then 

          j2 = rowcol(2,k) - j2lo + 1
          aij = a(k)

          if ( aij == 0.0D+00 ) then
            cycle
          end if

          nonzero = .true.

          if ( r8_is_int ( aij ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) aij
          else
            write ( ctemp(j2), '(g14.6)' ) aij
          end if
        end if

      end do

      if ( nonzero ) then
        write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )
      end if

    end do

  end do

  return
end
subroutine r8pbl_det ( n, mu, a_lu, det )

!*****************************************************************************80
!
!! R8PBL_DET computes the determinant of a matrix factored by R8PBL_FA.
!
!  Discussion:
!
!    The R8PBL storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and lower triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row 1 of the array.
!    The first subdiagonal in row 2, columns 1 through MU.
!    The second subdiagonal in row 3, columns 1 through MU-1.
!    The MU-th subdiagonal in row MU+1, columns 1 through 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 October 2003
!
!  Author:
!
!    Original FORTRAN77 version Dongarra, Bunch, Moler, Stewart.
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A_LU(MU+1,N), the LU factors from R8PBL_FA.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(mu+1,n)
  real ( kind = 8 ) det

  det = product ( a_lu(1,1:n)**2 )

  return
end
subroutine r8pbl_indicator ( n, mu, a )

!*****************************************************************************80
!
!! R8PBL_INDICATOR sets up an R8PBL indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8PBL storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and lower triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row 1 of the array.
!    The first subdiagonal in row 2, columns 1 through MU.
!    The second subdiagonal in row 3, columns 1 through MU-1.
!    The MU-th subdiagonal in row MU+1, columns 1 through 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2004
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
!    Input, integer ( kind = 4 ) MU, the number of subdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Output, real ( kind = 8 ) A(MU+1,N), the R8PBL matrix.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )
!
!  Zero out the "junk" entries.
!
  do j = n + 1 - mu, n
    do i = n + 1, j + mu
      a(i-j+1,j) = 0.0D+00
    end do
  end do
!
!  Set the meaningful values.
!
  do i = 1, n
    do j = max ( 1, i - mu ), i
      a(i-j+1,j) = real ( fac * i + j, kind = 8 )
    end do
  end do

  return
end
subroutine r8pbl_print ( n, mu, a, title )

!*****************************************************************************80
!
!! R8PBL_PRINT prints an R8PBL matrix.
!
!  Discussion:
!
!    The R8PBL storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and lower triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row 1 of the array.
!    The first subdiagonal in row 2, columns 1 through MU.
!    The second subdiagonal in row 3, columns 1 through MU-1.
!    The MU-th subdiagonal in row MU+1, columns 1 through 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
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
!    Input, integer ( kind = 4 ) MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8PBL matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  character ( len = * ) title

  call r8pbl_print_some ( n, mu, a, 1, 1, n, n, title )

  return
end
subroutine r8pbl_print_some ( n, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8PBL_PRINT_SOME prints some of an R8PBL matrix.
!
!  Discussion:
!
!    The R8PBL storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and lower triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row 1 of the array.
!    The first subdiagonal in row 2, columns 1 through MU.
!    The second subdiagonal in row 3, columns 1 through MU-1.
!    The MU-th subdiagonal in row MU+1, columns 1 through 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
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
!    Input, integer ( kind = 4 ) MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8PBL matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - mu )
    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + mu )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j .and. j <= i + mu ) then
          aij = a(j-i+1,i)
        else if ( j <= i .and. i <= j + mu ) then
          aij = a(i-j+1,j)
        else
          aij = 0.0D+00
        end if

        if ( mu < i-j .or. mu < j-i ) then
          ctemp(j2) = '              '
        else if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8pbl_random ( n, mu, seed, a )

!*****************************************************************************80
!
!! R8PBL_RANDOM randomizes an R8PBL matrix.
!
!  Discussion:
!
!    The R8PBL storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and lower triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row 1 of the array.
!    The first subdiagonal in row 2, columns 1 through MU.
!    The second subdiagonal in row 3, columns 1 through MU-1.
!    The MU-th subdiagonal in row MU+1, columns 1 through 1.
!
!    The matrix returned will be positive definite, but of limited
!    randomness.  The off diagonal elements are random values between
!    0 and 1, and the diagonal element of each row is selected to
!    ensure strict diagonal dominance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
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
!    Input, integer ( kind = 4 ) MU, the number of subdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(MU+1,N), the R8PBL matrix.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sum2
!
!  Zero out the "junk" entries.
!
  do j = n+1-mu, n
    a(2+n-j:mu+1,j) = 0.0D+00
  end do
!
!  Set the off diagonal values.
!
  do i = 1, n
    do j = max ( 1, i - mu ), i - 1
      a(i-j+1,j) = r8_uniform_01 ( seed )
    end do
  end do
!
!  Set the diagonal values.
!
  do i = 1, n

    sum2 = 0.0D+00

    jlo = max ( 1, i - mu )
    do j = jlo, i-1
      sum2 = sum2 + abs ( a(i-j+1,j) )
    end do

    jhi = min ( i + mu, n )
    do j = i+1, jhi
      sum2 = sum2 + abs ( a(j-i+1,i) )
    end do

    r = r8_uniform_01 ( seed )

    a(1,i) = ( 1.0D+00 + r ) * ( sum2 + 0.01D+00 )

  end do

  return
end
subroutine r8pbl_to_r8ge ( n, mu, a, b )

!*****************************************************************************80
!
!! R8PBL_TO_R8GE copies an R8PBL matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8PBL storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and lower triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row 1 of the array.
!    The first subdiagonal in row 2, columns 1 through MU.
!    The second subdiagonal in row 3, columns 1 through MU-1.
!    The MU-th subdiagonal in row MU+1, columns 1 through 1.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth of A.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8PBL matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n
      if ( i <= j .and. j <= i + mu ) then
        b(i,j) = a(j-i+1,i)
      else if ( i - mu <= j .and. j < i ) then
        b(i,j) = a(i-j+1,j)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8pbu_cg ( n, mu, a, b, x )

!*****************************************************************************80
!
!! R8PBU_CG uses the conjugate gradient method on an R8PBU system.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!    The matrix A must be a positive definite symmetric band matrix.
!
!    The method is designed to reach the solution after N computational
!    steps.  However, roundoff may introduce unacceptably large errors for
!    some problems.  In such a case, calling the routine again, using
!    the computed solution as the new starting estimate, should improve
!    the results.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    FS Beckman,
!    The Solution of Linear Equations by the Conjugate Gradient Method,
!    Mathematical Methods for Digital Computers, pages 62-72.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the number of superdiagonals.
!    MU must be at least 0, and no more than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8PBU matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side vector.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, an estimate for the solution, which may be 0.
!    On output, the approximate solution vector.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) ap(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) beta
  integer ( kind = 4 ) it
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) pap
  real ( kind = 8 ) pr
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rap
  real ( kind = 8 ) x(n)
!
!  Initialize
!    AP = A * x,
!    R  = b - A * x,
!    P  = b - A * x.
!
  call r8pbu_mxv ( n, mu, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP=A*P.
!
    call r8pbu_mxv ( n, mu, a, p, ap )
!
!  Compute the dot products
!    PAP = P*AP,
!    PR  = P*R
!  Set
!    ALPHA = PR / PAP.
!
    pap = sum ( p(1:n) * ap(1:n) )
    pr = sum ( p(1:n) * r(1:n) )

    if ( pap == 0.0D+00 ) then
      return
    end if

    alpha = pr / pap
!
!  Set
!    X = X + ALPHA * P
!    R = R - ALPHA * AP.
!
    x(1:n) = x(1:n) + alpha * p(1:n)
    r(1:n) = r(1:n) - alpha * ap(1:n)
!
!  Compute the vector dot product
!    RAP = R*AP
!  Set
!    BETA = - RAP / PAP.
!
    rap = sum ( r(1:n) * ap(1:n) )

    beta = - rap / pap
!
!  Update the perturbation vector
!    P = R + BETA * P.
!
    p(1:n) = r(1:n) + beta * p(1:n)

  end do

  return
end
subroutine r8pbu_det ( n, mu, a_lu, det )

!*****************************************************************************80
!
!! R8PBU_DET computes the determinant of a matrix factored by R8PBU_FA.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 1998
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the number of superdiagonals of the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real ( kind = 8 ) A_LU(MU+1,N), the LU factors from R8PBU_FA.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(mu+1,n)
  real ( kind = 8 ) det

  det = product ( a_lu(mu+1,1:n)**2 )

  return
end
subroutine r8pbu_fa ( n, mu, a, info )

!*****************************************************************************80
!
!! R8PBU_FA factors an R8PBU matrix.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!    The matrix A must be a positive definite symmetric band matrix.
!
!    Once factored, linear systems A*x=b involving the matrix can be solved
!    by calling R8PBU_SL.  No pivoting is performed.  Pivoting is not necessary
!    for positive definite symmetric matrices.  If the matrix is not positive
!    definite, the algorithm may behave correctly, but it is also possible
!    that an illegal divide by zero will occur.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 1998
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt
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
!    Input, integer ( kind = 4 ) MU, the number of superdiagonals of the matrix.
!    MU must be at least 0, and no more than N-1.
!
!    Input/output, real ( kind = 8 ) A(MU+1,N), the N by N matrix, stored 
!    in LINPACK positive definite symmetric band matrix storage.
!    On output, A contains information describing a factored form
!    of the matrix, that can be used to solve linear systems
!    A*x=b, using R8PBU_SL.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, the matrix is nonsingular.
!    nonzero, the matrix is singular.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mm
  real ( kind = 8 ) s

  info = 0

  do j = 1, n

    ik = mu + 1
    jk = max ( j - mu, 1 )
    mm = max ( mu + 2 - j, 1 )

    s = 0.0D+00

    do k = mm, mu

      a(k,j) = ( a(k,j) - sum ( a(ik:ik+k-mm-1,jk) * a(mm:k-1,j) ) ) &
        / a(mu+1,jk)

      s = s + a(k,j) * a(k,j)

      ik = ik - 1
      jk = jk + 1

    end do

    s = a(mu+1,j) - s

    if ( s <= 0.0D+00 ) then
      info = j
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8PBU_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Nonpositive pivot on step ', info
      stop
    end if

    a(mu+1,j) = sqrt ( s )

  end do

  return
end
subroutine r8pbu_indicator ( n, mu, a )

!*****************************************************************************80
!
!! R8PBU_INDICATOR sets up an R8PBU indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2004
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
!    Input, integer ( kind = 4 ) MU, the number of superdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Output, real ( kind = 8 ) A(MU+1,N), the R8PBU matrix.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )
!
!  Zero out the "junk" entries.
!
  do j = 1, mu
    do i = 1, mu + 1 - j
      a(i,j) = 0.0D+00
    end do
  end do
!
!  Set the meaningful values.
!
  do i = 1, n
    do j = i, min ( i+mu, n )
      a(mu+1+i-j,j) = real ( fac * i + j, kind = 8 )
    end do
  end do

  return
end
subroutine r8pbu_ml ( n, mu, a_lu, x, b )

!*****************************************************************************80
!
!! R8PBU_ML multiplies an R8VEC times a matrix that was factored by R8PBU_FA.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 1998
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
!    Input, integer ( kind = 4 ) MU, the number of superdiagonals of the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real ( kind = 8 ) A_LU(MU+1,N), the LU factors from R8PBU_FA.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  b(1:n) = x(1:n)
!
!  Multiply U * X = Y.
!
  do k = 1, n

    ilo = max ( 1, k - mu )
    do i = ilo, k - 1
      b(i) = b(i) + a_lu(mu+1+i-k,k) * b(k)
    end do

    b(k) = a_lu(mu+1,k) * b(k)

  end do
!
!  Multiply L * Y = B.
!
  do k = n, 1, -1

    jhi = min ( k + mu, n )
    do j = k + 1, jhi
      b(j) = b(j) + a_lu(mu+1+k-j,j) * b(k)
    end do

    b(k) = a_lu(mu+1,k) * b(k)

  end do

  return
end
subroutine r8pbu_mxv ( n, mu, a, x, b )

!*****************************************************************************80
!
!! R8PBU_MXV multiplies an R8PBU matrix by an R8VEC.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 October 1998
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
!    Input, integer ( kind = 4 ) MU, the number of superdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8PBU matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the result vector A * x.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ieqn
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)
!
!  Multiply X by the diagonal of the matrix.
!
  b(1:n) = a(mu+1,1:n) * x(1:n)
!
!  Multiply X by the superdiagonals of the matrix.
!
  do i = mu, 1, -1
    do j = mu + 2 - i, n
      ieqn = i + j - mu - 1
      b(ieqn) = b(ieqn) + a(i,j) * x(j)
      b(j) = b(j) + a(i,j) * x(ieqn)
    end do
  end do

  return
end
subroutine r8pbu_print ( n, mu, a, title )

!*****************************************************************************80
!
!! R8PBU_PRINT prints an R8PBU matrix.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
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
!    Input, integer ( kind = 4 ) MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8PBU matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  character ( len = * ) title

  call r8pbu_print_some ( n, mu, a, 1, 1, n, n, title )

  return
end
subroutine r8pbu_print_some ( n, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8PBU_PRINT_SOME prints some of an R8PBU matrix.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
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
!    Input, integer ( kind = 4 ) MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8PBU matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - mu )
    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + mu )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j .and. j <= i + mu ) then
          aij = a(mu+1+i-j,j)
        else if ( i - mu <= j .and. j <= i ) then
          aij = a(mu+1+j-i,i)
        else
          aij = 0.0D+00
        end if

        if ( mu < i-j .or. mu < j-i ) then
          ctemp(j2) = '              '
        else if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8pbu_random ( n, mu, seed, a )

!*****************************************************************************80
!
!! R8PBU_RANDOM randomizes an R8PBU matrix.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!    The matrix returned will be positive definite, but of limited
!    randomness.  The off diagonal elements are random values between
!    0 and 1, and the diagonal element of each row is selected to
!    ensure strict diagonal dominance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2003
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
!    Input, integer ( kind = 4 ) MU, the number of superdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(MU+1,N), the R8PBU matrix.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sum2
!
!  Zero out the "junk" entries.
!
  do j = 1, mu
    a(1:mu+1-j,j) = 0.0D+00
  end do
!
!  Set the off diagonal values.
!
  do i = 1, n
    do j = i+1, min ( i+mu, n )
      a(mu+1+i-j,j) = r8_uniform_01 ( seed )
    end do
  end do
!
!  Set the diagonal values.
!
  do i = 1, n

    sum2 = 0.0D+00

    jlo = max ( 1, i - mu )
    do j = jlo, i-1
      sum2 = sum2 + abs ( a(mu+1+j-i,i) )
    end do

    jhi = min ( i + mu, n )
    do j = i+1, jhi
      sum2 = sum2 + abs ( a(mu+1+i-j,j) )
    end do

    r = r8_uniform_01 ( seed )

    a(mu+1,i) = ( 1.0D+00 + r ) * ( sum2 + 0.01D+00 )

  end do

  return
end
subroutine r8pbu_sl ( n, mu, a_lu, b )

!*****************************************************************************80
!
!! R8PBU_SL solves an R8PBU system factored by R8PBU_FA.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 1998
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the number of superdiagonals of the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real ( kind = 8 ) A_LU(MU+1,N), the LU factors from R8PBU_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B contains the right hand side of the linear system
!    to be solved.
!    On output, B contains X, the solution vector.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(mu+1,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) k
!
!  Solve L * Y = B.
!
  do k = 1, n
    ilo = max ( 1, k - mu )
    b(k) = ( b(k) - sum ( b(ilo:k-1) * a_lu(mu+1+ilo-k:mu,k) ) ) &
      / a_lu(mu+1,k)
  end do
!
!  Solve U * X = Y.
!
  do k = n, 1, -1

    b(k) = b(k) / a_lu(mu+1,k)

    ilo = max ( 1, k - mu )
    do i = ilo, k - 1
      b(i) = b(i) - b(k) * a_lu(mu+1+i-k,k)
    end do

  end do

  return
end
subroutine r8pbu_sor ( n, mu, a, b, eps, itchk, itknt, itmax, omega, x )

!*****************************************************************************80
!
!! R8PBU_SOR uses SOR iteration to solve an R8PBU linear system.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!    The matrix A must be a positive definite symmetric band matrix.
!
!    A relaxation factor OMEGA may be used.
!
!    The iteration will proceed until a convergence test is met,
!    or the iteration limit is reached.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 1998
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
!    Input, integer ( kind = 4 ) MU, the number of superdiagonals in the 
!    matrix.  MU must be at least 0, and no more than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8PBU matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the system.
!
!    Input, real ( kind = 8 ) EPS, convergence tolerance for the system. 
!    The vector b - A * x is computed every ITCHK iterations, and if the 
!    maximum entry of this vector is of norm less than EPS, the program
!    will return.
!
!    Input, integer ( kind = 4 ) ITCHK, the interval between convergence checks.
!    ITCHK steps will be taken before any check is made on whether the iteration
!    has converged.  ITCHK should be at least 1 and no greater
!    than ITMAX.
!
!    Output, integer ( kind = 4 ) ITKNT, the number of iterations taken.
!
!    Input, integer ( kind = 4 ) ITMAX, the maximum number of iterations 
!    allowed.  The program will return to the user if this many iterations 
!    are taken without convergence.
!
!    Input, real ( kind = 8 ) OMEGA, the relaxation factor.  OMEGA must be
!    strictly between 0 and 2.  Use OMEGA = 1 for no relaxation, classical
!    Jacobi iteration.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, a starting vector for the iteration.
!    On output, the current approximation to the solution.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) eps
  real ( kind = 8 ) err
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itchk
  integer ( kind = 4 ) itknt
  integer ( kind = 4 ) itmax
  real ( kind = 8 ) omega
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtemp(n)

  if ( itchk <= 0 .or. itmax < itchk ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8PBU_SOR - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal ITCHK= ', itchk
    stop
  end if

  if ( itmax <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8PBU_SOR - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive ITMAX =', itmax
    stop
  end if

  if ( omega <= 0.0D+00 .or. 2.0D+00 <= omega ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8PBU_SOR - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal value of OMEGA = ', omega
    stop
  end if

  itknt = 0
!
!  Take ITCHK steps of the iteration before doing a convergence check.
!
  do while ( itknt <= itmax )

    do it = 1, itchk
!
!  Compute XTEMP(I) = B(I) + A(I,I) * X(I) - SUM ( J=1 to N ) A(I,J) * X(J).
!
      call r8pbu_mxv ( n, mu, a, x, xtemp )

      xtemp(1:n) = x(1:n) + ( b(1:n) - xtemp(1:n) ) / a(mu+1,1:n)
!
!  Compute the next iterate as a weighted combination of the
!  old iterate and the just computed standard Jacobi iterate.
!
      if ( omega /= 1.0D+00 ) then
        xtemp(1:n) = ( 1.0D+00 - omega ) * x(1:n) + omega * xtemp(1:n)
      end if

      itknt = itknt + 1
!
!  Copy the new result into the old result vector.
!
      x(1:n) = xtemp(1:n)

    end do
!
!  Compute the maximum residual, the greatest entry in the vector
!  RESID(I) = B(I) - A(I,J) * X(J).
!
    call r8pbu_mxv ( n, mu, a, x, xtemp )

    err = maxval ( abs ( b(1:n) - xtemp(1:n) ) )
!
!  Test to see if we can quit because of convergence,
!
    if ( err <= eps ) then
      return
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_SOR - Warning!'
  write ( *, '(a)' ) '  The iteration did not converge.'

  return
end
subroutine r8pbu_to_r8ge ( n, mu, a, b )

!*****************************************************************************80
!
!! R8PBU_TO_R8GE copies an R8PBU matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth of A1.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the R8PBU matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n
      if ( i <= j .and. j <= i+mu ) then
        b(i,j) = a(mu+1+i-j,j)
      else if ( i-mu <= j .and. j < i ) then
        b(i,j) = a(mu+1+j-i,i)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8po_det ( n, a_lu, det )

!*****************************************************************************80
!
!! R8PO_DET computes the determinant of a matrix factored by R8PO_FA.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A_LU(N,N), the LU factors from R8PO_FA.
!
!    Output, real ( kind = 8 ) DET, the determinant of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i

  det = 1.0D+00

  do i = 1, n
    det = det * a_lu(i,i)
  end do

  return
end
subroutine r8po_fa ( n, a, info )

!*****************************************************************************80
!
!! R8PO_FA factors an R8PO matrix.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!    The positive definite symmetric matrix A has a Cholesky factorization
!    of the form:
!
!      A = R' * R
!
!    where R is an upper triangular matrix with positive elements on
!    its diagonal.  This routine overwrites the matrix A with its
!    factor R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2003
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, the matrix in R8PO storage.
!    On output, the Cholesky factor R in R8GE storage.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal return.
!    K, error condition.  The principal minor of order K is not
!    positive definite, and the factorization was not completed.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) s

  do j = 1, n

    do k = 1, j - 1
      a(k,j) = ( a(k,j) - sum ( a(1:k-1,k) * a(1:k-1,j) ) ) / a(k,k)
    end do

    s = a(j,j) - sum ( a(1:j-1,j)**2 )

    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if

    a(j,j) = sqrt ( s )

  end do

  info = 0
!
!  Since the Cholesky factor is stored in R8GE format, be sure to
!  zero out the lower triangle.
!
  do i = 1, n
    do j = 1, i-1
      a(i,j) = 0.0D+00
    end do
  end do

  return
end
subroutine r8po_indicator ( n, a )

!*****************************************************************************80
!
!! R8PO_INDICATOR sets up an R8PO indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix.  N must be positive.
!
!    Output, real ( kind = 8 ) A(N,N), the R8PO matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, n
    do j = 1, i-1
      a(i,j) = 0.0D+00
    end do
    do j = i, n
      a(i,j) = real ( fac * i + j, kind = 8 )
    end do
  end do

  return
end
subroutine r8po_inverse ( n, a )

!*****************************************************************************80
!
!! R8PO_INVERSE computes the inverse of a matrix factored by R8PO_FA.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, the Cholesky factor, in R8GE storage, returned by R8PO_FA.
!    On output, the inverse matrix, in R8PO storage.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
!
!  Compute Inverse ( R ).
!
  do k = 1, n

    a(k,k) = 1.0D+00 / a(k,k)
    a(1:k-1,k) = -a(1:k-1,k) * a(k,k)

    do j = k + 1, n
      t = a(k,j)
      a(k,j) = 0.0D+00
      a(1:k,j) = a(1:k,j) + t * a(1:k,k)
    end do

  end do
!
!  Compute Inverse ( R ) * ( Inverse ( R ) )'.
!
  do j = 1, n

    do k = 1, j - 1
      t = a(k,j)
      a(1:k,k) = a(1:k,k) + t * a(1:k,j)
    end do

    a(1:j,j) = a(1:j,j) * a(j,j)

  end do

  return
end
subroutine r8po_ml ( n, a_lu, x, b )

!*****************************************************************************80
!
!! R8PO_ML computes A * x = b after A has been factored by R8PO_FA.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A_LU(N,N), the Cholesky factor from R8PO_FA.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
!
!  Compute R * x = y.
!
  do i = 1, n
    b(i) = a_lu(i,i) * x(i) + sum ( a_lu(i,i+1:n) * x(i+1:n) )
  end do
!
!  Compute R' * y = b.
!
  do i = n, 1, -1
    b(i) = a_lu(i,i) * b(i) + sum ( b(1:i-1) * a_lu(1:i-1,i) )
  end do

  return
end
subroutine r8po_mxm ( n, a, b, c )

!*****************************************************************************80
!
!! R8PO_MXM multiplies two R8PO matrices.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N,N), B(N,N), the R8PO factor matrices.
!
!    Output, real ( kind = 8 ) C(N,N), the R8PO product matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) aik
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) bkj
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  c(1:n,1:n) = 0.0D+00

  do i = 1, n

    do j = i, n
      do k = 1, n

        if ( i <= k ) then
          aik = a(i,k)
        else
          aik = a(k,i)
        end if

        if ( k <= j ) then
          bkj = b(k,j)
        else
          bkj = b(j,k)
        end if

        c(i,j) = c(i,j) + aik * bkj

      end do
    end do

  end do

  return
end
subroutine r8po_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! R8PO_MXV multiplies an R8PO matrix by an R8VEC.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the R8PO matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  do i = 1, n
    b(i) = 0.0D+00
    do j = 1, i-1
      b(i) = b(i) + a(j,i) * x(j)
    end do
    do j = i, n
      b(i) = b(i) + a(i,j) * x(j)
    end do
  end do

  return
end
subroutine r8po_print ( n, a, title )

!*****************************************************************************80
!
!! R8PO_PRINT prints an R8PO matrix.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the R8PO matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  character ( len = * ) title

  call r8po_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r8po_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8PO_PRINT_SOME prints some of an R8PO matrix.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the R8PO matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(i,j)
        else
          aij = a(j,i)
        end if

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8po_random ( n, seed, a )

!*****************************************************************************80
!
!! R8PO_RANDOM randomizes an R8PO matrix.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!    The matrix computed here is not simply a set of random numbers in
!    the nonzero slots of the R8PO array.  It is also a positive definite
!    matrix.  It is computed by setting a "random" upper triangular
!    Cholesky factor R, and then computing A = R'*R.
!    The randomness is limited by the fact that all the entries of
!    R will be between 0 and 1.  A truly random R is only required
!    to have positive entries on the diagonal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
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
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) A(N,N), the R8PO matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
!
!  Set the whole matrix to zero.
!
  a(1:n,1:n) = 0.0D+00

  do i = n, 1, -1
!
!  Set row I of R.
!
    do j = i, n
      a(i,j) = r8_uniform_01 ( seed )
    end do
!
!  Consider element J of row I, last to first.
!
    do j = n, i, -1
!
!  Add multiples of row I to lower elements of column J.
!
      a(i+1:j,j) = a(i+1:j,j) + a(i,i+1:j) * a(i,j)
!
!  Reset element J.
!
      a(i,j) = a(i,i) * a(i,j)

    end do
  end do

  return
end
subroutine r8po_sl ( n, a_lu, b )

!*****************************************************************************80
!
!! R8PO_SL solves an R8PO system factored by R8PO_FA.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A_LU(N,N), the Cholesky factor from R8PO_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side.
!    On output, the solution vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) k
!
!  Solve R' * y = b.
!
  do k = 1, n
    b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
  end do
!
!  Solve R * x = y.
!
  do k = n, 1, -1
    b(k) = b(k) / a_lu(k,k)
    b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
  end do

  return
end
subroutine r8po_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R8PO_TO_R8GE copies an R8PO matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the R8PO matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n
      if ( i <= j ) then
        b(i,j) = a(i,j)
      else
        b(i,j) = a(j,i)
      end if
    end do
  end do

  return
end
subroutine r8pp_det ( n, a_lu, det )

!*****************************************************************************80
!
!! R8PP_DET computes the determinant of an R8PP matrix factored by R8PP_FA.
!
!  Discussion:
!
!    The R8PP storage format is used for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A_LU((N*(N+1))/2), the LU factors from R8PO_FA.
!
!    Output, real ( kind = 8 ) DET, the determinant of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu((n*(n+1))/2)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k

  det = 1.0D+00

  k = 0
  do i = 1, n
    k = k + i
    det = det * a_lu(k)
  end do

  return
end
subroutine r8pp_fa ( n, a, info )

!*****************************************************************************80
!
!! R8PP_FA factors an R8PP matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 November 2008
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
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A((N*(N+1))/2).  On input, an R8PP matrix A.
!    On output, an upper triangular matrix R, stored in packed form, 
!    so that A = R'*R.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, for normal return.
!    K, if the leading minor of order K is not positive definite.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kj
  integer ( kind = 4 ) kk
  real ( kind = 8 ) s
  real ( kind = 8 ) t

  info = 0
  jj = 0

  do j = 1, n

    s = 0.0D+00
    kj = jj
    kk = 0

    do k = 1, j-1

      kj = kj + 1
      t = a(kj)
      do i = 1, k - 1
        t = t - a(kk+i) * a(jj+i)
      end do
      kk = kk + k
      t = t / a(kk)
      a(kj) = t
      s = s + t * t

    end do

    jj = jj + j
    s = a(jj) - s

    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if

    a(jj) = sqrt ( s )

  end do

  return
end
subroutine r8pp_indicator ( n, a )

!*****************************************************************************80
!
!! R8PP_INDICATOR sets up an R8PP indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2004
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
!    Output, real ( kind = 8 ) A((N*(N+1))/2), the R8PP matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      a(k) = real ( fac * i + j, kind = 8 )
    end do
  end do

  return
end
subroutine r8pp_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! R8PP_MXV multiplies an R8PP matrix by an R8VEC.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), the R8PP matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  do i = 1, n
    b(i) = 0.0D+00
    do j = 1, i-1
      k = j + ( i * ( i - 1 ) ) / 2
      b(i) = b(i) + a(k) * x(j)
    end do
    do j = i, n
      k = i + ( j * ( j - 1 ) ) / 2
      b(i) = b(i) + a(k) * x(j)
    end do
  end do

  return
end
subroutine r8pp_print ( n, a, title )

!*****************************************************************************80
!
!! R8PP_PRINT prints an R8PP matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
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
!    Input, real ( kind = 8 ) A((N*(N+1))/2), the R8PP matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  character ( len = * ) title

  call r8pp_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r8pp_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8PP_PRINT_SOME prints some of an R8PP matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
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
!    Input, real ( kind = 8 ) A((N*(N+1))/2), the R8PP matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(i+(j*(j-1))/2)
        else
          aij = a(j+(i*(i-1))/2)
        end if

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8pp_random ( n, seed, a )

!*****************************************************************************80
!
!! R8PP_RANDOM randomizes an R8PP matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!    The matrix is computed by setting a "random" upper triangular
!    Cholesky factor R, and then computing A = R'*R.
!    The randomness is limited by the fact that all the entries of
!    R will be between 0 and 1.  A truly random R is only required
!    to have positive entries on the diagonal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 October 1998
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
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) A((N*(N+1))/2), the R8PP matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kj
  integer ( kind = 4 ) seed

  a(1:(n*(n+1))/2) = 0.0D+00

  do i = n, 1, -1
!
!  Set row I of R.
!
    do j = i, n
      ij = i + ( j * ( j - 1 ) ) / 2
      a(ij) = r8_uniform_01 ( seed )
    end do
!
!  Consider element J of row I, last to first.
!
    do j = n, i, -1
!
!  Add multiples of row I to lower elements of column J.
!
      ij = i + ( j * ( j - 1 ) ) / 2

      do k = i+1, j
        kj = k + ( j * ( j - 1 ) ) / 2
        ik = i + ( k * ( k - 1 ) ) / 2
        a(kj) = a(kj) + a(ik) * a(ij)
      end do
!
!  Reset element J.
!
      ii = i + ( i * ( i - 1 ) ) / 2
      a(ij) = a(ii) * a(ij)

    end do
  end do

  return
end
subroutine r8pp_sl ( n, a_lu, b )

!*****************************************************************************80
!
!! R8PP_SL solves an R8PP system factored by R8PP_FA.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 November 2008
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
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A_LU((N*(N+1))/2), the LU factors from R8PP_FA.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu((n*(n+1))/2)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  real ( kind = 8 ) t

  kk = 0

  do k = 1, n
    t = 0.0D+00
    do i = 1, k - 1
      t = t + a_lu(kk+i) * b(i)
    end do
    kk = kk + k
    b(k) = ( b(k) - t ) / a_lu(kk)
  end do

  do k = n, 1, -1
    b(k) = b(k) / a_lu(kk)
    kk = kk - k
    t = -b(k)
    do i = 1, k - 1
      b(i) = b(i) + t * a_lu(kk+i)
    end do
  end do

  return
end
subroutine r8pp_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R8PP_TO_R8GE copies an R8PP matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), the R8PP matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n
      if ( i <= j ) then
        b(i,j) = a(i+(j*(j-1))/2)
      else
        b(i,j) = a(j+(i*(i-1))/2)
      end if
    end do
  end do

  return
end
subroutine r8row_swap ( m, n, a, row1, row2 )

!*****************************************************************************80
!
!! R8ROW_SWAP swaps two rows of a table.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Input, integer ( kind = 4 ) ROW1, ROW2, the two rows to swap.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) row1
  integer ( kind = 4 ) row2
  integer ( kind = 4 ) j

  if ( row1 < 1 .or. m < row1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  ROW1 is out of range.'
    stop
  end if

  if ( row2 < 1 .or. m < row2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  ROW2 is out of range.'
    stop
  end if

  if ( row1 == row2 ) then
    return
  end if

  do j = 1, n
    call r8_swap ( a(row1,j), a(row2,j) )
  end do

  return
end
subroutine r8s3_diagonal ( n, nz_num, isym, row, col, a )

!*****************************************************************************80
!
!! R8S3_DIAGONAL reorders a square R8S3 matrix so the diagonal entries are first.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!    This routine reorders the entries of A so that the first N entries
!    are exactly the diagonal entries of the matrix, in order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Input/output, real ( kind = 8 ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) found
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) t

  found = 0

  do k = 1, nz_num

    do while ( row(k) == col(k) )

      if ( row(k) == k ) then
        found = found + 1
        exit
      end if

      i = row(k)

      call i4_swap ( row(i), row(k) )
      call i4_swap ( col(i), col(k) )

      t    = a(i)
      a(i) = a(k)
      a(k) = t
 
      found = found + 1

      if ( n <= found ) then
        exit
      end if
     
    end do

    if ( n <= found ) then
      exit
    end if

  end do

  if ( found < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8S3_DIAGONAL - Warning!'
    write ( *, '(a,i8)' ) '  Number of diagonal entries expected was ', n
    write ( *, '(a,i8)' ) '  Number found was ', found
  end if

  return
end
subroutine r8s3_indicator ( n, nz_num, isym, row, col, a )

!*****************************************************************************80
!
!! R8S3_INDICATOR sets up an R8S3 indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2004
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
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the indicator matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do k = 1, nz_num
    i = row(k)
    j = col(k)
    a(k) = real ( fac * i + j, kind = 8 )
  end do

  return
end
subroutine r8s3_jac_sl ( n, nz_num, isym, row, col, a, b, x, tol, it_max, &
  job, it, diff )

!*****************************************************************************80
!
!! R8S3_JAC_SL solves an R8S3 system using Jacobi iteration.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!    This routine REQUIRES that the matrix be square, that the matrix
!    have nonzero diagonal entries, and that the first N entries of
!    the array A be exactly the diagonal entries of the matrix, in order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = 8 ) X(N), an approximate solution 
!    to the system.
!
!    Input, real ( kind = 8 ) TOL, a tolerance.  If the maximum change in
!    the solution is less than TOL, the iteration is terminated early.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations taken.
!
!    Output, real ( kind = 8 ) DIFF, the maximum change in the solution
!    on the last iteration.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) col(nz_num)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_new(n)
  real ( kind = 8 ) x_norm

  if ( job == 0 ) then

    do it_num = 1, it_max

      it = it_num
!
!  Initialize to right hand side.
!
      x_new(1:n) = b(1:n)
!
!  Subtract off-diagonal terms.
!
      do k = n+1, nz_num
        i = row(k)
        j = col(k)
        x_new(i) = x_new(i) - a(k) * x(j)
        if ( isym == 1 ) then
          x_new(j) = x_new(j) - a(k) * x(i)
        end if
      end do
!
!  Divide by diagonal terms.
!
      x_new(1:n) = x_new(1:n) / a(1:n)
!
!  Measure change:
!
      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( x_new(1:n) - x(1:n) ) )
!
!  Update.
!
      x(1:n) = x_new(1:n)
!
!  Test for early termination.
!
      if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
        exit
      end if

    end do

  else

    do it_num = 1, it_max

      it = it_num
!
!  Initialize to right hand side.
!
      x_new(1:n) = b(1:n)
!
!  Subtract off-diagonal terms.
!
      do k = n+1, nz_num
        i = col(k)
        j = row(k)
        x_new(i) = x_new(i) - a(k) * x(j)
        if ( isym == 1 ) then
          x_new(j) = x_new(j) - a(k) * x(i)
        end if
      end do
!
!  Divide by diagonal terms.
!
      x_new(1:n) = x_new(1:n) / a(1:n)
!
!  Measure change:
!
      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( x_new(1:n) - x(1:n) ) )
!
!  Update.
!
      x(1:n) = x_new(1:n)
!
!  Test for early termination.
!
      if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
        exit
      end if

    end do

  end if

  return
end
subroutine r8s3_mxv ( m, n, nz_num, isym, row, col, a, x, b )

!*****************************************************************************80
!
!! R8S3_MXV multiplies an R8S3 matrix by an R8VEC.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric, and 1
!    if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(M), the product A * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) x(n)

  b(1:m) = 0.0D+00

  do k = 1, nz_num
    i = row(k)
    j = col(k)
    b(i) = b(i) + a(k) * x(j)
  end do
!
!  Handle the symmetric option.
!
  if ( isym == 1 .and. m == n ) then
    do k = 1, nz_num
      i = col(k)
      j = row(k)
      if ( i /= j ) then
        b(i) = b(i) + a(k) * x(j)
      end if
    end do
  end if

  return
end
subroutine r8s3_print ( m, n, nz_num, isym, row, col, a, title )

!*****************************************************************************80
!
!! R8S3_PRINT prints an R8S3 matrix.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  The symmetric case only makes sense
!    if the matrix is also square, that is, M = N.  In this case, only
!    the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) title

  call r8s3_print_some ( m, n, nz_num, isym, row, col, a, 1, 1, m, &
    n, title )

  return
end
subroutine r8s3_print_some ( m, n, nz_num, isym, row, col, a, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! R8S3_PRINT_SOME prints some of an R8S3 matrix.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric, and 1
!    if the matrix is symmetric.  The symmetric case only makes sense
!    if the matrix is also square, that is, M = N.  In this case, only
!    the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) col(nz_num)
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  logical nonzero
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      nonzero = .false.

      aij = 0.0D+00
      do j2 = 1, inc
        write ( ctemp(j2), '(f8.0,6x)' ) aij
      end do

      do k = 1, nz_num

        if ( i == row(k) .and. j2lo <= col(k) .and. col(k) <= j2hi ) then

          j2 = col(k) - j2lo + 1
          aij = a(k)

          if ( aij == 0.0D+00 ) then
            cycle
          end if

          nonzero = .true.

          if ( r8_is_int ( aij ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) aij
          else
            write ( ctemp(j2), '(g14.6)' ) aij
          end if

        else if ( isym == 1 .and. m == n .and. &
          i == col(k) .and. j2lo <= row(k) .and. row(k) <= j2hi ) then

          j2 = row(k) - j2lo + 1
          aij = a(k)

          if ( aij == 0.0D+00 ) then
            cycle
          end if

          nonzero = .true.

          if ( r8_is_int ( aij ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) aij
          else
            write ( ctemp(j2), '(g14.6)' ) aij
          end if

        end if

      end do

      if ( nonzero ) then
        write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )
      end if

    end do

  end do

  return
end
subroutine r8s3_read ( input_file, n, nz_num, row, col, a )

!*****************************************************************************80
!
!! R8S3_READ reads a square R8S3 matrix from a file.
!
!  Discussion:
!
!    This routine needs the value of NZ_NUM, which can be determined
!    by a call to R8S3_READ_SIZE.
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file to be read.
!
!    Unused, integer N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Output, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  character ( len = * ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8S3_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( input_file ) // '".'
    stop
  end if

  do k = 1, nz_num

    read ( input_unit, *, iostat = ios ) row(k), col(k), a(k)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8S3_READ - Fatal error!'
      write ( *, '(a,i8)' ) '  I/O error while reading record ', k
      stop
    end if

  end do

  close ( unit = input_unit )

  return
end
subroutine r8s3_read_size ( input_file, n, nz_num )

!*****************************************************************************80
!
!! R8S3_READ_SIZE reads the size of a square R8S3 matrix from a file.
!
!  Discussion:
!
!    The value of NZ_NUM is simply the number of records in the input file.
!
!    The value of N is determined as the maximum entry in the row and column
!    vectors.
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file to 
!    be read.
!
!    Output, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements 
!    in the matrix.
!
  implicit none

  real ( kind = 8 ) a_k
  integer ( kind = 4 ) col_k
  character ( len = * ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ) row_k

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8S3_READ_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( input_file ) // '".'
    stop
  end if

  nz_num = 0
  n = 0

  do

    read ( input_unit, *, iostat = ios ) row_k, col_k, a_k

    if ( ios /= 0 ) then
      exit
    end if

    nz_num = nz_num + 1
    n = max ( n, row_k )
    n = max ( n, col_k )

  end do

  close ( unit = input_unit )

  return
end
subroutine r8s3_to_r8ge ( m, n, nz_num, isym, row, col, a, b )

!*****************************************************************************80
!
!! R8S3_TO_R8GE copies an R8S3 matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  The symmetric case only makes sense
!    if the matrix is also square, that is, M = N.  In this case, only
!    the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Output, real ( kind = 8 ) B(M,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(m,n)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)

  b(1:m,1:n) = 0.0D+00

  do k = 1, nz_num
    i = row(k)
    j = col(k)
    b(i,j) = b(i,j) + a(k)
    if ( isym == 1 .and. m == n .and. i /= j ) then
      b(j,i) = b(j,i) + a(k)
    end if
  end do

  return
end
subroutine r8s3_vxm ( m, n, nz_num, isym, row, col, a, x, b )

!*****************************************************************************80
!
!! R8S3_VXM multiplies an R8VEC times an R8S3 matrix.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric, and 1
!    if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) B(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) x(m)

  b(1:n) = 0.0D+00

  do k = 1, nz_num
    i = col(k)
    j = row(k)
    b(i) = b(i) + a(k) * x(j)
  end do
!
!  Handle the symmetric option.
!
  if ( isym == 1 .and. m == n ) then
    do k = 1, nz_num
      i = row(k)
      j = col(k)
      if ( i /= j ) then
        b(i) = b(i) + a(k) * x(j)
      end if
    end do
  end if

  return
end
subroutine r8s3_write ( n, nz_num, isym, row, col, a, output_file )

!*****************************************************************************80
!
!! R8S3_WRITE writes a square R8S3 matrix to a file.
!
!  Discussion:
!
!    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ISYM, is 0 if the matrix is not symmetric, and 1
!    if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
!    Input, character ( len = * ) OUTPUT_FILE, the name of the file to which
!    the information is to be written.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  character ( len = * ) output_file
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) row(nz_num)

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8S3_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' &
      // trim ( output_file ) // '".'
    stop
  end if

  do k = 1, nz_num
    write ( output_unit, '(2x,i8,2x,i8,2x,g16.8)' ) row(k), col(k), a(k)
  end do

  close ( unit = output_unit )

  return
end
subroutine r8sd_cg ( n, ndiag, offset, a, b, x )

!*****************************************************************************80
!
!! R8SD_CG uses the conjugate gradient method on an R8SD linear system.
!
!  Discussion:
!
!    The R8SD storage format is for symmetric matrices whose only nonzero entries
!    occur along a few diagonals, but for which these diagonals are not all
!    close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0, and 
!    each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!
!    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
!    we then create an array B that has N rows and NDIAG columns, and simply
!    "collapse" the matrix A to the left:
!
!    For the conjugate gradient method to be applicable, the matrix A must 
!    be a positive definite symmetric matrix.
!
!    The method is designed to reach the solution to the linear system
!      A * x = b
!    after N computational steps.  However, roundoff may introduce
!    unacceptably large errors for some problems.  In such a case,
!    calling the routine a second time, using the current solution estimate
!    as the new starting guess, should result in improved results.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    FS Beckman,
!    The Solution of Linear Equations by the Conjugate Gradient Method,
!    Mathematical Methods for Digital Computers, pages 62-72.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8SD matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side vector.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, an estimate for the solution, which may be 0.
!    On output, the approximate solution vector.  Note that repeated
!    calls to this routine, using the value of X output on the previous
!    call, MAY improve the solution.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) ap(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) beta
  integer ( kind = 4 ) it
  integer ( kind = 4 ) offset(ndiag)
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) pap
  real ( kind = 8 ) pr
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rap
  real ( kind = 8 ) x(n)
!
!  Initialize
!    AP = A * x,
!    R  = b - A * x,
!    P  = b - A * x.
!
  call r8sd_mxv ( n, ndiag, offset, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP = A*P.
!
    call r8sd_mxv ( n, ndiag, offset, a, p, ap )
!
!  Compute the dot products
!    PAP = P*AP,
!    PR  = P*R
!  Set
!    ALPHA = PR / PAP.
!
    pap = sum ( p(1:n) * ap(1:n) )
    pr = sum ( p(1:n) * r(1:n) )

    if ( pap == 0.0D+00 ) then
      return
    end if

    alpha = pr / pap
!
!  Set
!    X = X + ALPHA * P
!    R = R - ALPHA * AP.
!
    x(1:n) = x(1:n) + alpha * p(1:n)
    r(1:n) = r(1:n) - alpha * ap(1:n)
!
!  Compute the vector dot product
!    RAP = R*AP
!  Set
!    BETA = - RAP / PAP.
!
    rap = sum ( r(1:n) * ap(1:n) )

    beta = - rap / pap
!
!  Update the perturbation vector
!    P = R + BETA * P.
!
    p(1:n) = r(1:n) + beta * p(1:n)

  end do

  return
end
subroutine r8sd_indicator ( n, ndiag, offset, a )

!*****************************************************************************80
!
!! R8SD_INDICATOR sets up an R8SD indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8SD storage format is for symmetric matrices whose only nonzero 
!    entries occur along a few diagonals, but for which these diagonals are not 
!    all close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0, and 
!    each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!
!    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
!    we then create an array B that has N rows and NDIAG columns, and simply
!    "collapse" the matrix A to the left:
!
!  Example:
!
!    The "offset" value is printed above each column.
!
!    Original matrix               New Matrix
!
!       0   1   2   3   4   5       0   1   3   5
!
!      11  12   0  14   0  16      11  12  14  16
!      21  22  23   0  25   0      22  23  25  --
!       0  32  33  34   0  36      33  34  36  --
!      41   0  43  44  45   0      44  45  --  --
!       0  52   0  54  55  56      55  56  --  --
!      61   0  63   0  65  66      66  --  --  --
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2004
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Output, real ( kind = 8 ) A(N,NDIAG), the R8SD matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) offset(ndiag)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        a(i,diag) = real ( fac * i + j, kind = 8 )
      else
        a(i,diag) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8sd_mxv ( n, ndiag, offset, a, x, b )

!*****************************************************************************80
!
!! R8SD_MXV multiplies an R8SD matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SD storage format is for symmetric matrices whose only nonzero 
!    entries occur along a few diagonals, but for which these diagonals are not 
!    all close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0, and 
!    each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!
!    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
!    we then create an array B that has N rows and NDIAG columns, and simply
!    "collapse" the matrix A to the left:
!
!  Example:
!
!    The "offset" value is printed above each column.
!
!    Original matrix               New Matrix
!
!       0   1   2   3   4   5       0   1   3   5
!
!      11  12   0  14   0  16      11  12  14  16
!      21  22  23   0  25   0      22  23  25  --
!       0  32  33  34   0  36      33  34  36  --
!      41   0  43  44  45   0      44  45  --  --
!       0  52   0  54  55  56      55  56  --  --
!      61   0  63   0  65  66      66  --  --  --
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 October 1998
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8SD matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jdiag
  integer ( kind = 4 ) offset(ndiag)
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  do i = 1, n
    do jdiag = 1, ndiag
      j = i + offset(jdiag)
      if ( 1 <= j .and. j <= n ) then
        b(i) = b(i) + a(i,jdiag) * x(j)
        if ( offset(jdiag) /= 0 ) then
          b(j) = b(j) + a(i,jdiag) * x(i)
        end if
      end if
    end do
  end do

  return
end
subroutine r8sd_print ( n, ndiag, offset, a, title )

!*****************************************************************************80
!
!! R8SD_PRINT prints an R8SD matrix.
!
!  Discussion:
!
!    The R8SD storage format is for symmetric matrices whose only nonzero 
!    entries occur along a few diagonals, but for which these diagonals are not 
!    all close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0, and 
!    each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!
!    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
!    we then create an array B that has N rows and NDIAG columns, and simply
!    "collapse" the matrix A to the left:
!
!  Example:
!
!    The "offset" value is printed above each column.
!
!    Original matrix               New Matrix
!
!       0   1   2   3   4   5       0   1   3   5
!
!      11  12   0  14   0  16      11  12  14  16
!      21  22  23   0  25   0      22  23  25  --
!       0  32  33  34   0  36      33  34  36  --
!      41   0  43  44  45   0      44  45  --  --
!       0  52   0  54  55  56      55  56  --  --
!      61   0  63   0  65  66      66  --  --  --
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than N.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the 
!    diagonal storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8SD matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) offset(ndiag)
  character ( len = * ) title

  call r8sd_print_some ( n, ndiag, offset, a, 1, 1, n, n, title )

  return
end
subroutine r8sd_print_some ( n, ndiag, offset, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8SD_PRINT_SOME prints some of an R8SD matrix.
!
!  Discussion:
!
!    The R8SD storage format is for symmetric matrices whose only nonzero 
!    entries occur along a few diagonals, but for which these diagonals are not 
!    all close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0, and 
!    each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!
!    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
!    we then create an array B that has N rows and NDIAG columns, and simply
!    "collapse" the matrix A to the left:
!
!  Example:
!
!    The "offset" value is printed above each column.
!
!    Original matrix               New Matrix
!
!       0   1   2   3   4   5       0   1   3   5
!
!      11  12   0  14   0  16      11  12  14  16
!      21  22  23   0  25   0      22  23  25  --
!       0  32  33  34   0  36      33  34  36  --
!      41   0  43  44  45   0      44  45  --  --
!       0  52   0  54  55  56      55  56  --  --
!      61   0  63   0  65  66      66  --  --  --
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than N.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8SD matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jdiag
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) off
  integer ( kind = 4 ) offset(ndiag)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0D+00
        off = j - i
        do jdiag = 1, ndiag
          if ( off == offset(jdiag) ) then
            aij = a(i,jdiag)
          else if ( off == - offset(jdiag) ) then
            aij = a(j,jdiag)
          end if
        end do

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8sd_random ( n, ndiag, offset, seed, a )

!*****************************************************************************80
!
!! R8SD_RANDOM randomizes an R8SD matrix.
!
!  Discussion:
!
!    The R8SD storage format is for symmetric matrices whose only nonzero 
!    entries occur along a few diagonals, but for which these diagonals are not 
!    all close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0, and 
!    each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!
!    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
!    we then create an array B that has N rows and NDIAG columns, and simply
!    "collapse" the matrix A to the left:
!
!  Example:
!
!    The "offset" value is printed above each column.
!
!    Original matrix               New Matrix
!
!       0   1   2   3   4   5       0   1   3   5
!
!      11  12   0  14   0  16      11  12  14  16
!      21  22  23   0  25   0      22  23  25  --
!       0  32  33  34   0  36      33  34  36  --
!      41   0  43  44  45   0      44  45  --  --
!       0  52   0  54  55  56      55  56  --  --
!      61   0  63   0  65  66      66  --  --  --
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 1998
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) A(N,NDIAG), the R8SD matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) offset(ndiag)
  integer ( kind = 4 ) seed

  do i = 1, n
    do j = 1, ndiag
      jj = i + offset(j)
      if ( 1 <= jj .and. jj <= n ) then
        a(i,j) = r8_uniform_01 ( seed )
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8sd_to_r8ge ( n, ndiag, offset, a, b )

!*****************************************************************************80
!
!! R8SD_TO_R8GE copies an R8SD matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8SD storage format is for symmetric matrices whose only nonzero 
!    entries occur along a few diagonals, but for which these diagonals are not 
!    all close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0, and 
!    each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!
!    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
!    we then create an array B that has N rows and NDIAG columns, and simply
!    "collapse" the matrix A to the left:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Example:
!
!    The "offset" value is printed above each column.
!
!    Original matrix               New Matrix
!
!       0   1   2   3   4   5       0   1   3   5
!
!      11  12   0  14   0  16      11  12  14  16
!      21  22  23   0  25   0      22  23  25  --
!       0  32  33  34   0  36      33  34  36  --
!      41   0  43  44  45   0      44  45  --  --
!       0  52   0  54  55  56      55  56  --  --
!      61   0  63   0  65  66      66  --  --  --
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 1998
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
!    Input, integer ( kind = 4 ) NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer ( kind = 4 ) OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Input, real ( kind = 8 ) A(N,NDIAG), the R8SD matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) a(n,ndiag)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) offset(ndiag)

  b(1:n,1:n) = 0.0D+00

  do i = 1, n
    do j = 1, ndiag
      jj = i + offset(j)
      if ( 1 <= jj .and. jj <= n ) then
        b(i,jj) = a(i,j)
        if ( i /= jj ) then
          b(jj,i) = a(i,j)
        end if
      end if
    end do
  end do

  return
end
subroutine r8sm_ml ( n, a_lu, u, v, pivot, x, b, job )

!*****************************************************************************80
!
!! R8SM_ML multiplies a factored square R8SM matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
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
!    Input, real ( kind = 8 ) U(N), V(N), the Sherman Morrison vectors.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector computed by R8GE_FA.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(N), the result of the multiplication.
!
!    Input, integer ( kind = 4 ) JOB, specifies the operation to be done:
!    JOB = 0, compute (A-u*v') * x.
!    JOB nonzero, compute (A-u*v')' * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) job
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) x(n)

  call r8ge_ml ( n, a_lu, pivot, x, b, job )

  if ( job == 0 ) then

    b(1:n) = b(1:n) - u(1:n) * sum ( v(1:n) * x(1:n) )

  else

    b(1:n) = b(1:n) - v(1:n) * sum ( u(1:n) * x(1:n) )

  end if

  return
end
subroutine r8sm_mxv ( m, n, a, u, v, x, b )

!*****************************************************************************80
!
!! R8SM_MXV multiplies an R8SM matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the R8SM matrix.
!
!    Input, real ( kind = 8 ) U(M), V(N), the R8SM vectors U and V.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by (A-u*v').
!
!    Output, real ( kind = 8 ) B(M), the product (A-u*v') * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) x(n)

  b(1:m) = matmul ( a(1:m,1:n), x(1:n) ) &
    - u(1:m) * sum ( v(1:n) * x(1:n) )

  return
end
subroutine r8sm_print ( m, n, a, u, v, title )

!*****************************************************************************80
!
!! R8SM_PRINT prints an R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the R8SM matrix.
!
!    Input, real ( kind = 8 ) U(M), V(N), the R8SM vectors.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)

  call r8sm_print_some ( m, n, a, u, v, 1, 1, m, n, title )

  return
end
subroutine r8sm_print_some ( m, n, a, u, v, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8SM_PRINT_SOME prints some of an R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the R8SM matrix.
!
!    Input, real ( kind = 8 ) U(M), V(N), the R8SM vectors.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = a(i,j) - u(i) * v(j)

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8sm_random ( m, n, seed, a, u, v )

!*****************************************************************************80
!
!! R8SM_RANDOM randomizes an R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) A(M,N), the R8SM matrix.
!
!    Output, real ( kind = 8 ) U(M), V(N), the R8SM vectors.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)

  call r8mat_uniform_01 ( m, n, seed, a )
  
  call r8vec_uniform_01 ( m, seed, u )

  call r8vec_uniform_01 ( n, seed, v )

  return
end
subroutine r8sm_sl ( n, a_lu, u, v, b, ierror, pivot, job )

!*****************************************************************************80
!
!! R8SM_SL solves a square R8SM system that has been factored.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!    It is assumed that A has been decomposed into its LU factors
!    by R8GE_FA.  The Sherman Morrison formula allows
!    us to solve linear systems involving (A-u*v') by solving linear
!    systems involving A and adjusting the results.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash
!    Numerical Methods and Software,
!    Prentice Hall, 1989
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A_LU(N,N), the LU factors from R8GE_FA.
!
!    Input, real ( kind = 8 ) U(N), V(N), the R8SM vectors U and V.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error occurred.  The solution was successfully computed.
!    1, an error occurred.  1 - v' * Inverse(A) * u = 0.
!    The solution was not computed.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector produced by R8GE_FA.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve (A-u*v') * X = B.
!    nonzero, solve (A-u*v') * X = B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(n,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) beta
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) job_local
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) w(n)

  ierror = 0

  if ( job == 0 ) then
!
!  Solve A' * w = v.
!
    w(1:n) = v(1:n)

    job_local = 1
    call r8ge_sl ( n, a_lu, pivot, w, job_local )
!
!  Set beta = w' * b.
!
    beta = sum ( w(1:n) * b(1:n) )
!
!  Solve A * b = b.
!
    job_local = 0
    call r8ge_sl ( n, a_lu, pivot, b, job_local )
!
!  Solve A * w = u.
!
    w(1:n) = u(1:n)

    job_local = 0
    call r8ge_sl ( n, a_lu, pivot, w, job_local )
!
!  Set alpha = 1 / ( 1 - v' * w ).
!
    alpha = 1.0D+00 - sum ( v(1:n) * w(1:n) )

  else
!
!  Solve A * w = u.
!
    w(1:n) = u(1:n)

    job_local = 0
    call r8ge_sl ( n, a_lu, pivot, w, job_local )
!
!  Set beta = w' * b.
!
    beta = sum ( w(1:n) * b(1:n) )
!
!  Solve A' * b = b.
!
    job_local = 1
    call r8ge_sl ( n, a_lu, pivot, b, job_local )
!
!  Solve A' * w = v.
!
    w(1:n) = v(1:n)

    job_local = 1
    call r8ge_sl ( n, a_lu, pivot, w, job_local )
!
!  Set alpha = 1 / ( 1 - u' * w ).
!
    alpha = 1.0D+00 - sum ( u(1:n) * w(1:n) )

  end if

  if ( alpha == 0.0D+00 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SM_SL - Fatal error!'
    write ( *, '(a)' ) '  The divisor ALPHA is zero.'
    stop
  end if

  alpha = 1.0D+00 / alpha
!
!  Set b = b + alpha * beta * w.
!
  b(1:n) = b(1:n) + alpha * beta * w(1:n)

  return
end
subroutine r8sm_to_r8ge ( m, n, a, u, v, b )

!*****************************************************************************80
!
!! R8SM_TO_R8GE copies an R8SM matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the R8SM matrix.
!
!    Input, real ( kind = 8 ) U(M), V(N), the R8SM vectors.
!
!    Output, real ( kind = 8 ) B(M,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)

  do i = 1, m
    b(i,1:n) = a(i,1:n) - u(i) * v(1:n)
  end do

  return
end
subroutine r8sm_vxm ( m, n, a, u, v, x, b )

!*****************************************************************************80
!
!! R8SM_VXM multiplies an R8VEC by an R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the R8SM matrix.
!
!    Input, real ( kind = 8 ) U(M), V(N), the R8SM vectors.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(N), the product (A-u*v')' * X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) u(m)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) x(m)

  b(1:n) = matmul ( transpose ( a(1:m,1:n) ), x(1:m) ) &
    - v(1:n) * sum ( u(1:m) * x(1:m) )

  return
end
subroutine r8sp_check ( m, n, nz_num, row, col, check )

!*****************************************************************************80
!
!! R8SP_CHECK checks that an R8SP matrix data structure is properly sorted.
!
!  Discussion:
!
!    This routine assumes that the data structure has been sorted,
!    so that the entries of ROW are ascending sorted, and that the
!    entries of COL are ascending sorted, within the group of entries
!    that have a common value of ROW.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Output, logical CHECK, is TRUE if the matrix is properly defined.
!
  implicit none

  integer ( kind = 4 ) nz_num

  logical check
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)

  check = .true.
!
!  Check 1 <= ROW(*) <= M.
!
  do k = 1, nz_num

    if ( row(k) < 1 .or. m < row(k) ) then
      check = .false.
      return
    end if

  end do
!
!  Check 1 <= COL(*) <= N.
!
  do k = 1, nz_num

    if ( col(k) < 1 .or. n < col(k) ) then
      check = .false.
      return
    end if

  end do
!
!  Check that ROW(K) <= ROW(K+1).
!
  do k = 1, nz_num - 1

    if ( row(k+1) < row(k) ) then
      check = .false.
      return
    end if

  end do
!
!  Check that, if ROW(K) == ROW(K+1), that COL(K) < COL(K+1).
!
  do k = 1, nz_num - 1

    if ( row(k) == row(k+1) ) then
      if ( col(k+1) <= col(k) ) then
        check = .false.
        return
      end if
    end if

  end do

  return
end
subroutine r8sp_ij_to_k ( nz_num, row, col, i, j, k )

!*****************************************************************************80
!
!! R8SP_IJ_TO_K seeks the compressed index of the (I,J) entry of A.
!
!  Discussion:
!
!    If A(I,J) is nonzero, then its value is stored in location K.
!
!    This routine searches the R8SP storage structure for the index K
!    corresponding to (I,J), returning -1 if no such entry was found.
!
!    This routine assumes that the data structure has been sorted,
!    so that the entries of ROW are ascending sorted, and that the
!    entries of COL are ascending sorted, within the group of entries
!    that have a common value of ROW.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Input, integer ( kind = 4 ) I, J, the row and column indices of the
!    matrix entry.
!
!    Output, integer ( kind = 4 ) K, the R8SP index of the (I,J) entry.
!
  implicit none

  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) md
  integer ( kind = 4 ) row(nz_num)

  lo = 1
  hi = nz_num

  do

    if ( hi < lo ) then
      k = -1
      exit
    end if

    md = ( lo + hi ) / 2

    if ( row(md) < i .or. ( row(md) == i .and. col(md) < j ) ) then
      lo = md + 1
    else if ( i < row(md) .or. ( row(md) == i .and. j < col(md) ) ) then
      hi = md - 1
    else
      k = md
      exit
    end if

  end do

  return
end
subroutine r8sp_indicator ( m, n, nz_num, row, col, a )

!*****************************************************************************80
!
!! R8SP_INDICATOR sets up an R8SP indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do k = 1, nz_num

    i = row(k)
    j = col(k)
    a(k) = real ( fac * i + j, kind = 8 )

  end do

  return
end
subroutine r8sp_mxv ( m, n, nz_num, row, col, a, x, b )

!*****************************************************************************80
!
!! R8SP_MXV multiplies an R8SP matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(M), the product vector A*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) x(n)

  b(1:m) = 0.0D+00

  do k = 1, nz_num

    i = row(k)
    j = col(k)
    b(i) = b(i) + a(k) * x(j)

  end do

  return
end
subroutine r8sp_print ( m, n, nz_num, row, col, a, title )

!*****************************************************************************80
!
!! R8SP_PRINT prints an R8SP matrix.
!
!  Discussion:
!
!    This version of R8SP_PRINT has been specifically modified to allow,
!    and correctly handle, the case in which a single matrix location
!    A(I,J) is referenced more than once by the sparse matrix structure.
!    In such cases, the routine prints out the sum of all the values.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
! 
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) title

  call r8sp_print_some ( m, n, nz_num, row, col, a, 1, 1, m, n, title )

  return
end
subroutine r8sp_print_some ( m, n, nz_num, row, col, a, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! R8SP_PRINT_SOME prints some of an R8SP matrix.
!
!  Discussion:
!
!    This version of R8SP_PRINT_SOME has been specifically modified to allow,
!    and correctly handle, the case in which a single matrix location
!    A(I,J) is referenced more than once by the sparse matrix structure.
!    In such cases, the routine prints out the sum of all the values.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements
!    in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) aij(incx)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '
    write ( *, '(''  Col:  '',5(i7,7x))' ) ( j, j = j2lo, j2hi )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      aij(1:inc) = 0.0D+00
!
!  Is matrix entry K actually the value of A(I,J), with J2LO <= J <= J2HI?
!  Because MATLAB seems to allow for multiple (I,J,A) entries, we have
!  to sum up what we find.
! 
      do k = 1, nz_num

        if ( i == row(k) .and. &
             j2lo <= col(k) .and. &
             col(k) <= j2hi ) then 

          j2 = col(k) - j2lo + 1
          aij(j2) = aij(j2) + a(k)

        end if

      end do

      if ( any ( aij(1:inc) /= 0.0D+00 ) ) then
        write ( *, '(i5,1x,5g14.6)' ) i, aij(1:inc)
      end if

    end do

  end do

  return
end
subroutine r8sp_random ( m, n, nz_num, row, col, seed, a )

!*****************************************************************************80
!
!! R8SP_RANDOM sets a random R8SP matrix.
!
!  Discussion:
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements 
!    in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)
  integer ( kind = 4 ) seed

  call r8vec_uniform_01 ( nz_num, seed, a )

  return
end
subroutine r8sp_read ( input_file, m, n, nz_num, row, col, a )

!*****************************************************************************80
!
!! R8SP_READ reads an R8SP matrix from a file.
!
!  Discussion:
!
!    This routine needs the value of NZ_NUM, which can be determined
!    by a call to R8SP_READ_SIZE.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file to be read.
!
!    Unused, integer M, N, the number of rows and columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements 
!    in the matrix.
!
!    Output, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  character ( len = * ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SP_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( input_file ) // '".'
    stop
  end if

  do k = 1, nz_num

    read ( input_unit, *, iostat = ios ) row(k), col(k), a(k)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8SP_READ - Fatal error!'
      write ( *, '(a,i8)' ) '  I/O error while reading record ', k
      stop
    end if

  end do

  close ( unit = input_unit )

  return
end
subroutine r8sp_read_size ( input_file, m, n, nz_num )

!*****************************************************************************80
!
!! R8SP_READ_SIZE reads the size of an R8SP matrix from a file.
!
!  Discussion:
!
!    The value of NZ_NUM is simply the number of records in the input file.
!
!    The values of M and N are determined as the maximum entry in the row 
!    and column vectors.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file to 
!    be read.
!
!    Output, integer ( kind = 4 ) M, N, the number of rows and columns
!    of the matrix.
!
!    Output, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements 
!    in the matrix.
!
  implicit none

  real ( kind = 8 ) a_k
  integer ( kind = 4 ) col_k
  character ( len = * ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ) row_k

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SP_READ_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( input_file ) // '".'
    stop
  end if

  m = 0
  n = 0
  nz_num = 0

  do

    read ( input_unit, *, iostat = ios ) row_k, col_k, a_k

    if ( ios /= 0 ) then
      exit
    end if

    nz_num = nz_num + 1
    m = max ( m, row_k )
    n = max ( n, col_k )

  end do

  close ( unit = input_unit )

  return
end
subroutine r8sp_to_r8ge ( m, n, nz_num, row, col, a, b )

!*****************************************************************************80
!
!! R8SP_TO_R8GE converts an R8SP matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements
!    in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Output, real ( kind = 8 ) B(M,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(m,n)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)

  b(1:m,1:n) = 0.0D+00

  do k = 1, nz_num
    b(row(k),col(k)) = a(k)
  end do

  return
end
subroutine r8sp_to_r8ncf ( m, n, nz_num, row, col, a, rowcol )

!*****************************************************************************80
!
!! R8SP_TO_R8NCF converts an R8SP matrix to an R8NCF matrix.
!
!  Discussion:
!
!    The R8SP and R8NCF formats are essentially identical, except that
!    R8SP keeps separate ROW and COLUMN vectors, while R8NCF uses a single
!    ROWCOL array.  Therefore, the input values NZ_NUM and A used in
!    the R8SP representation can be regarded as part of the output
!    values used for the R8NCF representation.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Unused, integer M, N, the number of rows and columns of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements 
!    in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Output, integer ( kind = 4 ) ROWCOL(2,NZ_NUM), the R8NCF row and column
!    index vector.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)
  integer ( kind = 4 ) rowcol(2,nz_num)

  rowcol(1,1:nz_num) = row(1:nz_num)
  rowcol(2,1:nz_num) = col(1:nz_num)

  return
end
subroutine r8sp_vxm ( m, n, nz_num, row, col, a, x, b )

!*****************************************************************************80
!
!! R8SP_VXM multiplies an R8VEC times an R8SP matrix.
!
!  Discussion:
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product vector A'*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real ( kind = 8 ) x(m)

  b(1:n) = 0.0D+00

  do k = 1, nz_num

    i = row(k)
    j = col(k)
    b(j) = b(j) + a(k) * x(i)

  end do

  return
end
subroutine r8sp_write ( m, n, nz_num, row, col, a, output_file )

!*****************************************************************************80
!
!! R8SP_WRITE writes an R8SP matrix to a file.
!
!  Discussion:
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
!    Input, character ( len = * ) OUTPUT_FILE, the name of the file to which
!    the information is to be written.
!
  implicit none

  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = * ) output_file
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) row(nz_num)

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SP_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' &
      // trim ( output_file ) // '".'
    stop
  end if

  do k = 1, nz_num
    write ( output_unit, '(2x,i8,2x,i8,2x,g16.8)' ) row(k), col(k), a(k)
  end do

  close ( unit = output_unit )

  return
end
subroutine r8sp_zero ( m, n, nz_num, row, col, a )

!*****************************************************************************80
!
!! R8SP_ZERO zeros out an R8SP matrix.
!
!  Discussion:
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    Perhaps we must point out that when we say "nonzero" elements in the
!    discussion below, we mean "potentially nonzero" elements.  In other words,
!    the sparse matrix sets aside space for elements that are allowed to be
!    nonzero, but may of course take on zero values as well.  It's actually
!    really the entries for which we don't set aside space that we are sure
!    about.  Those entries are zero.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) row(nz_num)

  a(1:nz_num) = 0.0D+00

  return
end
subroutine r8sr_indicator ( n, nz, row, col, diag, off )

!*****************************************************************************80
!
!! R8SR_INDICATOR sets up an R8SR indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero elements
!    in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Output, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Output, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, n

    j = i
    diag(i) = real ( fac * i + j, kind = 8 )

    do k = row(i), row(i+1)-1
      j = col(k)
      off(k) = real ( fac * i + j, kind = 8 )
    end do

  end do

  return
end
subroutine r8sr_mxv ( n, nz, row, col, diag, off, x, b )

!*****************************************************************************80
!
!! R8SR_MXV multiplies an R8SR matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero 
!    elements in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix.
!
!    Input, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by the matrix.
!
!    Output, real ( kind = 8 ) B(N), the product A * X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)
  real ( kind = 8 ) x(n)

  b(1:n) = diag(1:n) * x(1:n)

  do i = 1, n
    do k = row(i), row(i+1)-1
      j = col(k)
      b(i) = b(i) + off(k) * x(j)
    end do
  end do

  return
end
subroutine r8sr_print ( n, nz, row, col, diag, off, title )

!*****************************************************************************80
!
!! R8SR_PRINT prints an R8SR matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero elements
!    in A.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of
!    the element in the corresponding position in A.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Input, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
!    Input, character ( len = * ) TITLE, a title.
! 
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)
  character ( len = * ) title

  call r8sr_print_some ( n, nz, row, col, diag, off, 1, 1, n, n, title )

  return
end
subroutine r8sr_print_some ( n, nz, row, col, diag, off, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! R8SR_PRINT_SOME prints some of an R8SR matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero elements
!    in A.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Input, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  real ( kind = 8 ) aij
  integer ( kind = 4 ) col(nz)
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
!  1) Assume everything is zero.
!
      aij = 0.0D+00
      do j2 = 1, inc
        write ( ctemp(j2), '(f8.0,6x)' ) aij
      end do
!
!  2) Insert the diagonal entry, if appropriate.
!
      if ( j2lo <= i .and. i <= j2hi ) then
        j2 = i - j2lo + 1
        aij = diag(i)
        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if
      end if
!
!  3) Now examine all the offdiagonal entries.
!
      do k = row(i), row(i+1)-1
        if ( j2lo <= col(k) .and. col(k) <= j2hi ) then 
          j2 = col(k) - j2lo + 1
          aij = off(k)
          if ( r8_is_int ( aij ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) aij
          else
            write ( ctemp(j2), '(g14.6)' ) aij
          end if
        end if
      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8sr_random ( n, nz, row, col, diag, off, seed )

!*****************************************************************************80
!
!! R8SR_RANDOM randomizes an R8SR matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero elements 
!    in A.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Output, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Output, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)
  integer ( kind = 4 ) seed

  do i = 1, n
    diag(i) = r8_uniform_01 ( seed )
    do j = row(i), row(i+1)-1
      off(j) = r8_uniform_01 ( seed )
    end do
  end do

  return
end
subroutine r8sr_to_r8ge ( n, nz, row, col, diag, off, b )

!*****************************************************************************80
!
!! R8SR_TO_R8GE converts an R8SR matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero 
!    elements in A.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Input, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SR_TO_R8GE - Fatal error!'
    write ( *, '(a,i8)' ) '  N is less than or equal to zero, N = ', n
    stop
  end if

  b(1:n,1:n) = 0.0D+00

  do i = 1, n
    b(i,i) = diag(i)
  end do

  do i = 1, n
    do j = row(i), row(i+1)-1
      b(i,col(j)) = off(j)
    end do
  end do

  return
end
subroutine r8sr_vxm ( n, nz, row, col, diag, off, x, b )

!*****************************************************************************80
!
!! R8SR_VXM multiplies an R8VEC times an R8SR matrix.
!
!  Discussion:
!
!    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
!    The off-diagonal entries of row I are stored in entries ROW(I)
!    through ROW(I+1)-1 of OFF.  COL(J) records the column index
!    of the entry in A(J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) NZ, the number of offdiagonal nonzero 
!    elements in A.
!
!    Input, integer ( kind = 4 ) ROW(N+1).  The nonzero offdiagonal elements 
!    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
!
!    Input, integer ( kind = 4 ) COL(NZ), contains the column index of the 
!    element in the corresponding position in A.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of A.
!
!    Input, real ( kind = 8 ) OFF(NZ), the off-diagonal elements of A.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplies by A.
!
!    Output, real ( kind = 8 ) B(N), the product A' * X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz

  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) col(nz)
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) off(nz)
  integer ( kind = 4 ) row(n+1)
  real ( kind = 8 ) x(n)

  b(1:n) = diag(1:n) * x(1:n)

  do i = 1, n
    do k = row(i), row(i+1)-1
      j = col(k)
      b(j) = b(j) + off(k) * x(i)
    end do
  end do

  return
end
subroutine r8ss_error ( n, na, diag, ierror )

!*****************************************************************************80
!
!! R8SS_ERROR checks dimensions for an R8SS matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 1998
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
!    Input, integer ( kind = 4 ) NA, the dimension of the array A.
!    NA must be at least N.
!
!    Input, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal 
!    elements.
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.
!    0, no error.
!    1, N is less than 1.
!    2, NA is less than N.
!    3, DIAG(1) is not 1.
!    4, the elements of DIAG are not strictly increasing.
!    5, DIAG(N) is greater than NA.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) na

  ierror = 0

  if ( n < 1 ) then
    ierror = 1
    return
  end if

  if ( na < n ) then
    ierror = 2
    return
  end if

  if ( diag(1) /= 1 ) then
    ierror = 3
    return
  end if

  do i = 1, n - 1
    if ( diag(i+1) <= diag(i) ) then
      ierror = 4
      return
    end if
  end do

  if ( na < diag(n) ) then
    ierror = 5
    return
  end if

  return
end
subroutine r8ss_indicator ( n, na, diag, a )

!*****************************************************************************80
!
!! R8SS_INDICATOR sets up an R8SS indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.   Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
!    although the actual storage needed will generally be about half of
!    that.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 2004
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
!    Output, integer ( kind = 4 ) NA, the dimension of the array A, which for
!    this special case will be the maximum, ( N * ( N + 1 ) ) / 2
!
!    Output, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Output, real ( kind = 8 ) A((N*(N+1))/2), the R8SS matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) na

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  na = 0

  do j = 1, n

    do i = 1, j
      na = na + 1
      a(na) = real ( fac * i + j, kind = 8 )
    end do

    diag(j) = na

  end do

  return
end
subroutine r8ss_mxv ( n, na, diag, a, x, b )

!*****************************************************************************80
!
!! R8SS_MXV multiplies an R8SS matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 1998
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
!    Input, integer ( kind = 4 ) NA, the dimension of the array A.
!    NA must be at least N.
!
!    Input, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Input, real ( kind = 8 ) A(NA), the R8SS matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A*x.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) na

  real ( kind = 8 ) a(na)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) diagold
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  b(1:n) = 0.0D+00

  diagold = 0
  k = 0

  do j = 1, n

    ilo = j + 1 + diagold - diag(j)

    do i = ilo, j - 1
      k = k + 1
      b(i) = b(i) + a(k) * x(j)
      b(j) = b(j) + a(k) * x(i)
    end do

    k = k + 1
    b(j) = b(j) + a(k) * x(j)

    diagold = diag(j)

  end do

  return
end
subroutine r8ss_print ( n, na, diag, a, title )

!*****************************************************************************80
!
!! R8SS_PRINT prints an R8SS matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
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
!    Input, integer ( kind = 4 ) NA, the dimension of the array A.
!
!    Input, integer ( kind = 4 ) DIAG(N), the indices in A of the N 
!    diagonal elements.
!
!    Input, real ( kind = 8 ) A(NA), the R8SS matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) na
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(na)
  integer ( kind = 4 ) diag(n)
  character ( len = * ) title

  call r8ss_print_some ( n, na, diag, a, 1, 1, n, n, title )

  return
end
subroutine r8ss_print_some ( n, na, diag, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8SS_PRINT_SOME prints some of an R8SS matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
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
!    Input, integer ( kind = 4 ) NA, the dimension of the array A.
!
!    Input, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Input, real ( kind = 8 ) A(NA), the R8SS matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) na
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(na)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ijm1
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0D+00

        if ( j < i ) then
          if ( i == 1 ) then
            ijm1 = 0
          else
            ijm1 = diag(i-1)
          end if
          ij = diag(i)
          if ( ijm1 < ij+j-i ) then
            aij = a(ij+j-i)
          end if
        else if ( j == i ) then
          ij = diag(j)
          aij = a(ij)
        else if ( i < j ) then
          if ( j == 1 ) then
            ijm1 = 0
          else
            ijm1 = diag(j-1)
          end if
          ij = diag(j)
          if ( ijm1 < ij+i-j ) then
            aij = a(ij+i-j)
          end if
        end if

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8ss_random ( n, na, diag, a, seed )

!*****************************************************************************80
!
!! R8SS_RANDOM randomizes an R8SS matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
!    although the actual storage needed will generally be about half of
!    that.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
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
!    Output, integer ( kind = 4 ) NA, the dimension of the array A.
!    NA will be at least N and no greater than ( N * ( N + 1 ) ) / 2.
!
!    Output, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Output, real ( kind = 8 ) A((N*(N+1))/2), the R8SS matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) na

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) diagold
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
!
!  Set the values of DIAG.
!
  diag(1) = 1
  na = 1
  do i = 2, n
    k = i4_uniform ( 1, i, seed )
    diag(i) = diag(i-1) + k
    na = na + k
  end do
!
!  Now set the values of A.
!
  diagold = 0
  k = 0

  do j = 1, n

    ilo = j + 1 + diagold - diag(j)

    do i = ilo, j
      k = k + 1
      a(k) = r8_uniform_01 ( seed )
    end do

    diagold = diag(j)

  end do

  return
end
subroutine r8ss_to_r8ge ( n, na, diag, a, b  )

!*****************************************************************************80
!
!! R8SS_TO_R8GE copies an R8SS matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8SS storage format is used for real symmetric skyline matrices.
!    This storage is appropriate when the nonzero entries of the
!    matrix are generally close to the diagonal, but the number
!    of nonzeroes above each diagonal varies in an irregular fashion.
!
!    In this case, the strategy is essentially to assign column J
!    its own bandwidth, and store the strips of nonzeros one after
!    another.  Note that what's important is the location of the
!    furthest nonzero from the diagonal.  A slot will be set up for
!    every entry between that and the diagonal, whether or not
!    those entries are zero.
!
!    A skyline matrix can be Gauss-eliminated without disrupting
!    the storage scheme, as long as no pivoting is required.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Example:
!
!    11   0  13  0 15
!     0  22  23  0  0
!    31  32  33 34  0
!     0   0  43 44  0
!    51   0   0  0 55
!
!    A = ( 11 | 22 | 13, 23, 33 | 34, 44 | 15, 0, 0, 0, 55 )
!    NA = 12
!    DIAG = ( 1, 2, 5, 7, 12 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 1998
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
!    Input, integer ( kind = 4 ) NA, the dimension of the array A.
!    NA must be at least N.
!
!    Input, integer ( kind = 4 ) DIAG(N), the indices in A of the N diagonal
!    elements.
!
!    Input, real ( kind = 8 ) A(NA), the R8SS matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) na

  real ( kind = 8 ) a(na)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) diagold
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  diagold = 0
  k = 0

  do j = 1, n

    ilo = j + 1 + diagold - diag(j)

    b(1:ilo-1,j) = 0.0D+00
    b(j,1:ilo-1) = 0.0D+00

    do i = ilo, j - 1
      k = k + 1
      b(i,j) = a(k)
      b(j,i) = a(k)
    end do

    k = k + 1
    b(j,j) = a(k)

    diagold = diag(j)

  end do

  return
end
subroutine r8sto_indicator ( n, a )

!*****************************************************************************80
!
!! R8STO_INDICATOR sets up an R8STO indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2004
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
!    Output, real ( kind = 8 ) A(N), the R8STO matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  i = 1
  k = 0
  do j = 1, n
    k = k + 1
    a(k) = real ( fac * i + j, kind = 8 )
  end do
  
  return
end
subroutine r8sto_inverse ( n, a, b )

!*****************************************************************************80
!
!! R8STO_INVERSE computes the inverse of an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!    For this routine, the matrix is also required to be positive definite.
!
!    The original implementation of the algorithm assumed that the
!    diagonal element was 1.  The algorithm has been modified so that
!    this is no longer necessary.
!
!    The inverse matrix is NOT guaranteed to be a Toeplitz matrix.  
!    It is guaranteed to be symmetric and persymmetric.
!    The inverse matrix is returned in general storage, that is,
!    as an "SGE" matrix.
!
!  Example:
!
!    To compute the inverse of
!
!     1.0 0.5 0.2
!     0.5 1.0 0.5
!     0.2 0.5 1.0
!
!    we input:
!
!      N = 3
!      A = (/ 1.0, 0.5, 0.2 /)
!
!    with output:
!
!      B(1:3,1:3) = (/ 75, -40,   5,
!                     -40,  96, -40,
!                       5, -40,  75 /) / 56
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Section 4.7.3, "Computing the Inverse",
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, real ( kind = 8 ) A(N), the R8STO matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the inverse of the matrix, in R8GE format.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a2(n-1)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) v(n)

  a2(1:n-1) = a(2:n) / a(1)

  call r8sto_yw_sl ( n - 1, a2, v )
!
!  Compute the N-th entry of V.
!
  v(n) = 1.0D+00 / ( 1.0D+00 + sum ( a2(1:n-1) * v(1:n-1) ) )
!
!  Reverse and scale entries 1 through N-1.
!
  v(1:n-1) = v(n-1:1:-1)

  v(1:n-1) = v(n) * v(1:n-1)
!
!  Set the boundaries of B.
!
  b(1,1:n) = v(n:1:-1)
  b(n,1:n) = v(1:n)
  b(2:n-1,1) = v(n-1:2:-1)
  b(2:n-1,n) = v(2:n-1)
!
!  Fill the interior.
!
  do i = 2, 1 + ( ( n - 1 ) / 2 )
    do j = i, n - i + 1
      b(i,j) = b(i-1,j-1) + ( v(n+1-j) * v(n+1-i) - v(i-1) * v(j-1) ) / v(n)
      b(j,i) = b(i,j)
      b(n+1-i,n+1-j) = b(i,j)
      b(n+1-j,n+1-i) = b(i,j)
    end do
  end do
!
!  Scale B.
!
  b(1:n,1:n) = b(1:n,1:n) / a(1)

  return
end
subroutine r8sto_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! R8STO_MXV multiplies an R8STO matrix by an R8VEC.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8STO matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n
    b(i) = sum ( a(i:2:-1) * x(1:i-1) ) + sum ( a(1:n+1-i) * x(i:n) )
  end do

  return
end
subroutine r8sto_print ( n, a, title )

!*****************************************************************************80
!
!! R8STO_PRINT prints an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
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
!    Input, real ( kind = 8 ) A(N), the R8STO matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  character ( len = * ) title

  call r8sto_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r8sto_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8STO_PRINT_SOME prints some of an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
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
!    Input, real ( kind = 8 ) A(N), the R8STO matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(1+j-i)
        else
          aij = a(1+i-j)
        end if

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8sto_random ( n, seed, a )

!*****************************************************************************80
!
!! R8STO_RANDOM randomizes an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2003
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
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) A(N), the R8STO matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) seed

  call r8vec_uniform_01 ( n, seed, a )

  return
end
subroutine r8sto_sl ( n, a, b, x )

!*****************************************************************************80
!
!! R8STO_SL solves an R8STO system.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!    For this routine, the matrix is also required to be positive definite.
!
!    This implementation of the algorithm assumes that the diagonal element
!    is 1.
!
!    The real symmetric Toeplitz matrix can be described by N numbers, which,
!    for convenience, we will label A(0:N-1).
!
!    Note that there is a typographical error in the presentation
!    of this algorithm in the reference, and another in the presentation
!    of a sample problem.  Both involve sign errors.  A minor error
!    makes the algorithm incorrect for the case N = 1.
!
!  Example:
!
!    To solve
!
!     1.0 0.5 0.2    x1    4.0
!     0.5 1.0 0.5 *  x2 = -1.0
!     0.2 0.5 1.0    x3    3.0
!
!    we input:
!
!      N = 3
!      A(0:N-1) = (/ 1.0, 0.5, 0.2 /)
!      B(1:3) = (/ 4.0, -1.0, 3.0 /)
!
!    with output:
!
!      X(1:3) = (/ 355, -376, 285 /) / 56
!             = (/ 6.339, -6.714, 5.089 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Section 4.7.3, "The General Right Hand Side Problem",
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, real ( kind = 8 ) A(N), the R8STO matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) beta
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  k = 0

  beta = 1.0D+00
  x(k+1) = b(k+1) / beta

  if ( k < n-1 ) then
    y(k+1) = -a(k+2) / beta
  end if

  do k = 1, n-1

    beta = ( 1.0D+00 - y(k) * y(k) ) * beta

    x(k+1) = ( b(k+1) - sum ( a(2:k+1) * x(k:1:-1) ) ) / beta

    x(1:k) = x(1:k) + x(k+1) * y(k:1:-1)

    if ( k < n - 1 ) then
      y(k+1) = ( -a(k+2) - sum ( a(2:k+1) * y(k:1:-1) ) ) / beta
      y(1:k) = y(1:k) + y(k+1) * y(k:1:-1)
    end if

  end do

  return
end
subroutine r8sto_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R8STO_TO_R8GE copies an R8STO matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8STO matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i

  do i = 1, n
    b(i,1:i-1) = a(i:2:-1)
    b(i,i:n) = a(1:n-i+1)
  end do

  return
end
subroutine r8sto_yw_sl ( n, b, x )

!*****************************************************************************80
!
!! R8STO_YW_SL solves the Yule-Walker equations for an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!    The matrix is also required to be positive definite.
!
!    This implementation of the algorithm assumes that the diagonal element
!    is 1.
!
!    The real symmetric Toeplitz matrix can be described by N numbers, which,
!    for convenience, we will label B(0:N-1).  We assume there is one more
!    number, B(N).  If we let A be the symmetric Toeplitz matrix whose first
!    row is B(0:N-1), then the Yule-Walker equations are:
!
!      A * X = -B(1:N)
!
!  Example:
!
!    To solve
!
!     1.0 0.5 0.2    x1   0.5
!     0.5 1.0 0.5 *  x2 = 0.2
!     0.2 0.5 1.0    x3   0.1
!
!    we input:
!
!      N = 3
!      B(1:3) = (/ 0.5, 0.2, 0.1 /)
!
!    with output:
!
!      X(1:3) = (/ -75, 12, -5 /) / 140
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Section 4.7.2, "Solving the Yule-Walker Equations",
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, real ( kind = 8 ) B(N), defines the linear system.  The first
!    entry of A is a 1, followed by B(1) through B(N-1).  The right hand
!    side of the system is -B(1:N).
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  x(1) = -b(1)
  beta = 1.0D+00
  alpha = -b(1)

  do i = 1, n-1
    beta = ( 1.0D+00 - alpha * alpha ) * beta
    alpha = - ( b(i+1) + sum ( b(i:1:-1) * x(1:i) ) ) / beta
    x(1:i) = x(1:i) + alpha * x(i:1:-1)
    x(i+1) = alpha
  end do

  return
end
subroutine r8to_indicator ( n, a )

!*****************************************************************************80
!
!! R8TO_INDICATOR sets up an R8TO indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8TO storage format is used for a real Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2004
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
!    Output, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  i = 1
  k = 0
  do j = 1, n
    k = k + 1
    a(k) = real ( fac * i + j, kind = 8 )
  end do

  j = 1
  do i = 2, n
    k = k + 1
    a(k) = real ( fac * i + j, kind = 8 )
  end do
  
  return
end
subroutine r8to_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! R8TO_MXV multiplies an R8TO matrix by an R8VEC.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  b(1) = sum ( a(1:n) * x(1:n) )

  do i = 2, n
    b(i) = sum ( a(n+i-1:n+1:-1) * x(1:i-1) ) &
         + sum ( a(1:n+1-i) * x(i:n) )
  end do

  return
end
subroutine r8to_print ( n, a, title )

!*****************************************************************************80
!
!! R8TO_PRINT prints an R8TO matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
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
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  character ( len = * ) title

  call r8to_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r8to_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8TO_PRINT_SOME prints some of an R8TO matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
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
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(j+1-i)
        else
          aij = a(n+i-j)
        end if

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8to_random ( n, seed, a )

!*****************************************************************************80
!
!! R8TO_RANDOM randomizes an R8TO matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2003
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
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) seed

  n2 = 2 * n - 1

  call r8vec_uniform_01 ( n2, seed, a )

  return
end
subroutine r8to_sl ( n, a, b, x, job )

!*****************************************************************************80
!
!! R8TO_SL solves an R8TO system.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 March 2001
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Input, real ( kind = 8 ) B(N) the right hand side vector.
!
!    Output, real ( kind = 8 ) X(N), the solution vector.  X and B may share the
!    same storage.
!
!    Input, integer ( kind = 4 ) JOB,
!    0 to solve A*X=B,
!    nonzero to solve A'*X=B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c1(n-1)
  real ( kind = 8 ) c2(n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) r5
  real ( kind = 8 ) r6
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    return
  end if
!
!  Solve the system with the principal minor of order 1.
!
  r1 = a(1)
  x(1) = b(1) / r1

  if ( n == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system with the Toeplitz matrix.
!
  do nsub = 2, n
!
!  Compute multiples of the first and last columns of the inverse of
!  the principal minor of order NSUB.
!
    if ( job == 0 ) then
      r5 = a(n+nsub-1)
      r6 = a(nsub)
    else
      r5 = a(nsub)
      r6 = a(n+nsub-1)
    end if

    if ( 2 < nsub ) then

      c1(nsub-1) = r2

      do i = 1, nsub - 2
        if ( job == 0 ) then
          r5 = r5 + a(n+i) * c1(nsub-i)
          r6 = r6 + a(i+1) * c2(i)
        else
          r5 = r5 + a(i+1) * c1(nsub-i)
          r6 = r6 + a(n+i) * c2(i)
        end if
      end do

    end if

    r2 = -r5 / r1
    r3 = -r6 / r1
    r1 = r1 + r5 * r3

    if ( 2 < nsub ) then

      r6 = c2(1)
      c2(nsub-1) = 0.0D+00

      do i = 2, nsub - 1
        r5 = c2(i)
        c2(i) = c1(i) * r3 + r6
        c1(i) = c1(i) + r6 * r2
        r6 = r5
      end do

    end if

    c2(1) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
    if ( job == 0 ) then
      r5 = sum ( a(n+1:n+nsub-1) * x(nsub-1:1:-1) )
    else
      r5 = sum ( a(2:nsub) * x(nsub-1:1:-1) )
    end if

    r6 = ( b(nsub) - r5 ) / r1

    x(1:nsub-1) = x(1:nsub-1) + c2(1:nsub-1) * r6
    x(nsub) = r6

  end do

  return
end
subroutine r8to_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R8TO_TO_R8GE copies an R8TO matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i

  do i = 1, n
    b(i,1:i-1) = a(n+i-1:n+1:-1)
    b(i,i:n) = a(1:n-i+1)
  end do

  return
end
subroutine r8to_vxm ( n, a, x, b )

!*****************************************************************************80
!
!! R8TO_VXM multiplies an R8VEC by an R8TO matrix.
!
!  Discussion:
!
!    The R8TO storage format is used for a Toeplitz matrix, which is constant
!    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
!    2*N-1 distinct entries.  The format stores the N elements of the first
!    row, followed by the N-1 elements of the first column (skipping the
!    entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(2*N-1), the R8TO matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A' * X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*n-1)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, n

    b(i) = sum ( a(i:1:-1) * x(1:i) ) + &
           sum ( a(n+1:2*n-i) * x(i+1:n) )

  end do

  return
end
subroutine r8ut_det ( n, a, det )

!*****************************************************************************80
!
!! R8UT_DET computes the determinant of an R8UT matrix.
!
!  Discussion:
!
!    The R8UT storage format is used for an M by N upper triangular 
!    matrix.  The format stores all M*N entries, even those which are zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 August 1999
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
!    Input, real ( kind = 8 ) A(N,N), the R8UT matrix.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det
  real ( kind = 8 ) diag(n)

  call r8mat_diag_get_vector ( n, a, diag )

  det = product ( diag(1:n) )

  return
end
subroutine r8ut_indicator ( m, n, a )

!*****************************************************************************80
!
!! R8UT_INDICATOR sets up an R8UT indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8UT storage format is used for an M by N upper triangular 
!    matrix.  The format stores all M*N entries, even those which are zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.  M and N must be positive.
!
!    Output, real ( kind = 8 ) A(M,N), the R8UT matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) fac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, m
    do j = 1, min ( i-1, n )
      a(i,j) = 0.0D+00
    end do
    do j = i, n
      a(i,j) = real ( fac * i + j, kind = 8 )
    end do
  end do

  return
end
subroutine r8ut_inverse ( n, a )

!*****************************************************************************80
!
!! R8UT_INVERSE computes the inverse of an R8UT matrix.
!
!  Discussion:
!
!    The R8UT storage format is used for an M by N upper triangular 
!    matrix.  The format stores all M*N entries, even those which are zero.
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
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, the upper triangular matrix to be inverted.
!    On output, the inverse of the upper triangular matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!
!  Check.
!
  do i = 1, n
    if ( a(i,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UT_INVERSE - Fatal error!'
      write ( *, '(a)' ) '  Zero diagonal element.'
      stop
    end if
  end do

  do j = n, 1, -1

    do i = n, 1, -1

      if ( j < i ) then

        a(i,j) = 0.0D+00

      else if ( i == j ) then

        a(i,j) = 1.0D+00 / a(i,j)

      else if ( i < j ) then

        a(i,j) = - sum ( a(i,i+1:j) * a(i+1:j,j) ) / a(i,i)

      end if

    end do
  end do

  return
end
subroutine r8ut_mxm ( n, a, b, c )

!*****************************************************************************80
!
!! R8UT_MXM multiplies two R8UT matrices.
!
!  Discussion:
!
!    The R8UT storage format is used for an M by N upper triangular 
!    matrix.  The format stores all M*N entries, even those which are zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N,N), B(N,N), the R8UT factor matrices.
!
!    Output, real ( kind = 8 ) C(N,N), the R8UT product matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, i - 1
      c(i,j) = 0.0D+00
    end do
    do j = i, n
      c(i,j) = sum ( a(i,i:j) * b(i:j,j) )
    end do
  end do

  return
end
subroutine r8ut_mxv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8UT_MXV multiplies an R8UT matrix by an R8VEC.
!
!  Discussion:
!
!    The R8UT storage format is used for an M by N upper triangular 
!    matrix.  The format stores all M*N entries, even those which are zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8UT matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(M), the product A * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  do i = 1, m
    b(i) = sum ( a(i,i:n) * x(i:n) )
  end do

  return
end
subroutine r8ut_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8UT_PRINT prints an R8UT matrix.
!
!  Discussion:
!
!    The R8UT storage format is used for an M by N upper triangular 
!    matrix.  The format stores all M*N entries, even those which are zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8UT matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8ut_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8ut_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8UT_PRINT_SOME prints some of an R8UT matrix.
!
!  Discussion:
!
!    The R8UT storage format is used for an M by N upper triangular 
!    matrix.  The format stores all M*N entries, even those which are zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8UT matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j < i ) then
          ctemp(j2) = '              '
        else if ( r8_is_int ( a(i,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8ut_random ( m, n, seed, a )

!*****************************************************************************80
!
!! R8UT_RANDOM randomizes an R8UT matrix.
!
!  Discussion:
!
!    The R8UT storage format is used for an M by N upper triangular 
!    matrix.  The format stores all M*N entries, even those which are zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.  M and N must be positive.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(M,N), the R8UT matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  do i = 1, m
    do j = 1, min ( i-1, n )
      a(i,j) = 0.0D+00
    end do
    do j = i, n
      a(i,j) = r8_uniform_01 ( seed )
    end do
  end do

  return
end
subroutine r8ut_sl ( n, a, b, job )

!*****************************************************************************80
!
!! R8UT_SL solves an R8UT system.
!
!  Discussion:
!
!    The R8UT storage format is used for an M by N upper triangular 
!    matrix.  The format stores all M*N entries, even those which are zero.
!
!    No factorization of the upper triangular matrix is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the R8UT matrix.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB, is 0 to solve the untransposed system,
!    nonzero to solve the transposed system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job

  if ( job == 0 ) then

    do j = n, 1, -1
      b(j) = b(j) / a(j,j)
      b(1:j-1) = b(1:j-1) - a(1:j-1,j) * b(j)
    end do

  else

    do j = 1, n
      b(j) = b(j) / a(j,j)
      b(j+1:n) = b(j+1:n) - a(j,j+1:n) * b(j)
    end do

  end if

  return
end
subroutine r8ut_vxm ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8UT_VXM multiplies an R8VEC by an R8UT matrix.
!
!  Discussion:
!
!    The R8UT storage format is used for an M by N upper triangular 
!    matrix.  The format stores all M*N entries, even those which are zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(M,N), the R8UT matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jhi
  real ( kind = 8 ) x(m)

  do i = 1, n
    jhi = min ( i, m )
    b(i) = sum ( x(1:jhi) * a(1:jhi,i) )
  end do

  return
end
subroutine r8vec_indicator ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR sets an R8VEC to the indicator vector.
!
!  Discussion:
!
!    A(1:N) = (/ 1 : N /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 1999
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

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) i, a(i)
  end do

  return
end
subroutine r8vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines 
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    if ( all ( a(1:n) == aint ( a(1:n) ) ) ) then
      do i = 1, n
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:n) ) < 1000000.0D+00 ) ) then
      do i = 1, n
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, n
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

  else if ( 3 <= max_print ) then

    if ( all ( a(1:max_print-2) == aint ( a(1:max_print-2) ) ) ) then
      do i = 1, max_print - 2
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-2) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print - 2
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print - 2
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

    write ( *, '(a)' ) '......  ..............'
    i = n

    if ( a(i) == real ( int ( a(i) ), kind = 8 ) ) then
      write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i8,2x,f14.6)' ) i, a(i)
    else
      write ( *, '(i8,2x,g14.6)' ) i, a(i)
    end if

  else

    if ( all ( a(1:max_print-1) == aint ( a(1:max_print-1) ) ) ) then
      do i = 1, max_print - 1
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-1) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print - 1
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print - 1
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

    i = max_print

    if ( a(i) == aint ( a(i) ) ) then
      write ( *, '(i8,2x,i8,a)' ) i, int ( a(i) ), '...more entries...'
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i8,2x,f14.6,a)' ) i, a(i), '...more entries...'
    else
      write ( *, '(i8,2x,g14.6,a)' ) i, a(i), '...more entries...'
    end if

  end if

  return
end
subroutine r8vec_read ( input_file, n, r )

!*****************************************************************************80
!
!! R8VEC_READ reads an R8VEC from a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file to be read.
!
!    Input, integer ( kind = 4 ) N, the size of the vector.
!
!    Output, real ( kind = 8 ) R(N), the vector.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) k
  real ( kind = 8 ) r(n)

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( input_file ) // '".'
    stop
  end if

  do k = 1, n

    read ( input_unit, *, iostat = ios ) r(k)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_READ - Fatal error!'
      write ( *, '(a,i8)' ) '  I/O error while reading record ', k
      stop
    end if

  end do

  close ( unit = input_unit )

  return
end
subroutine r8vec_read_size ( input_file, n )

!*****************************************************************************80
!
!! R8VEC_READ_SIZE reads the size of an R8VEC from a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file to 
!    be read.
!
!    Output, integer ( kind = 4 ) N, the size of the vector.
!
  implicit none

  character ( len = * ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) n
  real ( kind = 8 ) r

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_READ_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( input_file ) // '".'
    stop
  end if

  n = 0

  do

    read ( input_unit, *, iostat = ios ) r

    if ( ios /= 0 ) then
      exit
    end if

    n = n + 1

  end do

  close ( unit = input_unit )

  return
end
subroutine r8vec_to_r8cb ( m, n, ml, mu, x, a )

!*****************************************************************************80
!
!! R8VEC_TO_R8CB copies an R8VEC into an R8CB matrix.
!
!  Discussion:
!
!    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
!    a data item carries its dimensionality implicitly, and so cannot be
!    regarded sometimes as a vector and sometimes as an array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!
!    Input, real ( kind = 8 ) X((ML+MU+1)*N), the vector to be copied
!    into the array.
!
!    Output, real ( kind = 8 ) A(ML+MU+1,N), the array.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+mu+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x((ml+mu+1)*n)

  do j = 1, n
    do i = 1, ml + mu + 1

      if ( 1 <= i + j - mu - 1 .and. i + j - mu - 1 <= m ) then
        a(i,j) = x(i+(ml+mu+1)*(j-1))
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine r8vec_to_r8gb ( m, n, ml, mu, x, a )

!*****************************************************************************80
!
!! R8VEC_TO_R8GB copies an R8VEC into an R8GB matrix.
!
!  Discussion:
!
!    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
!    a data item carries its dimensionality implicitly, and so cannot be
!    regarded sometimes as a vector and sometimes as an array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!
!    Input, real ( kind = 8 ) X((2*ML+MU+1)*N), the vector to be copied
!    into the array.
!
!    Output, real ( kind = 8 ) A(2*ML+MU+1,N), the array.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x((2*ml+mu+1)*n)

  do j = 1, n
    do i = 1, 2 * ml + mu + 1

      if ( 1 <= i + j - ml - mu - 1 .and. i + j - ml - mu - 1 <= m ) then
        a(i,j) = x(i+(2*ml+mu+1)*(j-1))
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine r8vec_to_r8ge ( m, n, x, a )

!*****************************************************************************80
!
!! R8VEC_TO_R8GE copies an R8VEC into an R8GE matrix.
!
!  Discussion:
!
!    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
!    a data item carries its dimensionality implicitly, and so cannot be
!    regarded sometimes as a vector and sometimes as an array.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input, real ( kind = 8 ) X(M*N), the vector to be copied into the array.
!
!    Output, real ( kind = 8 ) A(M,N), the array.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(m*n)

  k = 0
  do j = 1, n
    do i = 1, m
      k = k + 1
      a(i,j) = x(k)
    end do
  end do

  return
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
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries 
!    in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, 
!    which should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
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
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_write ( n, r, output_file )

!*****************************************************************************80
!
!! R8VEC_WRITE writes an R8VEC to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) R(N), the vector to be written.
!
!    Input, character ( len = * ) OUTPUT_FILE, the name of the file to which
!    the information is to be written.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  character ( len = * ) output_file
  integer ( kind = 4 ) output_unit
  real ( kind = 8 ) r(n)

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' &
      // trim ( output_file ) // '".'
    stop
  end if

  do i = 1, n
    write ( output_unit, '(2x,g16.8)' ) r(i)
  end do

  close ( unit = output_unit )

  return
end
subroutine r8vec2_print_some ( n, x1, x2, max_print, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT_SOME prints "some" of a pair of R8VEC's.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vectors, is no more than MAX_PRINT, then
!    the entire vectors are printed, one entry of each per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vectors.
!
!    Input, real ( kind = 8 ) X1(N), X2(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines 
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title
  real ( kind = 8 ) x1(n)
  real ( kind = 8 ) x2(n)

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    write ( *, '(a)' ) '......  ..............  ..............'
    i = n
    write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    i = max_print
    write ( *, '(i8,2x,g14.6,2x,g14.6,2x,a)' ) i, x1(i), x2(i), &
      '...more entries...'

  end if

  return
end
subroutine r8vm_det ( n, a, det )

!*****************************************************************************80
!
!! R8VM_DET computes the determinant of an R8VM matrix.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8VM matrix.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  det = 1.0D+00
  do j = 1, n
    do i = j+1, n
      det = det * ( a(i) - a(j) )
    end do
  end do

  return
end
subroutine r8vm_mxv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8VM_MXV multiplies an R8VM matrix by an R8VEC.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8VM matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(M), the product A * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(n)

  do i = 1, m
    b(i) = 0.0D+00
    do j = 1, n
      if ( i == 1 ) then
        b(i) = b(i) + x(j)
      else
        b(i) = b(i) + a(j)**(i-1) * x(j)
      end if
    end do
  end do

  return
end
subroutine r8vm_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8VM_PRINT prints an R8VM matrix.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8VM matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) m
  character ( len = * ) title

  call r8vm_print_some ( m, n, a, 1, 1, n, n, title )

  return
end
subroutine r8vm_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8VM_PRINT_SOME prints some of an R8VM matrix.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8VM matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  logical r8_is_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) m
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i == 1 ) then
          aij = 1.0D+00
        else
          aij = a(j)**(i-1)
        end if

        if ( r8_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8vm_random ( m, n, seed, a )

!*****************************************************************************80
!
!! R8VM_RANDOM randomizes an R8VM matrix.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!    The parameter M is not actually needed by this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(N), the R8VM matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) seed

  call r8vec_uniform_01 ( n, seed, a )

  return
end
subroutine r8vm_sl ( n, a, b, x, job, info )

!*****************************************************************************80
!
!! R8VM_SL solves an R8VM linear system.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!    Vandermonde systems are very close to singularity.  The singularity
!    gets worse as N increases, and as any pair of values defining
!    the matrix get close.  Even a system as small as N = 10 will
!    involve the 9th power of the defining values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2003
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8VM matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, no error.
!    nonzero, at least two of the values in A are equal.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real ( kind = 8 ) x(n)
!
!  Check for explicit singularity.
!
  info = 0

  do j = 1, n - 1
    do i = j + 1, n
      if ( a(i) == a(j) ) then
        info = 1
        return
      end if
    end do
  end do

  x(1:n) = b(1:n)

  if ( job == 0 ) then

    do j = 1, n - 1
      do i = n, j + 1, -1
        x(i) = x(i) - a(j) * x(i-1)
      end do
    end do

    do j = n-1, 1, -1

      do i = j + 1, n
        x(i) = x(i) / ( a(i) - a(i-j) )
      end do

      do i = j, n - 1
        x(i) = x(i) - x(i+1)
      end do

    end do

  else

    do j = 1, n - 1
      do i = n, j + 1, -1
        x(i) = ( x(i) - x(i-1) ) / ( a(i) - a(i-j) )
      end do
    end do

    do j = n - 1, 1, -1
      do i = j, n - 1
        x(i) = x(i) - x(i+1) * a(j)
      end do
    end do

  end if

  return
end
subroutine r8vm_to_r8ge ( m, n, a, b )

!*****************************************************************************80
!
!! R8VM_TO_R8GE copies an R8VM matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8VM matrix.
!
!    Output, real ( kind = 8 ) B(M,N), the R8GE matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m
    do j = 1, n
      if ( i == 1 ) then
        b(i,j) = 1.0D+00
      else
        b(i,j) = b(i-1,j) * a(j)
      end if
    end do
  end do

  return
end
subroutine r8vm_vxm ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8VM_VXM multiplies an R8VEC by an R8VM matrix.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8VM matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m)

  do j = 1, n
    b(j) = 0.0D+00
    do i = 1, m
      if ( i == 1 ) then
        b(j) = b(j) + x(i)
      else
        b(j) = b(j) + a(j)**(i-1) * x(i)
      end if
    end do
  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I precedes J, ISGN = +1 if J precedes I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements 
!    I and J.  (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I precedes J;
!    ISGN => 0 means J precedes I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i = 2 * k1

    if ( i == n1 ) then
      j = k1
      k1 = i
      indx = - 1
      return
    else if ( i <= n1 ) then
      j = i + 1
      indx = - 2
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
  end if

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

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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
