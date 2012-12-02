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
subroutine c83_cr_sl ( n, a_cr, b, x )

!*****************************************************************************80
!
!! C83_CR_SL solves a linear system factored by C83_CR_FA.
!
!  Discussion:
!
!    The matrix A must be tridiagonal.  C83_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.  If
!    the factors computed by C83_CR_FA are passed to C83_CR_SL, then a
!    linear system involving the matrix A may be solved.
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
!    Input, complex ( kind = 8 ) A_CR(3,0:2*N), factorization information 
!    computed by C83_CR_FA.
!
!    Input, real ( kind = 8 ) B(N), the right hand sides.
!
!    Output, real ( kind = 8 ) X(N), the solutions of the linear systems.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a_cr(3,0:2*n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) iful
  integer ( kind = 4 ) ifulm
  integer ( kind = 4 ) ihaf
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ipnt
  integer ( kind = 4 ) ipntp
  integer ( kind = 4 ) ndiv
  complex ( kind = 8 ) rhs(0:2*n)
  complex ( kind = 8 ) x(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C83_CR_SL - Fatal error!'
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

  rhs(ihaf) = rhs(ihaf) * a_cr(2,ihaf)
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
subroutine c83_cr_sls ( n, a_cr, nb, b, x )

!*****************************************************************************80
!
!! C83_CR_SLS solves several linear systems factored by C83_CR_FA.
!
!  Discussion:
!
!    The matrix A must be tridiagonal.  C83_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.  If
!    the factors computed by C83_CR_FA are passed to C83_CR_SLS, then one or
!    many linear systems involving the matrix A may be solved.
!
!    Note that C83_CR_FA does not perform pivoting, and so the solution 
!    produced by C83_CR_SLS may be less accurate than a solution produced 
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

  do j = 1, n - 1
    i = j + 1
    a(3,j) = cmplx ( i, j, kind = 8 )
  end do
  a(3,n) = 0.0D+00

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

        if ( 1 < i - j .or. 1 < j - i ) then

          crtemp(j2) = ' '
          citemp(j2) = ' '

        else

          if ( j == i - 1 ) then
            xr = real ( a(1,i), kind = 8 )
            xi = imag ( a(1,i) )
          else if ( j == i ) then
            xr = real ( a(2,i), kind = 8 )
            xi = imag ( a(2,i) )
          else if ( j == i + 1 ) then
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

      write ( *, '(i5,a,6a12)' ) i, ':', ( crtemp(j2), citemp(j2), j2 = 1, inc )

    end do

  end do

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

      write ( *, '(i5,a,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

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
    a(i) = cmplx ( i, - i, kind = 8 )
  end do

  return
end
subroutine c8vec_print ( n, a, title )

!*****************************************************************************80
!
!! C8VEC_PRINT prints a C8VEC, with an optional title.
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
    write ( *, '(i8,a,2g14.6)' ) i, ':', a(i)
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
      write ( *, '(i8,a,1x,2g14.6)' ) i, ':', x(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print-2
      write ( *, '(i8,a,1x,2g14.6)' ) i, ':', x(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i8,a,1x,2g14.6)' ) i, ':', x(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i8,a,1x,2g14.6)' ) i, ':', x(i)
    end do
    i = max_print
    write ( *, '(i8,a,1x,2g14.6,2x,a)' ) i, ':', x(i), '...more entries...'

  end if

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
!    Note that R83_CR_FA does not perform pivoting, and so the solutions 
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

  rhs(ihaf,1:nb) = a_cr(2,ihaf) * rhs(ihaf,1:nb)

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

  real      ( kind = 8 ) a(3,n)
  character ( len = * )  title

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

        if ( 1 < i - j .or. 1 < j - i ) then
          ctemp(j2) = '              '
        else if ( j == i + 1 ) then
          write ( ctemp(j2), '(g14.6)' ) a(1,j)
        else if ( j == i ) then
          write ( ctemp(j2), '(g14.6)' ) a(2,j)
        else if ( j == i - 1 ) then
          write ( ctemp(j2), '(g14.6)' ) a(3,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2004
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
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 March 2005
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
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

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
!! R8VEC_PRINT prints an R8VEC, with an optional title.
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
    write ( *, '(i8,a,g14.6)' ) i, ':', a(i)
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
        write ( *, '(i8,a,1x,i8)' ) i, ':', int ( a(i) )
      end do
    else if ( all ( abs ( a(1:n) ) < 1000000.0D+00 ) ) then
      do i = 1, n
        write ( *, '(i8,a,1x,f14.6)' ) i, ':', a(i)
      end do
    else
      do i = 1, n
        write ( *, '(i8,a,1x,g14.6)' ) i, ':', a(i)
      end do
    end if

  else if ( 3 <= max_print ) then

    if ( all ( a(1:max_print-2) == aint ( a(1:max_print-2) ) ) ) then
      do i = 1, max_print - 2
        write ( *, '(i8,a,1x,i8)' ) i, ':', int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-2) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print - 2
        write ( *, '(i8,a,1x,f14.6)' ) i, ':', a(i)
      end do
    else
      do i = 1, max_print - 2
        write ( *, '(i8,a,1x,g14.6)' ) i, ':', a(i)
      end do
    end if

    write ( *, '(a)' ) '......  ..............'
    i = n

    if ( a(i) == real ( int ( a(i) ), kind = 8 ) ) then
      write ( *, '(i8,a,1x,i8)' ) i, ':', int ( a(i) )
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i8,a,1x,f14.6)' ) i, ':', a(i)
    else
      write ( *, '(i8,a,1x,g14.6)' ) i, ':', a(i)
    end if

  else

    if ( all ( a(1:max_print-1) == aint ( a(1:max_print-1) ) ) ) then
      do i = 1, max_print - 1
        write ( *, '(i8,a,1x,i8)' ) i, ':', int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-1) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print - 1
        write ( *, '(i8,a,1x,f14.6)' ) i, ':', a(i)
      end do
    else
      do i = 1, max_print - 1
        write ( *, '(i8,a,1x,g14.6)' ) i, ':', a(i)
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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
