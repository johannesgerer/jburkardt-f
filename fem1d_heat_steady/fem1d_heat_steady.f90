subroutine fem1d_heat_steady ( n, a, b, ua, ub, k, f, x, u )

!*****************************************************************************80
!
!! FEM1D_HEAT_STEADY solves the steady 1D heat equation with finite elements.
!
!  Discussion:
!
!    The program uses the finite element method, with piecewise linear basis
!    functions to solve the steady state heat equation in one dimension.
!
!    The problem is defined on the region A <= x <= B.
!
!    The following differential equation is imposed between A and B:
!
!      - d/dx k(x) du/dx = f(x)
!
!    where k(x) and f(x) are given functions.
!
!    At the boundaries, the following conditions are applied:
!
!      u(A) = UA
!      u(B) = UB
!
!    A set of N equally spaced nodes is defined on this
!    interval, with A = X(1) < X(2) < ... < X(N) = B.
!
!    At each node I, we associate a piecewise linear basis function V(I,X),
!    which is 0 at all nodes except node I.  This implies that V(I,X) is
!    everywhere 0 except that
!
!    for X(I-1) <= X <= X(I):
!
!      V(I,X) = ( X - X(I-1) ) / ( X(I) - X(I-1) ) 
!
!    for X(I) <= X <= X(I+1):
!
!      V(I,X) = ( X(I+1) - X ) / ( X(I+1) - X(I) )
!
!    We now assume that the solution U(X) can be written as a linear
!    sum of these basis functions:
!
!      U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X)
!
!    where U(X) on the left is the function of X, but on the right,
!    is meant to indicate the coefficients of the basis functions.
!
!    To determine the coefficient U(J), we multiply the original
!    differential equation by the basis function V(J,X), and use
!    integration by parts, to arrive at the I-th finite element equation:
!
!        Integral K(X) * U'(X) * V'(I,X) dx = Integral F(X) * V(I,X) dx
!
!    We note that the functions U(X) and U'(X) can be replaced by
!    the finite element form involving the linear sum of basis functions,
!    but we also note that the resulting integrand will only be nonzero
!    for terms where J = I - 1, I, or I + 1.
!
!    By writing this equation for basis functions I = 2 through N - 1,
!    and using the boundary conditions, we have N linear equations
!    for the N unknown coefficients U(1) through U(N), which can
!    be easily solved.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, real ( kind = 8 ) A, B, the left and right endpoints.
!
!    Input, real ( kind = 8 ) UA, UB, the prescribed value of U at A and B.
!
!    Input, external K, a function which evaluates k(x);
!
!    Input, external F, a function which evaluates f(x);
!
!    Input, real ( kind = 8 ) X(N), the mesh points.
!
!    Output, real ( kind = 8 ) U(N), the finite element coefficients, which 
!    are also the value of the computed solution at the mesh points.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: quad_num = 2

  real ( kind = 8 ) a
  real ( kind = 8 ) abscissa(quad_num)
  real ( kind = 8 ) al
  real ( kind = 8 ) am
  real ( kind = 8 ) ar
  real ( kind = 8 ) amat(n,n)
  real ( kind = 8 ) b
  real ( kind = 8 ) bm
  real ( kind = 8 ) bvec(n)
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fxq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  real ( kind = 8 ), external :: k
  real ( kind = 8 ) kxq
  integer ( kind = 4 ) q
  real ( kind = 8 ) weight(quad_num)
  real ( kind = 8 ) wq
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ua
  real ( kind = 8 ) ub
  real ( kind = 8 ) vl
  real ( kind = 8 ) vlp
  real ( kind = 8 ) vm
  real ( kind = 8 ) vmp
  real ( kind = 8 ) vr
  real ( kind = 8 ) vrp
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xm
  real ( kind = 8 ) xq
  real ( kind = 8 ) xr
!
!  Define a quadrature rule on the interval [-1,+1].
!
  abscissa(1) = -0.577350269189625764509148780502D+00
  abscissa(2) = +0.577350269189625764509148780502D+00
  weight(1) = 1.0D+00
  weight(2) = 1.0D+00
!
!  Zero out the matrix and right hand side.
!
  amat(1:n,1:n) = 0.0D+00
  bvec(1:n) = 0.0D+00
!
!  Equation 1 is the left boundary condition, U(A) = UA;
!
  amat(1,1) = 1.0D+00
  bvec(1) = ua
!
!  Equation I involves the basis function at node I.
!  This basis function is nonzero from X(I-1) to X(I+1).
!  Equation I looks like this:
!
!    Integral K(X) U'(X) V'(I,X) 
!           + C(X) * U(X) V(I,X) dx 
!  = Integral F(X) V(I,X) dx
!
!  Then, we realize that U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X), 
!  (U(X) means the function; U(J) is the coefficient of V(J,X) ).
!
!  The only V functions that are nonzero when V(I,X) is nonzero are
!  V(I-1,X) and V(I+1,X). 
!
!  Let's use the shorthand 
!
!    VL(X) = V(I-1,X)
!    VM(X) = V(I,X)
!    VR(X) = V(I+1,X)
!
!  So our equation becomes
!
!    Integral K(X) [ VL'(X) U(I-1) + VM'(X) U(I) + VR'(X) U(I+1) ] * VM'(X) dx
!  = Integral F(X) VM(X) dx.
!
!  
!
!  This is actually a set of N-2 linear equations for the N coefficients U.
!
!  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1), 
!  and so on.
!
  do i = 2, n - 1
!
!  Get the left, right and middle coordinates.
!
    xl = x(i-1)
    xm = x(i)
    xr = x(i+1)
!
!  Make temporary variables for A(I,I-1), A(I,I), A(I,I+1) and B(I).
!
    al = 0.0D+00
    am = 0.0D+00
    ar = 0.0D+00
    bm = 0.0D+00
!
!  We approximate the integrals by using a weighted sum of
!  the integrand values at quadrature points.
!
    do q = 1, quad_num
!
!  Integrate over the LEFT interval, between XL and XM, where:
!
!  VL(X) = ( XM - X       ) / ( XM - XL )
!  VM(X) = (      X  - XL ) / ( XM - XL )
!  VR(X) = 0
!
!  VL'(X) =             - 1 / ( XM - XL )
!  VM'(X) =             + 1 / ( XM - XL ) 
!  VR'(X) = 0
!
      xq = ( ( 1.0D+00 - abscissa(q) ) * xl   &
           + ( 1.0D+00 + abscissa(q) ) * xm ) &
           /   2.0D+00

      wq = weight(q) * ( xm - xl ) / 2.0D+00

      vl =  ( xm - xq ) / ( xm - xl )
      vlp =  - 1.0D+00  / ( xm - xl )

      vm =  ( xq - xl ) / ( xm - xl )
      vmp =  + 1.0D+00  / ( xm - xl )

      vr =  0.0D+00
      vrp = 0.0D+00

      kxq = k ( xq )
      fxq = f ( xq )

      al = al + wq * ( kxq * vlp * vmp )
      am = am + wq * ( kxq * vmp * vmp )
      ar = ar + wq * ( kxq * vrp * vmp )
      bm = bm + wq * ( fxq * vm )
!
!  Integrate over the RIGHT interval, between XM and XR, where:
!
!  VL(X) = 0
!  VM(X) = ( XR - X       ) / ( XR - XM )
!  VR(X) = (      X  - XM ) / ( XR - XM )
!
!  VL'(X) = 0
!  VM'(X) =             - 1 / ( XR - XM )
!  VR'(X) =             + 1 / ( XR - XM ) 
!
      xq = ( ( 1.0D+00 - abscissa(q) ) * xm   &
           + ( 1.0D+00 + abscissa(q) ) * xr ) &
           /   2.0D+00

      wq = weight(q) * ( xr - xm ) / 2.0D+00

      vl = 0.0D+00
      vlp = 0.0D+00

      vm = ( xr - xq ) / ( xr - xm )
      vmp = - 1.0D+00  / ( xr - xm )

      vr = ( xq - xm ) / ( xr - xm )
      vrp =  1.0D+00   / ( xr - xm )

      kxq = k ( xq )
      fxq = f ( xq )

      al = al + wq * ( kxq * vlp * vmp )
      am = am + wq * ( kxq * vmp * vmp )
      ar = ar + wq * ( kxq * vrp * vmp )
      bm = bm + wq * ( fxq * vm )

    end do

    amat(i,i-1) = al
    amat(i,i)   = am
    amat(i,i+1) = ar
    bvec(i)     = bm

  end do
!
!  Equation N is the right boundary condition, U(B) = UB;
!
  amat(n,n) = 1.0D+00
  bvec(n) = ub
!
!  Solve the linear system.
!
  call r8mat_print ( n, n, amat, '  Matrix A:' )

  call r8mat_solve2 ( n, amat, bvec, u, ierror )

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

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = * )  title

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
!    10 September 2009
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
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

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
subroutine r8mat_solve2 ( n, a, b, x, ierror )

!*****************************************************************************80
!
!! R8MAT_SOLVE2 computes the solution of an N by N linear system.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!    The linear system may be represented as
!
!      A*X = B
!
!    If the linear system is singular, but consistent, then the routine will
!    still produce a solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, A is the coefficient matrix to be inverted.
!    On output, A has been overwritten.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B is the right hand side of the system.
!    On output, B has been overwritten.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error detected.
!    1, consistent singularity.
!    2, inconsistent singularity.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) amax
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) ipiv(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)

  ierror = 0

  ipiv(1:n) = 0
  x(1:n) = 0.0D+00
!
!  Process the matrix.
!
  do k = 1, n
!
!  In column K:
!    Seek the row IMAX with the properties that:
!      IMAX has not already been used as a pivot;
!      A(IMAX,K) is larger in magnitude than any other candidate.
!
    amax = 0.0D+00
    imax = 0
    do i = 1, n
      if ( ipiv(i) == 0 ) then
        if ( amax < abs ( a(i,k) ) ) then
          imax = i
          amax = abs ( a(i,k) )
        end if
      end if
    end do
!
!  If you found a pivot row IMAX, then,
!    eliminate the K-th entry in all rows that have not been used for pivoting.
!
    if ( imax /= 0 ) then

      ipiv(imax) = k
      a(imax,k+1:n) = a(imax,k+1:n) / a(imax,k)
      b(imax) = b(imax) / a(imax,k)
      a(imax,k) = 1.0D+00

      do i = 1, n

        if ( ipiv(i) == 0 ) then
          a(i,k+1:n) = a(i,k+1:n) - a(i,k) * a(imax,k+1:n)
          b(i) = b(i) - a(i,k) * b(imax)
          a(i,k) = 0.0D+00
        end if

      end do

    end if

  end do
!
!  Now, every row with nonzero IPIV begins with a 1, and
!  all other rows are all zero.  Begin solution.
!
  do j = n, 1, -1

    imax = 0
    do k = 1, n
      if ( ipiv(k) == j ) then
        imax = k
      end if
    end do

    if ( imax == 0 ) then

      x(j) = 0.0D+00

      if ( b(j) == 0.0D+00 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_SOLVE2 - Warning:'
        write ( *, '(a,i8)' ) '  Consistent singularity, equation = ', j
      else
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_SOLVE2 - Error:'
        write ( *, '(a,i8)' ) '  Inconsistent singularity, equation = ', j
      end if

    else

      x(j) = b(imax)

      do i = 1, n
        if ( i /= imax ) then
          b(i) = b(i) - a(i,j) * x(j)
        end if
      end do

    end if

  end do

  return
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns an R8VEC of evenly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If N is 1, then the midpoint is returned.
!
!    Otherwise, the two endpoints are returned, and N-2 evenly
!    spaced points between them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) ALO, AHI, the low and high values.
!
!    Output, real ( kind = 8 ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) ahi
  real ( kind = 8 ) alo
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * alo   &
             + real (     i - 1, kind = 8 ) * ahi ) &
             / real ( n     - 1, kind = 8 )
    end do

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
