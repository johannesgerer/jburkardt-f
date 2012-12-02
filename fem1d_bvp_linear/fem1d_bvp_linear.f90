subroutine compute_l2_error ( n, x, u, exact, l2_error )

!*****************************************************************************80
!
!! COMPUTE_L2_ERROR estimates the L2 error norm of a finite element solution.
!
!  Discussion:
!
!    We assume the finite element method has been used, over an interval [A,B]
!    involving N nodes, with piecewise linear elements used for the basis.
!    The coefficients U(1:N) have been computed, and a formula for the
!    exact solution is known.
!
!    This function estimates the L2 norm of the error:
!
!      L2_NORM = Integral ( A <= X <= B ) ( U(X) - EXACT(X) )^2 dX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, real ( kind = 8 ) X(N), the mesh points.
!
!    Input, real ( kind = 8 ) U(N), the finite element coefficients.
!
!    Input, function EQ = EXACT ( X ), returns the value of the exact
!    solution at the point X.
!
!    Output, real ( kind = 8 ) L2_ERROR, the estimated L2 norm of the error.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) abscissa(2)
  real ( kind = 8 ) eq
  real ( kind = 8 ), external :: exact
  integer ( kind = 4 ) i
  real ( kind = 8 ) l2_error
  integer ( kind = 4 ) q
  integer ( kind = 4 ) quad_num
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) uq
  real ( kind = 8 ) weight(2)
  real ( kind = 8 ) wq
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xq
  real ( kind = 8 ) xr

  l2_error = 0.0D+00
!
!  Quadrature definitions.
!
  quad_num = 2
  abscissa(1) = -0.577350269189625764509148780502D+00
  abscissa(2) = +0.577350269189625764509148780502D+00
  weight(1) = 1.0D+00
  weight(2) = 1.0D+00
!
!  Integrate over each interval.
!
  do i = 1, n - 1

    xl = x(i)
    xr = x(i+1)
    ul = u(i)
    ur = u(i+1)

    do q = 1, quad_num

      xq = ( ( 1.0D+00 - abscissa(q) ) * xl   &
           + ( 1.0D+00 + abscissa(q) ) * xr ) &
           /   2.0D+00

      wq = weight(q) * ( xr - xl ) / 2.0D+00
!
!  Use the fact that U is a linear combination of piecewise linears.
!
      uq = ( ( xr - xq      ) * ul &
           + (      xq - xl ) * ur ) &
           / ( xr      - xl )

      eq = exact ( xq )

      l2_error = l2_error + wq * ( uq - eq )**2

    end do

  end do

  l2_error = sqrt ( l2_error )

  return
end
subroutine compute_seminorm_error ( n, x, u, exact_ux, seminorm_error )

!*****************************************************************************80
!
!! COMPUTE_SEMINORM_ERROR estimates the seminorm error of a finite element solution.
!
!  Discussion:
!
!    We assume the finite element method has been used, over an interval [A,B]
!    involving N nodes, with piecewise linear elements used for the basis.
!    The coefficients U(1:N) have been computed, and a formula for the
!    exact derivative is known.
!
!    This function estimates the seminorm of the error:
!
!      SEMINORM = Integral ( A <= X <= B ) ( dU(X)/dx - EXACT_UX(X) )^2 dX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, real ( kind = 8 ) X(N), the mesh points.
!
!    Input, real ( kind = 8 ) U(N), the finite element coefficients.
!
!    Input, function EQ = EXACT_UX ( X ), returns the value of the exact
!    derivative at the point X.
!
!    Output, real ( kind = 8 ) SEMINORM_ERROR, the estimated seminorm of 
!    the error.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) abscissa(2)
  real ( kind = 8 ), external :: exact_ux
  real ( kind = 8 ) exq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) q
  integer ( kind = 4 ) quad_num
  real ( kind = 8 ) seminorm_error
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) uxq
  real ( kind = 8 ) weight(2)
  real ( kind = 8 ) wq
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xq
  real ( kind = 8 ) xr

  seminorm_error = 0.0D+00
!
!  Quadrature definitions.
!
  quad_num = 2
  abscissa(1) = -0.577350269189625764509148780502D+00
  abscissa(2) = +0.577350269189625764509148780502D+00
  weight(1) = 1.0D+00
  weight(2) = 1.0D+00
!
!  Integrate over each interval.
!
  do i = 1, n - 1

    xl = x(i)
    xr = x(i+1)
    ul = u(i)
    ur = u(i+1)

    do q = 1, quad_num

      xq = ( ( 1.0D+00 - abscissa(q) ) * xl   &
           + ( 1.0D+00 + abscissa(q) ) * xr ) &
           /   2.0D+00

      wq = weight(q) * ( xr - xl ) / 2.0D+00
!
!  The piecewise linear derivative is a constant in the interval.
!
      uxq = ( ur - ul ) / ( xr - xl )

      exq = exact_ux ( xq )
 
      seminorm_error = seminorm_error + wq * ( uxq - exq )**2

    end do

  end do

  seminorm_error = sqrt ( seminorm_error )

  return
end
subroutine fem1d_bvp_linear ( n, a, c, f, x, u )

!*****************************************************************************80
!
!! FEM1D_BVP_LINEAR solves a two point boundary value problem.
!
!  Discussion:
!
!    The program uses the finite element method, with piecewise linear basis
!    functions to solve a boundary value problem in one dimension.
!
!    The problem is defined on the region 0 <= x <= 1.
!
!    The following differential equation is imposed between 0 and 1:
!
!      - d/dx a(x) du/dx + c(x) * u(x) = f(x)
!
!    where a(x), c(x), and f(x) are given functions.
!
!    At the boundaries, the following conditions are applied:
!
!      u(0.0) = 0.0
!      u(1.0) = 0.0
!
!    A set of N equally spaced nodes is defined on this
!    interval, with 0 = X(1) < X(2) < ... < X(N) = 1.0.
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
!        Integral A(X) * U'(X) * V'(I,X) + C(X) * U(X) * V(I,X) dx 
!      = Integral F(X) * V(I,X) dx
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
!    20 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, function A ( X ), evaluates a(x);
!
!    Input, function C ( X ), evaluates c(x);
!
!    Input, function F ( X ), evaluates f(x);
!
!    Input, real ( kind = 8 ) X(N), the mesh points.
!
!    Output, real ( kind = 8 ) U(N), the finite element coefficients, which 
!    are also the value of the computed solution at the mesh points.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: quad_num = 2

  real ( kind = 8 ), external :: a
  real ( kind = 8 ) abscissa(quad_num)
  real ( kind = 8 ) al
  real ( kind = 8 ) am
  real ( kind = 8 ) ar
  real ( kind = 8 ) amat(n,n)
  real ( kind = 8 ) axq
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) bm
  real ( kind = 8 ), external :: c
  real ( kind = 8 ) cxq
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fxq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) q
  real ( kind = 8 ) weight(quad_num)
  real ( kind = 8 ) wq
  real ( kind = 8 ) u(n)
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
!  Quadrature definitions.
!
  abscissa(1) = -0.577350269189625764509148780502D+00
  abscissa(2) = +0.577350269189625764509148780502D+00
  weight(1) = 1.0D+00
  weight(2) = 1.0D+00
!
!  Zero out the matrix and right hand side.
!
  amat(1:n,1:n) = 0.0D+00
  b(1:n) = 0.0D+00
!
!  Equation 1 is the left boundary condition, U(0.0) = 0.0;
!
  amat(1,1) = 1.0D+00
  b(1) = 0.0D+00
!
!  Equation I involves the basis function at node I.
!  This basis function is nonzero from X(I-1) to X(I+1).
!  Equation I looks like this:
!
!    Integral A(X) U'(X) V'(I,X) 
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
!    Integral A(X) [ VL'(X) U(I-1) + VM'(X) U(I) + VR'(X) U(I+1) ] * VM'(X)
!           + C(X) [ VL(X)  U(I-1) + VM(X)  U(I) + VR(X)  U(I+1) ] * VM(X) dx
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

      axq = a ( xq )
      cxq = c ( xq )
      fxq = f ( xq )

      al = al + wq * ( axq * vlp * vmp + cxq * vl * vm )
      am = am + wq * ( axq * vmp * vmp + cxq * vm * vm )
      ar = ar + wq * ( axq * vrp * vmp + cxq * vr * vm )
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

      axq = a ( xq )
      cxq = c ( xq )
      fxq = f ( xq )

      al = al + wq * ( axq * vlp * vmp + cxq * vl * vm )
      am = am + wq * ( axq * vmp * vmp + cxq * vm * vm )
      ar = ar + wq * ( axq * vrp * vmp + cxq * vr * vm )
      bm = bm + wq * ( fxq * vm )

    end do

    amat(i,i-1) = al
    amat(i,i)   = am
    amat(i,i+1) = ar
    b(i)        = bm

  end do
!
!  Equation N is the right boundary condition, U(1.0) = 0.0;
!
  amat(n,n) = 1.0D+00
  b(n) = 0.0D+00
!
!  Solve the linear system.
!
  call r8mat_solve2 ( n, amat, b, u, ierror )

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
