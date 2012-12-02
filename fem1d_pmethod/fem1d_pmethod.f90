program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM1D_PMETHOD.
!
!  Discussion:
!
!    FEM1D_PMETHOD implements the P-version of the finite element method.
!
!    Program to solve the one dimensional problem:
!
!      - d/dX (P dU/dX) + Q U  =  F
!
!    by the finite-element method using a sequence of polynomials
!    which satisfy the boundary conditions and are orthogonal
!    with respect to the inner product:
!
!      (U,V)  =  Integral (-1 to 1) P U' V' + Q U V dx
!
!    Here U is an unknown scalar function of X defined on the
!    interval [-1,1], and P, Q and F are given functions of X.
!
!    The boundary values are U(-1) = U(1)=0.
!
!    Sample problem #1:
!
!      U=1-x**4,        P=1, Q=1, F=1.0+12.0*x**2-x**4
!
!    Sample problem #2:
!
!      U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x)
!
!    The program should be able to get the exact solution for
!    the first problem, using NP = 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    real ( kind = 8 ) A(0:NP), the squares of the norms of the
!    basis functions.
!
!    real ( kind = 8 ) ALPHA(NP).
!    ALPHA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    real ( kind = 8 ) BETA(NP).
!    BETA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    real ( kind = 8 ) F(0:NP).
!    F contains the basis function coefficients that form the
!    representation of the solution U.  That is,
!      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
!    where "BASIS(I)(X)" means the I-th basis function
!    evaluated at the point X.
!
!    integer NP.
!    The highest degree polynomial to use.
!
!    integer NPRINT.
!    The number of points at which the computed solution
!    should be printed out at the end of the computation.
!
!    integer PROBLEM, indicates the problem being solved.
!    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
!    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
!
!    integer QUAD_NUM, the order of the quadrature rule.
!
!    real ( kind = 8 ) QUAD_W(QUAD_NUM), the quadrature weights.
!
!    real ( kind = 8 ) QUAD_X(QUAD_NUM), the quadrature abscissas.
!
  implicit none

  integer, parameter :: np = 2
  integer, parameter :: quad_num = 10

  real ( kind = 8 ) a(0:np)
  real ( kind = 8 ) alpha(np)
  real ( kind = 8 ) beta(np)
  real ( kind = 8 ) f(0:np)
  integer, parameter :: nprint = 10
  integer, parameter :: problem = 2
  real ( kind = 8 ) quad_w(quad_num)
  real ( kind = 8 ) quad_x(quad_num)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_PMETHOD'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the two-point boundary value problem'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  - d/dX (P dU/dX) + Q U  =  F'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  on the interval [-1,1], with'
  write ( *, '(a)' ) '  U(-1) = U(1) = 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The P method is used, which represents U as'
  write ( *, '(a)' ) '  a weighted sum of orthogonal polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Highest degree polynomial to use is ', np
  write ( *, '(a,i8)' ) '  Number of points to be used for output = ', nprint

  if ( problem == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Problem #1:'
    write ( *, '(a)' ) '  U=1-x**4,'
    write ( *, '(a)' ) '  P=1,'
    write ( *, '(a)' ) '  Q=1,'
    write ( *, '(a)' ) '  F=1 + 12 * x**2 - x**4'
  else if ( problem == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Problem #2:'
    write ( *, '(a)' ) '  U=cos(0.5*pi*x),'
    write ( *, '(a)' ) '  P=1,'
    write ( *, '(a)' ) '  Q=0,'
    write ( *, '(a)' ) '  F=0.25*pi*pi*cos(0.5*pi*x)'
  end if
!
!  Get quadrature abscissas and weights for interval [-1,1].
!
  call quad ( quad_num, quad_w, quad_x )
!
!  Compute the constants for the recurrence relationship
!  that defines the basis functions.
!
  call alpbet ( a, alpha, beta, np, problem, quad_num, quad_w, quad_x )
!
!  Test the orthogonality of the basis functions.
!
  call ortho ( a, alpha, beta, np, problem, quad_num, quad_w, quad_x )
!
!  Solve for the solution of the problem, in terms of coefficients
!  of the basis functions.
!
  call sol ( a, alpha, beta, f, np, problem, quad_num, quad_w, quad_x )
!
!  Print out the solution, evaluated at each of the NPRINT points.
!
  call output ( alpha, beta, f, np, nprint )
!
!  Compare the computed and exact solutions.
!
  call compare ( alpha, beta, f, np, nprint, problem, quad_num, quad_w, quad_x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PMETHOD'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine alpbet ( a, alpha, beta, np, problem, quad_num, quad_w, quad_x )

!*****************************************************************************80
!
!! ALPBET calculates the coefficients in the recurrence relationship.
!
!  Discussion:
!
!    ALPHA and BETA are the coefficients in the three
!    term recurrence relation for the orthogonal basis functions
!    on [-1,1].
!
!    The routine also calculates A, the square of the norm of each basis
!    function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A(0:NP), the squares of the norms of the
!    basis functions.
!
!    Output, real ( kind = 8 ) ALPHA(NP).
!    ALPHA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Output, real ( kind = 8 ) BETA(NP).
!    BETA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Input, integer NP.
!    The highest degree polynomial to use.
!
!    Input, integer PROBLEM, indicates the problem being solved.
!    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
!    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
!
!    Input, integer QUAD_NUM, the order of the quadrature rule.
!
!    Input, real ( kind = 8 ) QUAD_W(QUAD_NUM), the quadrature weights.
!
!    Input, real ( kind = 8 ) QUAD_X(QUAD_NUM), the quadrature abscissas.
!
  implicit none

  integer np
  integer quad_num

  real ( kind = 8 ) a(0:np)
  real ( kind = 8 ) alpha(np)
  real ( kind = 8 ) beta(np)
  integer i
  integer iq
  integer k
  real ( kind = 8 ) pp
  integer problem
  real ( kind = 8 ) q
  real ( kind = 8 ) qm1
  real ( kind = 8 ) qm1x
  real ( kind = 8 ) qm2
  real ( kind = 8 ) qm2x
  real ( kind = 8 ) qq
  real ( kind = 8 ) quad_w(quad_num)
  real ( kind = 8 ) quad_x(quad_num)
  real ( kind = 8 ) qx
  real ( kind = 8 ) s
  real ( kind = 8 ) ss
  real ( kind = 8 ) su
  real ( kind = 8 ) sv
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x

  ss = 0.0D+00
  su = 0.0D+00

  do iq = 1, quad_num

    x = quad_x(iq)

    s = 4.0D+00 * pp ( x, problem ) * x * x &
      + qq ( x, problem ) * ( 1.0D+00 - x * x )**2

    u = 2.0D+00 * pp ( x, problem ) * x * ( 3.0D+00 * x * x - 1.0D+00 ) &
      + x * qq ( x, problem ) * ( 1.0D+00 - x * x )**2

    ss = ss + s * quad_w(iq)
    su = su + u * quad_w(iq)

  end do

  a(0) = ss
  alpha(1) = su / ss
  beta(1) = 0.0D+00

  do i = 2, np+1

    ss = 0.0D+00
    su = 0.0D+00
    sv = 0.0D+00

    do iq = 1, quad_num

      x = quad_x(iq)
      q = 1.0D+00
      qm1 = 0.0D+00
      qx = 0.0D+00
      qm1x = 0.0D+00

      do k = 1, i - 1
        qm2 = qm1
        qm1 = q
        qm2x = qm1x
        qm1x = qx
        q = ( x - alpha(k) ) * qm1 - beta(k) * qm2
        qx = qm1 + ( x - alpha(k) ) * qm1x - beta(k) * qm2x
      end do

      t = 1.0D+00 - x * x

      s = pp ( x, problem ) * ( t * qx - 2.0D+00 * x * q )**2 &
        + qq ( x, problem ) * t * t * q * q

      u = pp ( x, problem ) &
        * ( x * t * qx + ( 1.0D+00 - 3.0D+00 * x * x ) * q ) &
        * ( t * qx - 2.0D+00 * x * q ) + x * qq ( x, problem ) &
        * t * t * q * q

      v = pp ( x, problem ) &
        * ( x * t * qx + ( 1.0D+00 - 3.0D+00 * x * x ) * q ) &
        * ( t * qm1x - 2.0D+00 * x * qm1 ) &
        + x * qq ( x, problem ) * t * t * q * qm1

      ss = ss + s * quad_w(iq)
      su = su + u * quad_w(iq)
      sv = sv + v * quad_w(iq)

    end do

    a(i-1) = ss

    if ( i <= np ) then
      alpha(i) = su / ss
      beta(i) = sv / a(i-2)
    end if

  end do

  return
end
subroutine compare ( alpha, beta, f, np, nprint, problem, quad_num, &
  quad_w, quad_x )

!*****************************************************************************80
!
!! COMPARE compares the computed and exact solutions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA(NP).
!    ALPHA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Input, real ( kind = 8 ) BETA(NP).
!    BETA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Input, real ( kind = 8 ) F(0:NP).
!    F contains the basis function coefficients that form the
!    representation of the solution U.  That is,
!      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
!    where "BASIS(I)(X)" means the I-th basis function
!    evaluated at the point X.
!
!    Input, integer NP.
!    The highest degree polynomial to use.
!
!    Input, integer NPRINT.
!    The number of points at which the computed solution
!    should be printed out at the end of the computation.
!
!    Input, integer PROBLEM, indicates the problem being solved.
!    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
!    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
!
!    Input, integer QUAD_NUM, the order of the quadrature rule.
!
!    Input, real ( kind = 8 ) QUAD_W(QUAD_NUM), the quadrature weights.
!
!    Input, real ( kind = 8 ) QUAD_X(QUAD_NUM), the quadrature abscissas.
!
  implicit none

  integer np
  integer quad_num

  real ( kind = 8 ) alpha(np)
  real ( kind = 8 ) beta(np)
  real ( kind = 8 ) big_l2
  real ( kind = 8 ) error
  real ( kind = 8 ) f(0:np)
  integer i
  integer ip
  integer j
  integer k
  integer nprint
  integer, parameter :: nsub = 10
  real ( kind = 8 ) phii
  real ( kind = 8 ) phiix
  integer problem
  real ( kind = 8 ) quad_w(quad_num)
  real ( kind = 8 ) quad_x(quad_num)
  real ( kind = 8 ) ue
  real ( kind = 8 ) u_exact
  real ( kind = 8 ) up
  real ( kind = 8 ) x
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Comparison of computed and exact solutions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X        U computed    U exact     Difference'
  write ( *, '(a)' ) ' '

  do i = 0, nprint
    x = real ( 2 * i - nprint, kind = 8 ) / real ( nprint, kind = 8 )
    ue = u_exact ( x, problem )
    up = 0.0D+00
    do j = 0, np
      call phi ( alpha, beta, j, np, phii, phiix, x )
      up = up + phii * f(j)
    end do
    write(*,'(f8.4,3g14.6)') x, up, ue, ue - up
  end do
!
!  Compute the big L2 error.
!
  big_l2 = 0.0D+00

  do i = 1, nsub

    xl = real ( 2 * i - nsub - 1, kind = 8 ) / real ( nsub, kind = 8 )
    xr = real ( 2 * i - nsub,     kind = 8 ) / real ( nsub, kind = 8 )

    do j = 1, quad_num

      x = ( xl * ( 1.0D+00 - quad_x(j) ) &
          + xr * ( 1.0D+00 + quad_x(j) ) ) / 2.0D+00

      up = 0.0D+00
      do k = 0, np
        call phi ( alpha, beta, k, np, phii, phiix, x )
        up = up + phii * f(k)
      end do

      big_l2 = big_l2 + ( up - u_exact ( x, problem ) )**2 * quad_w(j) &
        * ( xr - xl ) / 2.0D+00

    end do

  end do

  big_l2 = sqrt ( big_l2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) 'Big L2 error = ', big_l2

  return
end
function ff ( x, problem )

!*****************************************************************************80
!
!! FF evaluates the right hand side function F(X) at any point X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Input, integer PROBLEM, indicates the problem being solved.
!    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
!    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
!
!    Output, real ( kind = 8 ) FF, the value of F(X).
!
  implicit none

  real ( kind = 8 ) ff
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer problem
  real ( kind = 8 ) x
!
!  Test problem 1
!
  if ( problem == 1 ) then

    ff = 1.0D+00 + 12.0D+00 * x**2 - x**4
!
!  Test problem 2
!
  else if ( problem == 2 ) then

    ff = 0.25D+00 * pi**2 * cos ( 0.5D+00 * pi * x )

  end if

  return
end
subroutine ortho ( a, alpha, beta, np, problem, quad_num, quad_w, quad_x )

!*****************************************************************************80
!
!! ORTHO tests the basis functions for orthogonality.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(0:NP), the squares of the norms of the
!    basis functions.
!
!    Input, real ( kind = 8 ) ALPHA(NP).
!    ALPHA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Input, real ( kind = 8 ) BETA(NP).
!    BETA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Input, integer NP.
!    The highest degree polynomial to use.
!
!    Input, integer PROBLEM, indicates the problem being solved.
!    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
!    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
!
!    Input, integer QUAD_NUM, the order of the quadrature rule.
!
!    Input, real ( kind = 8 ) QUAD_W(QUAD_NUM), the quadrature weights.
!
!    Input, real ( kind = 8 ) QUAD_X(QUAD_NUM), the quadrature abscissas.
!
  implicit none

  integer np
  integer quad_num

  real ( kind = 8 ) a(0:np)
  real ( kind = 8 ) alpha(np)
  real ( kind = 8 ) b(0:np,0:np)
  real ( kind = 8 ) beta(np)
  real ( kind = 8 ) bij
  integer i
  integer iq
  integer j
  real ( kind = 8 ) phii
  real ( kind = 8 ) phiix
  real ( kind = 8 ) phij
  real ( kind = 8 ) phijx
  real ( kind = 8 ) pp
  integer problem
  real ( kind = 8 ) qq
  real ( kind = 8 ) quad_w(quad_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) quad_x(quad_num)
!
!  Zero out the B array, so we can start summing up the dot products.
!
  b(0:np,0:np) = 0.0D+00
!
!  Approximate the integral of the product of basis function
!  I and basis function J over the interval [-1,1].
!
!  We expect to get zero, except when I and J are equal,
!  when we should get A(I).
!
  do iq = 1, quad_num
    x = quad_x(iq)
    do i = 0, np
      call phi ( alpha, beta, i, np, phii, phiix, x )
      do j = 0, np
        call phi ( alpha, beta, j, np, phij, phijx, x )
        bij = pp ( x, problem ) * phiix * phijx &
            + qq ( x, problem ) * phii * phij
        b(i,j) = b(i,j) + bij * quad_w(iq)
      end do
    end do
  end do
!
!  Print out the results of the test.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Basis function orthogonality test:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i   j     b(i,j)/a(i)'
  write ( *, '(a)' ) ' '
  do i = 0, np
    write ( *, '(a)' ) ' '
    do j = 0, np
      write ( *, '(2i4,g14.6)' ) i, j, b(i,j) / a(i)
    end do
  end do

  return
end
subroutine output ( alpha, beta, f, np, nprint )

!*****************************************************************************80
!
!! OUTPUT prints out the computed solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA(NP).
!    ALPHA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Input, real ( kind = 8 ) BETA(NP).
!    BETA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Input, real ( kind = 8 ) F(0:NP).
!    F contains the basis function coefficients that form the
!    representation of the solution U.  That is,
!      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
!    where "BASIS(I)(X)" means the I-th basis function
!    evaluated at the point X.
!
!    Input, integer NP.
!    The highest degree polynomial to use.
!
!    Input, integer NPRINT.
!    The number of points at which the computed solution
!    should be printed out at the end of the computation.
!
  implicit none

  integer np

  real ( kind = 8 ) alpha(np)
  real ( kind = 8 ) beta(np)
  real ( kind = 8 ) f(0:np)
  integer i
  integer ip
  integer nprint
  real ( kind = 8 ) phii
  real ( kind = 8 ) phiix
  real ( kind = 8 ) up
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Representation of solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Basis function coefficients:'
  write ( *, '(a)' ) ' '
  do i = 0, np
    write ( *, '(i4,g14.6)' ) i, f(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X     Approximate Solution'
  write ( *, '(a)' ) ' '
  do ip = 0, nprint
    x = real ( 2 * ip - nprint, kind = 8 ) / real ( nprint, kind = 8 )
    up = 0.0D+00
    do i = 0, np
      call phi ( alpha, beta, i, np, phii, phiix, x )
      up = up + phii * f(i)
    end do
    write ( *, '(2g14.6)' ) x, up
  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine phi ( alpha, beta, i, np, phii, phiix, x )

!*****************************************************************************80
!
!! PHI evaluates the I-th basis function at the point X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA(NP).
!    ALPHA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Input, real ( kind = 8 ) BETA(NP).
!    BETA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Input, integer I, the index of the basis function.
!
!    Input, integer NP.
!    The highest degree polynomial to use.
!
!    Output, real ( kind = 8 ) PHII, PHIIX, the value of the basis
!    function and its derivative.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
  implicit none

  integer np

  real ( kind = 8 ) alpha(np)
  real ( kind = 8 ) beta(np)
  integer i
  integer j
  real ( kind = 8 ) phii
  real ( kind = 8 ) phiix
  real ( kind = 8 ) q
  real ( kind = 8 ) qm1
  real ( kind = 8 ) qm1x
  real ( kind = 8 ) qm2
  real ( kind = 8 ) qm2x
  real ( kind = 8 ) qx
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  qm1 = 0.0D+00
  q = 1.0D+00
  qm1x = 0.0D+00
  qx = 0.0D+00

  do j = 1, i
    qm2 = qm1
    qm1 = q
    qm2x = qm1x
    qm1x = qx
    t = x - alpha(j)
    q = t * qm1 - beta(j) * qm2
    qx = qm1 + t * qm1x - beta(j) * qm2x
  end do

  t = 1.0D+00 - x**2
  phii = t * q
  phiix = t * qx - 2.0D+00 * x * q

  return
end
function pp ( x, problem )

!*****************************************************************************80
!
!! PP returns the value of the coefficient function P(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Input, integer PROBLEM, indicates the problem being solved.
!    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
!    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
!
!    Output, real ( kind = 8 ) PP, the value of P(X).
!
  implicit none

  real ( kind = 8 ) pp
  integer problem
  real ( kind = 8 ) x
!
!  Test problem 1
!
  if ( problem == 1 ) then

    pp = 1.0D+00
!
!  Test problem 2
!
  else if ( problem == 2 ) then

    pp = 1.0D+00

  end if

  return
end
subroutine quad ( quad_num, quad_w, quad_x )

!*****************************************************************************80
!
!! QUAD returns the abscissas and weights for gaussian quadrature on [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer QUAD_NUM, the order of the quadrature rule.
!
!    Output, real ( kind = 8 ) QUAD_W(QUAD_NUM), the quadrature weights.
!
!    Output, real ( kind = 8 ) QUAD_X(QUAD_NUM), the quadrature abscissas.
!
  implicit none

  integer quad_num

  real ( kind = 8 ) quad_w(quad_num)
  real ( kind = 8 ) quad_x(quad_num)
!
!  Quadrature points on [-1,1]
!
  quad_x(1) = -0.973906528517172D+00
  quad_x(2) = -0.865063366688985D+00
  quad_x(3) = -0.679409568299024D+00
  quad_x(4) = -0.433395394129247D+00
  quad_x(5) = -0.148874338981631D+00
  quad_x(6) =  0.148874338981631D+00
  quad_x(7) =  0.433395394129247D+00
  quad_x(8) =  0.679409568299024D+00
  quad_x(9) =  0.865063366688985D+00
  quad_x(10) = 0.973906528517172D+00
!
!  Weight factors
!
  quad_w(1) =  0.066671344308688D+00
  quad_w(2) =  0.149451349150581D+00
  quad_w(3) =  0.219086362515982D+00
  quad_w(4) =  0.269266719309996D+00
  quad_w(5) =  0.295524224714753D+00
  quad_w(6) =  0.295524224714753D+00
  quad_w(7) =  0.269266719309996D+00
  quad_w(8) =  0.219086362515982D+00
  quad_w(9) =  0.149451349150581D+00
  quad_w(10) = 0.066671344308688D+00

  return
end
function qq ( x, problem )

!*****************************************************************************80
!
!! QQ returns the value of the coefficient function Q(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Input, integer PROBLEM, indicates the problem being solved.
!    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
!    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
!
!    Output, real ( kind = 8 ) QQ, the value of Q(X).
!
  implicit none

  integer problem
  real ( kind = 8 ) qq
  real ( kind = 8 ) x
!
!  Test problem 1
!
  if ( problem == 1 ) then

    qq = 1.0D+00
!
!  Test problem 2
!
  else if ( problem == 2 ) then

    qq = 0.0D+00

  end if

  return
end
subroutine sol ( a, alpha, beta, f, np, problem, quad_num, quad_w, quad_x )

!*****************************************************************************80
!
!! SOL solves a linear system for the finite element coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(0:NP), the squares of the norms of the
!    basis functions.
!
!    Input, real ( kind = 8 ) ALPHA(NP).
!    ALPHA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Input, real ( kind = 8 ) BETA(NP).
!    BETA(I) contains one of the coefficients of a recurrence
!    relationship that defines the basis functions.
!
!    Output, real ( kind = 8 ) F(0:NP).
!    F contains the basis function coefficients that form the
!    representation of the solution U.  That is,
!      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
!    where "BASIS(I)(X)" means the I-th basis function
!    evaluated at the point X.
!
!    Input, integer NP.
!    The highest degree polynomial to use.
!
!    Input, integer PROBLEM, indicates the problem being solved.
!    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
!    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
!
!    Input, integer QUAD_NUM, the order of the quadrature rule.
!
!    Input, real ( kind = 8 ) QUAD_W(QUAD_NUM), the quadrature weights.
!
!    Input, real ( kind = 8 ) QUAD_X(QUAD_NUM), the quadrature abscissas.
!
  implicit none

  integer np
  integer quad_num

  real ( kind = 8 ) a(0:np)
  real ( kind = 8 ) alpha(np)
  real ( kind = 8 ) beta(np)
  real ( kind = 8 ) f(0:np)
  real ( kind = 8 ) ff
  integer i
  integer iq
  real ( kind = 8 ) phii
  real ( kind = 8 ) phiix
  integer problem
  real ( kind = 8 ) t
  real ( kind = 8 ) quad_w(quad_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) quad_x(quad_num)

  f(0:np) = 0.0D+00

  do iq = 1, quad_num
    x = quad_x(iq)
    t = ff ( x, problem ) * quad_w(iq)
    do i = 0, np
      call phi ( alpha, beta, i, np, phii, phiix, x )
      f(i) = f(i) + phii * t
    end do
  end do

  f(0:np) = f(0:np) / a(0:np)

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
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
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
function u_exact ( x, problem )

!*****************************************************************************80
!
!! U_EXACT returns the value of the exact solution at a point X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2006
!
!  Author:
!
!    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Input, integer PROBLEM, indicates the problem being solved.
!    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
!    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
!
!    Output, real ( kind = 8 ) U_EXACT, the exact value of U(X).
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer problem
  real ( kind = 8 ) u_exact
  real ( kind = 8 ) x
!
!  Test problem 1
!
  if ( problem == 1 ) then

    u_exact = 1.0D+00 - x**4
!
!  Test problem 2
!
  else if ( problem == 2 ) then

    u_exact = cos ( 0.5D+00 * pi * x )

  end if

  return
end
