program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM1D_ADAPTIVE.
!
!  Discussion:
!
!    FEM1D_ADAPTIVE solves a 1D problem using an adaptive finite element method.
!
!    The equation to be treated is:
!
!      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
!
!    by the finite-element method using piecewise linear basis
!    functions.
!
!    An adaptive method is used to try to reduce the maximum
!    error by refining the mesh in certain places.
!
!    Here U is an unknown scalar function of X defined on the
!    interval [XL,XR], and P, Q and F are given functions of X.
!
!    The values of U at XL and XR are also specified.
!
!    The interval [XL,XR] is "meshed" with N+1 points,
!
!      XN(0) = XL, XN(1) = XL+H, XN(2) = XL+2*H, ..., XN(N) = XR.
!
!    This creates N subintervals, with interval I having endpoints
!    XN(I-1) and XN(I).
!
!
!    The algorithm tries to guarantee a certain amount
!    of accuracy by examining the current solution, estimating the error
!    in each subinterval, and, if necessary, subdividing one or more
!    subintervals and repeating the calculation.
!
!    We can think of the adaptive part of the algorithm as a refined
!    problem.  The program re-solves the problem on the pair of
!    intervals J and J+1, which extend from node J-1 to node J+1.
!    The values of U that were just computed at nodes J-1 and J+1
!    will be used as the boundary values for this refined problem.
!    The intervals J and J+1 will each be evenly divided into NY
!    smaller subintervals.  This boundary value problem is solved,
!    and the derivatives of the original and refined solutions are
!    then compared to get an estimate of the error.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) ADIAG(NU).
!    ADIAG(I) is the "diagonal" coefficient of the I-th
!    equation in the linear system.  That is, ADIAG(I) is
!    the coefficient of the I-th unknown in the I-th equation.
!
!    real ( kind = 8 ) ALEFT(NU).
!    ALEFT(I) is the "left hand" coefficient of the I-th
!    equation in the linear system.  That is, ALEFT(I) is the
!    coefficient of the (I-1)-th unknown in the I-th equation.
!    There is no value in ALEFT(1), since the first equation
!    does not refer to a "0-th" unknown.
!
!    real ( kind = 8 ) ARITE(NU).
!    ARITE(I) is the "right hand" coefficient of the I-th
!    equation in the linear system.  ARITE(I) is the coefficient
!    of the (I+1)-th unknown in the I-th equation.  There is
!    no value in ARITE(NU) because the NU-th equation does not
!    refer to an "NU+1"-th unknown.
!
!    real ( kind = 8 ) ETA(N).
!    ETA(I) is the error estimate for interval I.  It is computed
!    as the sum of two quantities, one associated with the left
!    and one with the right node of the interval.
!
!    real ( kind = 8 ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    real ( kind = 8 ) FY(M).
!    FY is the right hand side of the linear system of the refined
!    problem.
!
!    real ( kind = 8 ) H(N)
!    H(I) is the length of subinterval I.  This code uses
!    equal spacing for all the subintervals.
!
!    real ( kind = 8 ) HY(M).
!    HY(I) is the length of subinterval I in the refined problem.
!
!    integer IBC.
!    IBC declares what the boundary conditions are.
!    1, at the left endpoint, U has the value UL,
!       at the right endpoint, U' has the value UR.
!    2, at the left endpoint, U' has the value UL,
!       at the right endpoint, U has the value UR.
!    3, at the left endpoint, U has the value UL,
!       and at the right endpoint, U has the value UR.
!    4, at the left endpoint, U' has the value UL,
!       at the right endpoint U' has the value UR.
!
!    integer IBCY.
!    IBCY declares the boundary conditions for the refined problem
!    which should always be that the value of U is specified at
!    both the left and right endpoints.  This corresponds to a
!    value of IBCY = 3.
!
!    integer INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    integer INDY(0:M).
!    INDY(I) records the index of the unknown associated with
!    node I for the refined problem.
!
!    integer JADD(N).
!    JADD(I) is 1 if the error estimates show that interval I
!    should be subdivided.
!
!    integer KOUNT, the number of adaptive steps that have been taken.
!
!    integer M.
!    M is the number of subintervals used in the refined problem.
!    M is equal to NY for computations centered at node 0 or node N,
!    and otherwise, M is equal to 2*NY.
!
!    integer N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    integer NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    integer NMAX, the maximum number of unknowns that can be handled.
!
!    integer NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    integer NODEY(NL,M).
!    NODEY performs the same function for the refined problem that
!    NODE performs for the full problem, recording the node numbers
!    associated with a particular subinterval.
!
!    integer NQUAD
!    The number of quadrature points used in a subinterval.
!
!    integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    integer NUY.
!    The number of unknowns in the refined problem.
!
!    integer NY.
!    NY is the number of subintervals into which a given interval
!    will be subdivided, before solving the refined probelm.
!
!    integer PROBLEM, chooses the problem to be solved.
!    The user must choose this value by setting it in routine GETPRB.
!    * 1, u = x, p = 1, q = 0, f = 0, ibc = 3, ul = 0, ur = 1.
!    The program should find the solution exactly, and the
!    adaptive code should find that there is no reason to
!    subdivide any interval.
!    * 2, u = x*x, p = 1, q = 0, f = -2, ibc = 3, ul = 0, ur = 1.
!    This problem should find the solution exactly, and
!    the adaptive code should again find there is nothing
!    to do.
!    *3, u = sin(pi*x/2), p = 1, q = 0, ibc = 3, f = 0.25*pi*pi*sin(pi*x/2),
!    ul = 0, ur = 1.
!    *4, u = cos(pi*x/2), p = 1, q = 0, ibc = 3, f = 0.25*pi*pi*cos(pi*x/2),
!    ul = 1, ur = 0.
!    *5: u = x**(beta+2)/((beta+2)*(beta+1)), p = 1, q = 1, ibc = 3,
!    f = -x**beta + (x**(beta+2))/((beta+2)*(beta+1)),
!    ul = 0, ur = 1/((beta+2)*(beta+1))
!    (beta must be greater than -2, and not equal to -1)
!    *6: u = atan((x-0.5)/alpha), p = 1, q = 0, ibc = 3,
!    f =  2*alpha*(x-0.5) / (alpha**2 + (x-0.5)**2) **2,
!    ul = u(0), ur = u(1)
!
!    integer STATUS, reports status of subdivision.
!    0, a new subdivision was carried out.
!    1, no more subdivisions are needed.
!    -1, no more subdivisions can be carried out.
!
!    real ( kind = 8 ) TOL.
!    A tolerance that is used to determine whether the estimated
!    error in an interval is so large that it should be subdivided
!    and the problem solved again.
!
!    real ( kind = 8 ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    real ( kind = 8 ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    real ( kind = 8 ) WQUAD(NQUAD).
!    WQUAD(I) is the weight associated with the I-th point
!    of an NQUAD point Gaussian quadrature rule.
!
!    real ( kind = 8 ) XL.
!    XL is the left endpoint of the interval over which the
!    differential equation is being solved.
!
!    real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    real ( kind = 8 ) XQUAD(NQUAD,NMAX), the I-th quadrature point
!    in interval J.
!
!    real ( kind = 8 ) XQUADY(NQUAD,NMAY ), the I-th quadrature point
!    in subinterval J of the refined problem.
!
!    real ( kind = 8 ) XR.
!    XR is the right endpoint of the interval over which the
!    differential equation is being solved.
!
!    Workspace, double precision XT(0:NMAX), used to compute a new
!    set of nodes.
!
!    real ( kind = 8 ) YN(0:M).
!    YN(I) is the location of the I-th node in the refined
!    problem.
!
  implicit none

  integer, parameter :: nl = 2
  integer, parameter :: nmax = 30
  integer, parameter :: nquad = 2

  real ( kind = 8 ) adiag(nmax)
  real ( kind = 8 ) aleft(nmax)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) arite(nmax)
  real ( kind = 8 ) beta
  real ( kind = 8 ) eta(nmax)
  real ( kind = 8 ) f(nmax)
  real ( kind = 8 ) h(nmax)
  integer ibc
  integer indx(0:nmax)
  integer jadd(nmax)
  integer kount
  integer n
  integer node(nl,nmax)
  integer nu
  integer problem
  integer status
  real ( kind = 8 ) tol
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) wquad(nquad)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xn(0:nmax)
  real ( kind = 8 ) xquad(nquad,nmax)
  real ( kind = 8 ) xr
  real ( kind = 8 ) xt(0:nmax)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_ADAPTIVE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Solve the two-point boundary value problem:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'on the interval [0,1], specifying the value'
  write ( *, '(a)' ) 'of U at each endpoint.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) 'The number of basis functions per element is ', nl
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) 'The number of quadrature points per element is ', nquad

  call get_problem ( problem )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Problem index = ', problem
  write ( *, '(a)' ) ' '

  if ( problem == 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "Linear" problem:'
    write ( *, '(a)' ) '  (No refinement needed)'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  X'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  0.0'
    write ( *, '(a)' ) '  F(X) =  0.0'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  0.0'
    write ( *, '(a)' ) '  UR   =  1.0'

  else if ( problem == 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "Quadratic" problem:'
    write ( *, '(a)' ) '  (No refinement needed)'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  X*X'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  0.0'
    write ( *, '(a)' ) '  F(X) = -2.0'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  0.0'
    write ( *, '(a)' ) '  UR   =  1.0'

  else if ( problem == 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "SINE" problem:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  SIN(PI*X/2)'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  0.0'
    write ( *, '(a)' ) '  F(X) =  PI*PI*SIN(PI*X/2)/4'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  0.0'
    write ( *, '(a)' ) '  UR   =  1.0'

  else if ( problem == 4 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "COSINE" problem:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  COS(PI*X/2)'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  0.0'
    write ( *, '(a)' ) '  F(X) =  PI*PI*COS(PI*X/2)/4'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  0.0'
    write ( *, '(a)' ) '  UR   =  1.0'

  else if ( problem == 5 ) then

    call get_beta ( beta )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "RHEINBOLDT" problem:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  X**(B+2)/((B+2)*(B+1))'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  1.0'
    write ( *, '(a)' ) '  F(X) =  -X**B+(X**B+2))/((B+2)*(B+1))'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  0.0'
    write ( *, '(a)' ) '  UR   =  1/((B+2)*(B+1))'
    write ( *, '(a,g14.6)' ) '  B    = ', beta

  else if ( problem == 6 ) then

    call get_alpha ( alpha )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "ARCTAN" problem:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  U(X) =  ATAN((X-0.5)/A)'
    write ( *, '(a)' ) '  P(X) =  1.0'
    write ( *, '(a)' ) '  Q(X) =  0.0'
    write ( *, '(a)' ) '  F(X) =  2*A*(X-0.5)/(A**2+(X-0.5)**2)**2'
    write ( *, '(a)' ) '  IBC  =  3'
    write ( *, '(a)' ) '  UL   =  ATAN(-0.5/A)'
    write ( *, '(a)' ) '  UR   =  ATAN( 0.5/A)'
    write ( *, '(a,g14.6)' ) '  A    = ', alpha

  end if
!
!  Start out with just 4 subintervals.
!
  n = 4
!
!  Initialize values that define the problem.
!
  call init ( ibc, n, tol, ul, ur, xl, xn, xr )
!
!  Start the iteration counter off at 0.
!
  kount = 0
!
!  Begin the next iteration.
!
  do

    kount = kount + 1

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) 'Begin new iteration with ', n, ' nodes.'
    write ( *, '(a)' ) ' '
!
!  Solve the regular problem.
!
    call solvex ( f, h, ibc, kount, n, nl, nmax, nu, ul, ur, xn )
!
!  Solve N subproblems to get the error estimators.
!
    call solvey ( eta, f, h, n, nu, ul, ur, xn )
!
!  Examine the error estimators, and see how many intervals should
!  be subdivided.
!
    call subdiv ( eta, kount, n, nmax, tol, xn, status )

    if ( status /= 0 ) then
      exit
    end if
!
!  Solve the problem again, with the new nodes.
!
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_ADAPTIVE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine assemble ( adiag, aleft, arite, f, h, n, indx, node, nu, nl, &
  nquad, nmax, ul, ur, wquad, xn, xquad )

!*****************************************************************************80
!
!! ASSEMBLE assembles the global matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ADIAG(NU).
!    ADIAG(I) is the "diagonal" coefficient of the I-th
!    equation in the linear system.  That is, ADIAG(I) is
!    the coefficient of the I-th unknown in the I-th equation.
!
!    Output, real ( kind = 8 ) ALEFT(NU).
!    ALEFT(I) is the "left hand" coefficient of the I-th
!    equation in the linear system.  That is, ALEFT(I) is the
!    coefficient of the (I-1)-th unknown in the I-th equation.
!    There is no value in ALEFT(1), since the first equation
!    does not refer to a "0-th" unknown.
!
!    Output, real ( kind = 8 ) ARITE(NU).
!    ARITE(I) is the "right hand" coefficient of the I-th
!    equation in the linear system.  ARITE(I) is the coefficient
!    of the (I+1)-th unknown in the I-th equation.  There is
!    no value in ARITE(NU) because the NU-th equation does not
!    refer to an "NU+1"-th unknown.
!
!    Output, real ( kind = 8 ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    Input, real ( kind = 8 ) H(N)
!    H(I) is the length of subinterval I.  This code uses
!    equal spacing for all the subintervals.
!
!    Input, integer N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    Input, integer NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    Input, integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Input, integer NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    Input, integer NQUAD
!    The number of quadrature points used in a subinterval.
!
!    Input, integer NMAX, the maximum number of unknowns that can be handled.
!
!    Input, real ( kind = 8 ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Input, real ( kind = 8 ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Input, real ( kind = 8 ) WQUAD(NQUAD).
!    WQUAD(I) is the weight associated with the I-th point
!    of an NQUAD point Gaussian quadrature rule.
!
!    Input, real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    Input, real ( kind = 8 ) XQUAD(NQUAD,NMAX), the I-th quadrature point
!    in interval J.
!
  implicit none

  integer n
  integer nl
  integer nmax
  integer nquad
  integer nu

  real ( kind = 8 ) adiag(nu)
  real ( kind = 8 ) aleft(nu)
  real ( kind = 8 ) arite(nu)
  real ( kind = 8 ) aij
  real ( kind = 8 ) f(nu)
  real ( kind = 8 ), external :: ff
  real ( kind = 8 ) h(n)
  real ( kind = 8 ) he
  integer i
  integer ie
  integer ig
  integer il
  integer indx(0:n)
  integer iq
  integer iu
  integer jg
  integer jl
  integer ju
  integer node(nl,nmax)
  real ( kind = 8 ) phii
  real ( kind = 8 ) phiix
  real ( kind = 8 ) phij
  real ( kind = 8 ) phijx
  real ( kind = 8 ), external :: pp
  real ( kind = 8 ), external :: qq
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) wquad(nquad)
  real ( kind = 8 ) wquade
  real ( kind = 8 ) x
  real ( kind = 8 ) xleft
  real ( kind = 8 ) xn(0:n)
  real ( kind = 8 ) xquad(nquad,nmax)
  real ( kind = 8 ) xquade
  real ( kind = 8 ) xrite
!
!  Zero out the entries.
!
  f(1:nu) = 0.0D+00
  aleft(1:nu) = 0.0D+00
  arite(1:nu) = 0.0D+00
  adiag(1:nu) = 0.0D+00
!
!  For each interval,
!
  do ie = 1, n

    he = h(ie)
    xleft = xn(node(1,ie))
    xrite = xn(node(2,ie))
!
!  For each quadrature point in the interval,
!
    do iq = 1, nquad

      xquade = xquad(iq,ie)
      wquade = wquad(iq)
!
!  Pick a basis function which defines the equation,
!
      do il = 1, nl

        ig = node(il,ie)
        iu = indx(ig)

        if ( 0 < iu ) then

          call phi ( il, xquade, phii, phiix, xleft, xrite )
          f(iu) = f(iu) + he * wquade * ff ( xquade ) * phii
!
!  Take care of boundary conditions specifying the value of U'.
!
          if ( ig == 0 ) then
            x = 0.0D+00
            f(iu) = f(iu) - pp ( x ) * ul
          else if ( ig == n ) then
            x = 1.0D+00
            f(iu) = f(iu) + pp ( x ) * ur
          end if
!
!  Pick a basis function which defines the coefficient
!  being computed.
!
          do jl = 1, nl

            jg = node(jl,ie)
            ju = indx(jg)
            call phi ( jl, xquade, phij, phijx, xleft, xrite )

            aij = he * wquade * &
               ( pp ( xquade ) * phiix * phijx &
               + qq ( xquade ) * phii * phij )
!
!  Decide where the coefficient is to be added.
!
            if ( ju <= 0 ) then

              if ( jg == 0 ) then
                f(iu) = f(iu) - aij * ul
              else if ( jg == n ) then
                f(iu) = f(iu) - aij * ur
              end if

            else if ( iu == ju ) then
              adiag(iu) = adiag(iu) + aij
            else if ( ju < iu ) then
              aleft(iu) = aleft(iu) + aij
            else
              arite(iu) = arite(iu) + aij
            end if

          end do

        end if

      end do

    end do

  end do

  return
end
function ff ( x )

!*****************************************************************************80
!
!! FF evaluates the function F in the differential equation.
!
!  Discussion:
!
!    This is the function F(X) that appears on the right hand
!    side of the equation:
!
!      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) FF, the value of F(X).
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) ff
  integer problem
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
!
!  Find out which problem we're working on.
!
  call get_problem ( problem )

  if ( problem == 1 ) then

    ff = 0.0D+00

  else if ( problem == 2 ) then

    ff = -2.0D+00 * x

  else if ( problem == 3 ) then

    ff = 0.25D+00 * pi**2 * sin ( 0.5D+00 * pi * x )

  else if ( problem == 4 ) then

    ff = 0.25D+00 * pi**2 * cos ( 0.5D+00 * pi * x )

  else if ( problem == 5 ) then

    call get_beta ( beta )

    ff = - ( x**beta ) + ( x**( beta + 2.0D+00 ) ) &
      / ( ( beta + 2.0D+00 ) * ( beta + 1.0D+00 ) )

  else if ( problem == 6 ) then

    call get_alpha ( alpha )
    ff = 2.0D+00 * alpha * ( x - 0.5D+00 ) &
      / ( alpha**2 + ( x - 0.5D+00 )**2 )**2

  end if

  return
end
subroutine geometry ( h, ibc, indx, n, nl, nmax, node, nquad, nu, wquad, xn, &
  xquad )

!*****************************************************************************80
!
!! GEOMETRY sets up some of the geometric information for the problem.
!
!  Discussion:
!
!    Note, however, that the location of the nodes
!    is done outside of this routine, and, in fact, before this
!    routine is called.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) H(N)
!    H(I) is the length of subinterval I.  This code uses
!    equal spacing for all the subintervals.
!
!    Input, integer IBC.
!    IBC declares what the boundary conditions are.
!    1, at the left endpoint, U has the value UL,
!       at the right endpoint, U' has the value UR.
!    2, at the left endpoint, U' has the value UL,
!       at the right endpoint, U has the value UR.
!    3, at the left endpoint, U has the value UL,
!       and at the right endpoint, U has the value UR.
!    4, at the left endpoint, U' has the value UL,
!       at the right endpoint U' has the value UR.
!
!    Output, integer INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    Input, integer N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    Input, integer NMAX, the maximum number of unknowns that can be handled.
!
!    Output, integer NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    Input, integer NQUAD
!    The number of quadrature points used in a subinterval.
!
!    Output, integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Output, real ( kind = 8 ) WQUAD(NQUAD).
!    WQUAD(I) is the weight associated with the I-th point
!    of an NQUAD point Gaussian quadrature rule.
!
!    Input, real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    Output, real ( kind = 8 ) XQUAD(NQUAD,NMAX), the I-th quadrature point
!    in interval J.
!
  implicit none

  integer n
  integer nl
  integer nmax
  integer nquad

  real ( kind = 8 ) alfa
  real ( kind = 8 ) h(n)
  integer i
  integer ibc
  integer igl
  integer igr
  integer indx(0:n)
  integer node(nl,nmax)
  integer nu
  real ( kind = 8 ) wquad(nquad)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xn(0:n)
  real ( kind = 8 ) xquad(nquad,nmax)
  real ( kind = 8 ) xr
!
!  Store in NODE the fact that interval I has node I-1
!  as its left endpoint, and node I as its right endpoint.
!
  do i = 1, n
    node(1,i) = i-1
    node(2,i) = i
  end do
!
!  For every node that is associated with an unknown, we
!  record the number of the unknown in INDX.
!
  nu = 0
  do i = 0, n

    if ( i == 0 .and. ( ibc == 1 .or. ibc == 3 ) ) then
      indx(i) = -1
    else if ( i == n .and. ( ibc == 2 .or. ibc == 3 ) ) then
      indx(i) = -1
    else
      nu = nu+1
      indx(i) = nu
    end if

  end do
!
!  We compute the width of each interval.
!
  do i = 1, n
    igl = node(1,i)
    igr = node(2,i)
    h(i) = xn(igr) - xn(igl)
  end do
!
!  We compute the location of the quadrature points in each
!  interval.
!
  do i = 1, n

    xl = xn(node(1,i))
    xr = xn(node(2,i))

    if ( nquad == 1 ) then
      xquad(1,i) = 0.5D+00 * ( xl + xr )
    else if ( nquad == 2 ) then
      alfa = -0.577350D+00
      xquad(1,i) = ( ( 1.0D+00 - alfa ) * xl   &
                   + ( 1.0D+00 + alfa ) * xr ) &
                   /   2.0D+00
      alfa = +0.577350D+00
      xquad(2,i) = ( ( 1.0D+00 - alfa ) * xl   &
                   + ( 1.0D+00 + alfa ) * xr ) &
                   /   2.0D+00
    else if ( nquad == 3 ) then
      alfa = -0.774597D+00
      xquad(1,i) = ( ( 1.0D+00 - alfa ) * xl   &
                   + ( 1.0D+00 + alfa ) * xr ) &
                   /   2.0D+00
      xquad(2,i) = 0.5D+00 * ( xl + xr )
      alfa = +0.774597D+00
      xquad(3,i) = ( ( 1.0D+00 - alfa ) * xl   &
                   + ( 1.0D+00 + alfa ) * xr ) &
                   /   2.0D+00
    end if

  end do
!
!  Store the weights for the quadrature rule.
!
  if ( nquad == 1 ) then
    wquad(1) = 1.0D+00
  else if ( nquad == 2 ) then
    wquad(1) = 0.5D+00
    wquad(2) = 0.5D+00
  else if ( nquad == 3 ) then
    wquad(1) = 4.0D+00 / 9.0D+00
    wquad(2) = 5.0D+00 / 18.0D+00
    wquad(3) = 4.0D+00 / 9.0D+00
  end if

  return
end
subroutine get_alpha ( alpha )

!*****************************************************************************80
!
!! GET_ALPHA returns the value of ALPHA, for use by problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ALPHA, the value of ALPHA.
!
  implicit none

  real ( kind = 8 ) alpha

  alpha = 0.01D+00

  return
end
subroutine get_beta ( beta )

!*****************************************************************************80
!
!! GET_BETA returns the value of BETA, for use by problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) BETA, the value of BETA.
!
  implicit none

  real ( kind = 8 ) beta

  beta = -0.9D+00

  return
end
subroutine get_problem ( problem )

!*****************************************************************************80
!
!! GETPRB returns the value of the current problem number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, integer PROBLEM, the index of the problem.
!
  implicit none

  integer problem

  problem = 6

  return
end
subroutine init ( ibc, n, tol, ul, ur, xl, xn, xr )

!*****************************************************************************80
!
!! INIT initializes some parameters that define the problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, integer IBC.
!    IBC declares what the boundary conditions are.
!    1, at the left endpoint, U has the value UL,
!       at the right endpoint, U' has the value UR.
!    2, at the left endpoint, U' has the value UL,
!       at the right endpoint, U has the value UR.
!    3, at the left endpoint, U has the value UL,
!       and at the right endpoint, U has the value UR.
!    4, at the left endpoint, U' has the value UL,
!       at the right endpoint U' has the value UR.
!
!    Input, integer N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Output, real ( kind = 8 ) TOL.
!    A tolerance that is used to determine whether the estimated
!    error in an interval is so large that it should be subdivided
!    and the problem solved again.
!
!    Output, real ( kind = 8 ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Output, real ( kind = 8 ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Output, real ( kind = 8 ) XL.
!    XL is the left endpoint of the interval over which the
!    differential equation is being solved.
!
!    Output, real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    Output, real ( kind = 8 ) XR.
!    XR is the right endpoint of the interval over which the
!    differential equation is being solved.
!
  implicit none

  integer n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer i
  integer ibc
  integer problem
  real ( kind = 8 ) tol
  real ( kind = 8 ) u_exact
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) xl
  real ( kind = 8 ) xn(0:n)
  real ( kind = 8 ) xr

  tol = 0.01D+00
!
!  Find out which problem we're working on.
!
  call get_problem ( problem )
!
!  Set the boundary conditions for the problem, and
!  print out its title.
!
  if ( problem == 1 ) then

    ibc = 3
    ul = 0.0D+00
    ur = 1.0D+00
    xl = 0.0D+00
    xr = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Exact solution is U = X'

  else if ( problem == 2 ) then

    ibc = 3
    ul = 0.0D+00
    ur = 1.0D+00
    xl = 0.0D+00
    xr = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Exact solution is U = X*X'

  else if ( problem == 3 ) then

    ibc = 3
    ul = 0.0D+00
    ur = 1.0D+00
    xl = 0.0D+00
    xr = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Exact solution is U = SIN(PI*X/2)'

  else if ( problem == 4 ) then

    ibc = 3
    ul = 1.0D+00
    ur = 0.0D+00
    xl = 0.0D+00
    xr = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Exact solution is U = COS(PI*X/2)'

  else if ( problem == 5 ) then

    ibc = 3
    call get_beta ( beta )
    ul = 0.0D+00
    ur = 1.0D+00 / ( ( beta + 2.0D+00 ) * ( beta + 1.0D+00 ) )
    xl = 0.0D+00
    xr = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Rheinboldt problem'

  else if ( problem == 6 ) then

    ibc = 3
    call get_alpha ( alpha )
    xl = 0.0D+00
    xr = 1.0D+00
    ul = u_exact ( xl )
    ur = u_exact ( xr )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Arctangent problem'

  end if
!
!  The nodes are defined here, and not in the geometry routine.
!  This is because each new iteration chooses the location
!  of the new nodes in a special way.
!
  do i = 0, n
    xn(i) = ( real ( n - i, kind = 8 ) * xl   &
            + real (     i, kind = 8 ) * xr ) &
            / real ( n,     kind = 8 )
  end do

  write ( *, '(a)' ) 'The equation is to be solved for '
  write ( *, '(a,g14.6)' ) 'X greater than ', xl
  write ( *, '(a,g14.6)' ) ' and less than ', xr
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The boundary conditions are:'
  write ( *, '(a)' ) ' '

  if ( ibc == 1 .or. ibc == 3 ) then
    write ( *, '(a,g14.6)' ) '  At X = XL, U=', ul
  else
    write ( *, '(a,g14.6)' ) '  At X = XL, U''=', ul
  end if

  if ( ibc == 2 .or. ibc == 3 ) then
    write ( *, '(a,g14.6)' ) '  At X = XR, U=', ur
  else
    write ( *, '(a,g14.6)' ) '  At X = XR, U''=', ur
  end if

  return
end
subroutine output ( f, ibc, indx, n, nu, ul, ur, xn )

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
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    Input, integer IBC.
!    IBC declares what the boundary conditions are.
!    1, at the left endpoint, U has the value UL,
!       at the right endpoint, U' has the value UR.
!    2, at the left endpoint, U' has the value UL,
!       at the right endpoint, U has the value UR.
!    3, at the left endpoint, U has the value UL,
!       and at the right endpoint, U has the value UR.
!    4, at the left endpoint, U' has the value UL,
!       at the right endpoint U' has the value UR.
!
!    Input, integer INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    Input, integer N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Input, real ( kind = 8 ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Input, real ( kind = 8 ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Input, real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
  implicit none

  integer n
  integer nu

  real ( kind = 8 ) error
  real ( kind = 8 ) f(nu)
  integer i
  integer ibc
  integer indx(0:n)
  real ( kind = 8 ) u
  real ( kind = 8 ) uex
  real ( kind = 8 ) u_exact
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) xn(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Node    X(I)        U(X(I))        U exact       Error'
  write ( *, '(a)' ) ' '

  do i = 0, n

    if ( i == 0 ) then

      if ( ibc == 1 .or. ibc == 3 ) then
        u = ul
      else
        u = f(indx(i))
      end if

    else if ( i == n ) then

      if ( ibc == 2 .or. ibc == 3 ) then
        u = ur
      else
        u = f(indx(i))
      end if

    else

      u = f(indx(i))

    end if

    uex = u_exact ( xn(i) )
    error = u - uex

    write(*,'(i4,4g14.6)') i, xn(i), u, uex, error

  end do

  return
end
subroutine phi ( il, x, phii, phiix, xleft, xrite )

!*****************************************************************************80
!
!! PHI evaluates a linear basis function and its derivative.
!
!  Discussion:
!
!    The functions are evaluated at a point X in an interval.  In any
!    interval, there are just two basis functions.  The first
!    basis function is a line which is 1 at the left endpoint
!    and 0 at the right.  The second basis function is 0 at
!    the left endpoint and 1 at the right.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer IL, the local index of the basis function.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) PHII, PHIIX, the value of the basis function
!    and its derivative.
!
!    Input, real ( kind = 8 ) XLEFT, XRITE, the endpoints of the interval.
!
  implicit none

  integer il
  real ( kind = 8 ) phii
  real ( kind = 8 ) phiix
  real ( kind = 8 ) x
  real ( kind = 8 ) xleft
  real ( kind = 8 ) xrite

  if ( xleft <= x .and. x <= xrite ) then

    if ( il == 1 ) then
      phii = ( xrite - x ) / ( xrite - xleft )
      phiix = -1.0D+00 / ( xrite - xleft )
    else
      phii = ( x - xleft ) / ( xrite - xleft )
      phiix = 1.0D+00 / ( xrite - xleft )
    end if
!
!  If X is outside of the interval, then the basis function
!  is always zero.
!
  else
    phii = 0.0D+00
    phiix = 0.0D+00
  end if

  return
end
function pp ( x )

!*****************************************************************************80
!
!! PP evaluates the function P in the differential equation.
!
!  Discussion:
!
!    The function P(X) occurs in the differential equation:
!
!      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) PP, the value of P(X).
!
  implicit none

  integer problem
  real ( kind = 8 ) pp
  real ( kind = 8 ) x
!
!  Find out which problem we're working on.
!
  call get_problem ( problem )

  if ( problem == 1 ) then
    pp = 1.0D+00
  else if ( problem == 2 ) then
    pp = 1.0D+00
  else if ( problem == 3 ) then
    pp = 1.0D+00
  else if ( problem == 4 ) then
    pp = 1.0D+00
  else if ( problem == 5 ) then
    pp = 1.0D+00
  else if ( problem == 6 ) then
    pp = 1.0D+00
  end if

  return
end
subroutine prsys ( adiag, aleft, arite, f, nu )

!*****************************************************************************80
!
!! PRSYS prints out the tridiagonal linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ADIAG(NU).
!    ADIAG(I) is the "diagonal" coefficient of the I-th
!    equation in the linear system.  That is, ADIAG(I) is
!    the coefficient of the I-th unknown in the I-th equation.
!
!    Input, real ( kind = 8 ) ALEFT(NU).
!    ALEFT(I) is the "left hand" coefficient of the I-th
!    equation in the linear system.  That is, ALEFT(I) is the
!    coefficient of the (I-1)-th unknown in the I-th equation.
!    There is no value in ALEFT(1), since the first equation
!    does not refer to a "0-th" unknown.
!
!    Input, real ( kind = 8 ) ARITE(NU).
!    ARITE(I) is the "right hand" coefficient of the I-th
!    equation in the linear system.  ARITE(I) is the coefficient
!    of the (I+1)-th unknown in the I-th equation.  There is
!    no value in ARITE(NU) because the NU-th equation does not
!    refer to an "NU+1"-th unknown.
!
!    Input, real ( kind = 8 ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
  implicit none

  integer nu

  real ( kind = 8 ) adiag(nu)
  real ( kind = 8 ) aleft(nu)
  real ( kind = 8 ) arite(nu-1)
  real ( kind = 8 ) f(nu)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Printout of tridiagonal linear system:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Equation  A-Left  A-Diag  A-Rite  RHS'
  write ( *, '(a)' ) ' '

  do i = 1, nu
    if ( i == 1 ) then
      write(*,'(i3,14x,3g14.6)')i,adiag(i),arite(i),f(i)
    else if ( i < nu ) then
      write(*,'(i3,4g14.6)')i,aleft(i),adiag(i),arite(i),f(i)
    else
      write(*,'(i3,2g14.6,14x,g14.6)')i,aleft(i),adiag(i),f(i)
    end if
  end do

  return
end
function qq ( x )

!*****************************************************************************80
!
!! QQ evaluates the function Q in the differential equation.
!
!  Discussion:
!
!    The function Q(X) occurs in the differential equation:
!
!      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) QQ, the value of Q(X).
!
  implicit none

  integer problem
  real ( kind = 8 ) qq
  real ( kind = 8 ) x
!
!  Find out which problem we're working on.
!
  call get_problem ( problem )

  if ( problem == 1 ) then
    qq = 0.0D+00
  else if ( problem == 2 ) then
    qq = 0.0D+00
  else if ( problem == 3 ) then
    qq = 0.0D+00
  else if ( problem == 4 ) then
    qq = 0.0D+00
  else if ( problem == 5 ) then
    qq = 1.0D+00
  else if ( problem == 6 ) then
    qq = 0.0D+00
  end if

  return
end
subroutine solve ( adiag, aleft, arite, f, nu )

!*****************************************************************************80
!
!! SOLVE solves a tridiagonal matrix system of the form A*x = b.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) ADIAG(NU), ALEFT(NU), ARITE(NU).
!    On input, ADIAG, ALEFT, and ARITE contain the diagonal,
!    left and right entries of the equations.
!    On output, ADIAG and ARITE have been changed in order
!    to compute the solution.
!    Note that for the first equation, there is no ALEFT
!    coefficient, and for the last, there is no ARITE.
!    So there is no need to store a value in ALEFT(1), nor
!    in ARITE(NU).
!
!    Input/output, real ( kind = 8 ) F(NU).
!    On input, F contains the right hand side of the linear
!    system to be solve
!    On output, F contains the solution of the linear system.
!
!    Input, INTEGER NU, the number of equations to be solved.
!
  implicit none

  integer nu

  real ( kind = 8 ) adiag(nu)
  real ( kind = 8 ) aleft(nu)
  real ( kind = 8 ) arite(*)
  real ( kind = 8 ) f(nu)
  integer i
!
!  Handle the special case of a single equation.
!
  if ( nu == 1 ) then
    f(1) = f(1) / adiag(1)
!
!  The general case, when NU is greater than 1.
!
  else
    arite(1) = arite(1) / adiag(1)
    do i = 2, nu-1
      adiag(i) = adiag(i) - aleft(i) * arite(i-1)
      arite(i) = arite(i) / adiag(i)
    end do
    adiag(nu) = adiag(nu) - aleft(nu) * arite(nu-1)

    f(1) = f(1) / adiag(1)
    do i = 2, nu
      f(i) = ( f(i) - aleft(i) * f(i-1) ) / adiag(i)
    end do

    do i = nu-1, 1, -1
      f(i) = f(i) - arite(i) * f(i+1)
    end do

  end if

  return
end
subroutine solvex ( f, h, ibc, kount, n, nl, nmax, nu, ul, ur, xn )

!*****************************************************************************80
!
!! SOLVEX discretizes and solves a differential equation given the nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    Output, real ( kind = 8 ) H(N)
!    H(I) is the length of subinterval I.  This code uses
!    equal spacing for all the subintervals.
!
!    Input, integer IBC.
!    IBC declares what the boundary conditions are.
!    1, at the left endpoint, U has the value UL,
!       at the right endpoint, U' has the value UR.
!    2, at the left endpoint, U' has the value UL,
!       at the right endpoint, U has the value UR.
!    3, at the left endpoint, U has the value UL,
!       and at the right endpoint, U has the value UR.
!    4, at the left endpoint, U' has the value UL,
!       at the right endpoint U' has the value UR.
!
!    Input, integer KOUNT, the number of adaptive steps that have been taken.
!
!    Input, integer N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    Input, integer NMAX, the maximum number of unknowns that can be handled.
!
!    Output, integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Input, real ( kind = 8 ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Input, real ( kind = 8 ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Input, real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) ADIAG(NU).
!    ADIAG(I) is the "diagonal" coefficient of the I-th
!    equation in the linear system.  That is, ADIAG(I) is
!    the coefficient of the I-th unknown in the I-th equation.
!
!    Local, real ( kind = 8 ) ALEFT(NU).
!    ALEFT(I) is the "left hand" coefficient of the I-th
!    equation in the linear system.  That is, ALEFT(I) is the
!    coefficient of the (I-1)-th unknown in the I-th equation.
!    There is no value in ALEFT(1), since the first equation
!    does not refer to a "0-th" unknown.
!
!    Local, real ( kind = 8 ) ARITE(NU).
!    ARITE(I) is the "right hand" coefficient of the I-th
!    equation in the linear system.  ARITE(I) is the coefficient
!    of the (I+1)-th unknown in the I-th equation.  There is
!    no value in ARITE(NU) because the NU-th equation does not
!    refer to an "NU+1"-th unknown.
!
!    Local, integer INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    Local, integer NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    Local, integer NQUAD
!    The number of quadrature points used in a subinterval.
!
!    Local, real ( kind = 8 ) WQUAD(NQUAD).
!    WQUAD(I) is the weight associated with the I-th point
!    of an NQUAD point Gaussian quadrature rule.
!
!    Local, real ( kind = 8 ) XQUAD(NQUAD,NMAX), the I-th quadrature point
!    in interval J.
!
  implicit none

  integer n
  integer nl
  integer nmax
  integer, parameter :: nquad = 2

  real ( kind = 8 ) adiag(nmax)
  real ( kind = 8 ) aleft(nmax)
  real ( kind = 8 ) arite(nmax)
  real ( kind = 8 ) f(nmax)
  real ( kind = 8 ) h(n)
  integer ibc
  integer indx(0:n)
  integer kount
  integer node(nl,nmax)
  integer nu
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) wquad(nquad)
  real ( kind = 8 ) xn(0:n)
  real ( kind = 8 ) xquad(nquad,nmax)
!
!  Given a set of N nodes (where N increases on each iteration),
!  compute the other geometric information.
!
  call geometry ( h, ibc, indx, n, nl, nmax, node, nquad, nu, wquad, xn, xquad )
!
!  Assemble the linear system.
!
  call assemble ( adiag, aleft, arite, f, h, n, indx, node, nu, nl, &
    nquad, nmax, ul, ur, wquad, xn, xquad )
!
!  Print out the linear system, just once.
!
  if ( kount == 1 ) then
    call prsys ( adiag, aleft, arite, f, nu )
  end if
!
!  Solve the linear system.
!
  call solve ( adiag, aleft, arite, f, nu )
!
!  Print out the solution.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Basic solution'

  call output ( f, ibc, indx, n, nu, ul, ur, xn )

  return
end
subroutine solvey ( eta, f, h, n, nu, ul, ur, xn )

!*****************************************************************************80
!
!! SOLVEY computes error estimators for a finite element solution.
!
!  Discussion:
!
!    SOLVEY accepts information about the solution of a finite element
!    problem on a grid of nodes with coordinates XN.  It then starts
!    at node 0, and for each node, computes two "error estimators",
!    one for the left, and one for the right interval associated with the
!    node.  These estimators are found by solving a finite element problem
!    over the two intervals, using the known values of the original
!    solution as boundary data, and using a mesh that is "slightly"
!    refined over the original one.
!
!    Note that the computations at the 0-th and N-th nodes only involve
!    a single interval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ETA(N).
!    ETA(I) is the error estimate for interval I.  It is computed
!    as the sum of two quantities, one associated with the left
!    and one with the right node of the interval.
!
!    Input, real ( kind = 8 ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    Input, real ( kind = 8 ) H(N)
!    H(I) is the length of subinterval I.  This code uses
!    equal spacing for all the subintervals.
!
!    Input, integer N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Input, real ( kind = 8 ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Input, real ( kind = 8 ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Input, real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
  implicit none

  integer, parameter :: nl = 2
  integer, parameter :: ny = 2
  integer, parameter :: nquad = 2

  integer, parameter :: nmay = 2 * ny

  integer n
  integer nu

  real ( kind = 8 ) adiag(nmay)
  real ( kind = 8 ) aleft(nmay)
  real ( kind = 8 ) arite(nmay)
  real ( kind = 8 ) eta(n)
  real ( kind = 8 ) f(nu)
  real ( kind = 8 ) fy(nmay)
  real ( kind = 8 ) h(n)
  real ( kind = 8 ) hy(nmay)
  integer i
  integer ibcy
  integer indy(0:nmay)
  integer j
  integer jhi
  integer jlo
  integer jmid
  integer k
  integer m
  integer nodey(nl,nmay)
  integer nuy
  real ( kind = 8 ) pp
  real ( kind = 8 ) qq
  real ( kind = 8 ) total
  real ( kind = 8 ) ul
  real ( kind = 8 ) uleft
  real ( kind = 8 ) ulval
  real ( kind = 8 ) uly
  real ( kind = 8 ) uprime
  real ( kind = 8 ) ur
  real ( kind = 8 ) urite
  real ( kind = 8 ) urval
  real ( kind = 8 ) ury
  real ( kind = 8 ) uval
  real ( kind = 8 ) vlval
  real ( kind = 8 ) vprime
  real ( kind = 8 ) vrval
  real ( kind = 8 ) vval
  real ( kind = 8 ) wquad(nquad)
  real ( kind = 8 ) xn(0:n)
  real ( kind = 8 ) xquady(nquad,nmay)
  real ( kind = 8 ) y
  real ( kind = 8 ) yl
  real ( kind = 8 ) ym
  real ( kind = 8 ) yn(0:nmay)
  real ( kind = 8 ) yr
!
!  Initialize the error estimators to zero.
!
  eta(1:n) = 0.0D+00
!
!  Set the boundary conditions for each subproblem to be
!  known values of U at the left and right.
!
!
!  For each node, subdivide its left and right hand intervals
!  into NY subintervals.
!
!  Set up and solve the differential equation again on this
!  smaller region.
!
!  The 0-th and N-th nodes are special cases.
!
  ibcy = 3

  do j = 0, n

    if ( j == 0 ) then
      m = ny
      jlo = j
      jmid = j + 1
      jhi = j + 1
    else if ( j == n ) then
      m = ny
      jlo = j - 1
      jmid = j
      jhi = j
    else
      m = 2 * ny
      jlo = j - 1
      jmid = j
      jhi = j + 1
    end if
!
!  Set the location of the nodes in the subintervals.
!
    yl = xn(jlo)
    ym = xn(jmid)
    yr = xn(jhi)

    do i = 0, ny
      yn(i) = ( real ( ny - i, kind = 8 ) * yl   &
              + real (      i, kind = 8 ) * ym ) &
              / real ( ny,     kind = 8 )
    end do

    do i = ny+1, m
      yn(i) = ( real ( m - i,      kind = 8 ) * ym   &
              + real (     i - ny, kind = 8 ) * yr ) &
              / real ( m -     ny, kind = 8 )
    end do
!
!  Set up the geometry of the sub-problem.
!
    call geometry ( hy, ibcy, indy, m, nl, nmay, nodey, nquad, nuy, &
      wquad, yn, xquady )
!
!  Set the boundary values for the sub-problem.
!
    if ( j <= 1 ) then
      uly = ul
    else
      uly = f(j-1)
    end if

    if ( n - 1 <= j ) then
      ury = ur
    else
      ury = f(j+1)
    end if
!
!  Assemble the matrix for the sub-problem.
!
    call assemble ( adiag, aleft, arite, fy, hy, m, indy, nodey, nuy, nl, &
      nquad, nmay, uly, ury, wquad, yn, xquady )
!
!  Solve the system.
!
    call solve ( adiag, aleft, arite, fy, nuy )
!
!  Compute the weighted sum of the squares of the differences
!  of the original computed slope and the refined computed slopes.
!
!  Calculation for left interval.
!
    if ( 1 <= j ) then

      if ( j <= 1 ) then
        uleft = ul
        urite = f(1)
      else if ( j == n ) then
        uleft = f(j-1)
        urite = ur
      else
        uleft = f(j-1)
        urite = f(j)
      end if

      uprime = ( urite - uleft ) / h(j)

      total = 0.0D+00
      do i = 1, ny

        yl = yn(i-1)
        yr = yn(i)

        if ( i == 1 ) then
          vlval = uly
          vrval = fy(i)
        else if ( i == m ) then
          vlval = fy(i-1)
          vrval = ury
        else
          vlval = fy(i-1)
          vrval = fy(i)
        end if

        vprime = ( vrval - vlval ) / hy(i)

        ulval = ( real ( ny - i + 1, kind = 8 ) * uleft   &
                + real (      i - 1, kind = 8 ) * urite ) &
                / real ( ny,         kind = 8 )

        urval = ( real ( ny - i, kind = 8 ) * uleft   &
                + real (      i, kind = 8 ) * urite ) &
                / real ( ny,     kind = 8 )
!
!  Compute the integral of
!
!    p(x)*(u'(x)-v'(x))**2 + q(x)*(u(x)-v(x))**2
!
        do k = 1, nquad

          y  =  xquady(k,i)

          uval = ( ( yl - y      ) * urval   &
                 + (      y - yr ) * ulval ) &
                 / ( yl     - yr )

          vval = ( ( yl - y      ) * vrval   &
                 + (      y - yr ) * vlval ) &
                 / ( yl     - yr )

          total = total + 0.5D+00 * wquad(k) * hy(i) * &
            ( pp ( y ) * ( uprime - vprime )**2 &
            + qq ( y ) * ( uval - vval )**2 )

        end do

      end do

      eta(j) = eta(j) + 0.5D+00 * sqrt ( total )

    end if
!
!  Calculation for right interval.
!
    if ( j <= n - 1 ) then

      if ( j == 0 ) then
        uleft = ul
        urite = f(j+1)
      else if ( n - 1 <= j ) then
        uleft = f(j)
        urite = ur
      else
        uleft = f(j)
        urite = f(j+1)
      end if

      uprime = ( urite - uleft ) / h(j+1)

      total = 0.0D+00
      do i = m+1-ny, m

        yl = yn(i-1)
        yr = yn(i)

        if ( i == 1 ) then
          vlval = uly
          vrval = fy(i)
        else if ( i == m ) then
          vlval = fy(i-1)
          vrval = ury
        else
          vlval = fy(i-1)
          vrval = fy(i)
        end if

        vprime = ( vrval - vlval ) / hy(i)

        ulval = ( real (      m - i + 1, kind = 8 ) * uleft   &
                + real ( ny - m + i - 1, kind = 8 ) * urite ) &
                / real ( ny,             kind = 8 )

        urval = ( real (      m - i, kind = 8 ) * uleft   &
                + real ( ny - m + i, kind = 8 ) * urite ) &
                / real ( ny,         kind = 8 )
!
!  Compute the integral of
!
!    p(x)*(u'(x)-v'(x))**2 + q(x)*(u(x)-v(x))**2
!
        do k = 1, nquad

          y  =  xquady(k,i)

          uval = ( ( yl - y      ) * urval   &
                 + (      y - yr ) * ulval ) &
                 / ( yl     - yr )

          vval = ( ( yl - y      ) * vrval   &
                 + (      y - yr ) * vlval ) &
                 / ( yl     - yr )

          total = total + 0.5D+00 * wquad(k) * hy(i) * &
            ( pp ( y ) * ( uprime - vprime )**2 &
            + qq ( y ) * ( uval - vval )**2 )

        end do

      end do

      eta(j+1) = eta(j+1) + 0.5D+00 * sqrt ( total )

    end if

  end do
!
!  Print out the error estimators.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ETA'
  write ( *, '(a)' ) ' '
  do j = 1, n
    write ( *, '(g14.6)' ) eta(j)
  end do

  return
end
subroutine subdiv ( eta, kount, n, nmax, tol, xn, status )

!*****************************************************************************80*
!
!! SUBDIV decides which intervals should be subdivided.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ETA(N).
!    ETA(I) is the error estimate for interval I.  It is computed
!    as the sum of two quantities, one associated with the left
!    and one with the right node of the interval.
!
!    Input, integer KOUNT, the number of adaptive steps that have been taken.
!
!    Input/output, integer N
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer NMAX, the maximum number of unknowns that can be handled.
!
!    Input, real ( kind = 8 ) TOL.
!    A tolerance that is used to determine whether the estimated
!    error in an interval is so large that it should be subdivided
!    and the problem solved again.
!
!    Input/output, real ( kind = 8 ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    Output, integer STATUS, reports status of subdivision.
!    0, a new subdivision was carried out.
!    1, no more subdivisions are needed.
!    -1, no more subdivisions can be carried out.
!
!  Local Parameters:
!
!    Local, integer JADD(N).
!    JADD(I) is 1 if the error estimates show that interval I
!    should be subdivided.
!
  implicit none

  integer n
  integer nmax

  real ( kind = 8 ) ave
  real ( kind = 8 ) eta(n)
  integer j
  integer jadd(n)
  integer k
  integer kount
  integer status
  real ( kind = 8 ) temp
  real ( kind = 8 ) tol
  real ( kind = 8 ) xn(0:nmax)
  real ( kind = 8 ) xt(0:nmax)

  status = 0
!
!  Add up the ETA's, and get their average.
!
  ave = sum ( eta(1:n) ) / real ( n, kind = 8 )
!
!  Look for intervals whose ETA value is relatively large,
!  and note in JADD that these intervals should be subdivided.
!
  k = 0
  temp = max ( 1.2D+00 * ave + 0.00001D+00, tol**2 / real ( n, kind = 8 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) 'Tolerance = ', temp
  write ( *, '(a)' ) ' '

  do j = 1, n

    if ( temp < eta(j) ) then
      k = k + 1
      jadd(j) = 1
      write ( *, '(a,i8)' ) 'Subdivide interval ', j
    else
      jadd(j) = 0
    end if

  end do
!
!  If no subdivisions needed, we're done.
!
  if ( k <= 0 ) then
    write ( *, '(a,i8)' ) 'Success on step ', kount
    status = 1
    return
  end if
!
!  See if we're about to go over our limit.
!
  if ( nmax < n + k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'The iterations did not reach their goal.'
    write ( *, '(a,i8)' ) 'The next value of N is ', n + k
    write ( *, '(a,i8)' ) 'which exceeds NMAX = ', nmax
    status = -1
    return
  end if
!
!  Insert new nodes where needed.
!
  k = 0
  xt(0) = xn(0)
  do j = 1, n

    if ( 0 < jadd(j) ) then
      xt(j+k) = 0.5D+00 * ( xn(j) + xn(j-1) )
      k = k + 1
    end if

    xt(j+k) = xn(j)

  end do
!
!  Update the value of N, and copy the new nodes into XN.
!
  n = n + k

  xn(0:n) = xt(0:n)

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
function u_exact ( x )

!*****************************************************************************80
!
!! U_EXACT returns the value of the exact solution at any point X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) U_EXACT, the value of the exact solution at X.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer problem
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u_exact
  real ( kind = 8 ) x
!
!  Find out which problem we're working on.
!
  call get_problem ( problem )

  if ( problem == 1 ) then
    u_exact = x
  else if ( problem == 2 ) then
    u_exact = x**2
  else if ( problem == 3 ) then
    u_exact = sin ( pi * x / 2.0D+00 )
  else if ( problem == 4 ) then
    u_exact = cos ( pi * x / 2.0D+00 )
  else if ( problem == 5 ) then
    call get_beta ( beta )
    u_exact = ( x**( beta + 2.0D+00 ) ) &
      / ( ( beta + 2.0D+00 ) * ( beta + 1.0D+00 ) )
  else if ( problem == 6 ) then
    call get_alpha ( alpha )
    u_exact = atan ( ( x - 0.5D+00 ) / alpha )
  else
    u_exact = 0.0D+00
  end if

  return
end
