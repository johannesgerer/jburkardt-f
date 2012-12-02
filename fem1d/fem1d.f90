program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM1D.
!
!  Discussion:
!
!    FEM1D solves a one dimensional ODE using the finite element method.
!
!    The differential equation solved is
!
!      - d/dX (P dU/dX) + Q U  =  F
!
!    The finite-element method uses piecewise linear basis functions.
!
!    Here U is an unknown scalar function of X defined on the
!    interval [XL,XR], and P, Q and F are given functions of X.
!
!    The values of U or U' at XL and XR are also specified.
!
!
!    The interval [XL,XR] is "meshed" with NSUB+1 points,
!
!    XN(0) = XL, XN(1)=XL+H, XN(2)=XL+2*H, ..., XN(NSUB)=XR.
!
!    This creates NSUB subintervals, with interval number 1
!    having endpoints XN(0) and XN(1), and so on up to interval
!    NSUB, which has endpoints XN(NSUB-1) and XN(NSUB).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) ADIAG(NU), the "diagonal" coefficients. 
!    That is, ADIAG(I) is the coefficient of the I-th unknown in 
!    the I-th equation.
!
!    real ( kind = 8 ) ALEFT(NU), the "left hand" coefficients.  That is, 
!    ALEFT(I) is the coefficient of the (I-1)-th unknown in the I-th equation.
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
!    real ( kind = 8 ) F(NSUB+1) or F(NU).
!    ASSEMB stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    real ( kind = 8 ) H(NSUB)
!    H(I) is the length of subinterval I.  This code uses
!    equal spacing for all the subintervals.
!
!    integer ( kind = 4 ) IBC.
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
!    integer ( kind = 4 ) INDX(0:NSUB).
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
!    at node NSUB, which will be unknown NSUB or NSUB+1,
!    depending on whether there was an unknown at node 0.
!
!    integer ( kind = 4 ) NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    integer ( kind = 4 ) NODE(NL,NSUB).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    integer ( kind = 4 ) NQUAD.
!    The number of quadrature points used in a subinterval.
!    This code uses NQUAD = 1.
!
!    integer ( kind = 4 ) NSUB.
!    The number of subintervals into which the interval [XL,XR] is broken.
!
!    integer ( kind = 4 ) NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be NSUB-1,
!    NSUB, or NSUB+1 unknown values, which are the coefficients
!    of basis functions.
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
!    real ( kind = 8 ) XL.
!    XL is the left endpoint of the interval over which the
!    differential equation is being solved.
!
!    real ( kind = 8 ) XN(0:NSUB).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(NSUB) is XR.
!
!    real ( kind = 8 ) XQUAD(NSUB)
!    XQUAD(I) is the location of the single quadrature point
!    in interval I.
!
!    real ( kind = 8 ) XR.
!    XR is the right endpoint of the interval over which the
!    differential equation is being solved.
!
  implicit none
!
!  Set the number of subintervals.
!
  integer ( kind = 4 ), parameter :: nsub = 5
!
!  Set the number of basis functions per interval.
!  2 means linear functions are used.
!
  integer ( kind = 4 ), parameter :: nl = 2

  real ( kind = 8 ) adiag(nsub+1)
  real ( kind = 8 ) aleft(nsub+1)
  real ( kind = 8 ) arite(nsub+1)
  real ( kind = 8 ) f(nsub+1)
  real ( kind = 8 ) h(nsub)
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) indx(0:nsub)
  integer ( kind = 4 ) node(nl,nsub)
  integer ( kind = 4 ) nquad
  integer ( kind = 4 ) nu
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) xl
  real ( kind = 8 ) xn(0:nsub)
  real ( kind = 8 ) xquad(nsub)
  real ( kind = 8 ) xr

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solve the two-point boundary value problem'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  - d/dX (P dU/dX) + Q U  =  F'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  on the interval [XL,XR], specifying'
  write ( *, '(a)' ) '  the value of U or U'' at each end.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  The interval [XL,XR] is broken into NSUB = ', &
    nsub, ' subintervals'
  write ( *, '(a,i8)' ) '  Number of basis functions per element is NL = ', nl
!
!  Initialize the data that defines the problem.
!
  call init ( ibc, nquad, ul, ur, xl, xr )
!
!  Compute the quantities which define the geometry of the
!  problem.
!
  call geometry ( h, ibc, indx, nl, node, nsub, nu, xl, xn, xquad, xr )
!
!  Assemble the linear system.
!
  call assemble ( adiag, aleft, arite, f, h, indx, nl, node, nu, nquad, &
    nsub, ul, ur, xn, xquad )
!
!  Print out the linear system.
!
  call prsys ( adiag, aleft, arite, f, nu )
!
!  Solve the linear system.
!
  call solve ( adiag, aleft, arite, f, nu )
!
!  Print out the solution.
!
  call output ( f, ibc, indx, nsub, nu, ul, ur, xn )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine assemble ( adiag, aleft, arite, f, h, indx, nl, node, nu, &
  nquad, nsub, ul, ur, xn, xquad )

!*****************************************************************************80
!
!! ASSEMBLE assembles the matrix and right-hand-side of the linear system.
!
!  Discussion:
!
!    The linear system has the form:
!
!      K * C  =  F
!
!    that is to be solved for the coefficients C.
!
!    Numerical integration is used to compute the entries of K and F.
!
!    Note that a 1 point quadrature rule, which is sometimes used to
!    assemble the matrix and right hand side, is just barely accurate
!    enough for simple problems.  If you want better results, you
!    should use a quadrature rule that is more accurate.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ADIAG(NU), the "diagonal" coefficients. 
!    That is, ADIAG(I) is the coefficient of the I-th unknown in 
!    the I-th equation.
!
!    Output, real ( kind = 8 ) ALEFT(NU), the "left hand" coefficients.  That is, 
!    ALEFT(I) is the coefficient of the (I-1)-th unknown in the I-th equation.
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
!    Output, real ( kind = 8 ) F(NSUB+1) or F(NU).
!    ASSEMB stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    Input, real ( kind = 8 ) H(NSUB)
!    H(I) is the length of subinterval I.  This code uses
!    equal spacing for all the subintervals.
!
!    Input, integer ( kind = 4 ) INDX(0:NSUB).
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
!    at node NSUB, which will be unknown NSUB or NSUB+1,
!    depending on whether there was an unknown at node 0.
!
!    Input, integer ( kind = 4 ) NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    Input, integer ( kind = 4 ) NODE(NL,NSUB).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    Input, integer ( kind = 4 ) NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be NSUB-1,
!    NSUB, or NSUB+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Input, integer ( kind = 4 ) NQUAD.
!    The number of quadrature points used in a subinterval.
!    This code uses NQUAD = 1.
!
!    Input, integer ( kind = 4 ) NSUB.
!    The number of subintervals into which the interval [XL,XR] is broken.
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
!    Input, real ( kind = 8 ) XN(0:NSUB).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(NSUB) is XR.
!
!    Input, real ( kind = 8 ) XQUAD(NSUB)
!    XQUAD(I) is the location of the single quadrature point
!    in interval I.
!
  implicit none

  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) nu

  real ( kind = 8 ) aij
  real ( kind = 8 ) adiag(nu)
  real ( kind = 8 ) aleft(nu)
  real ( kind = 8 ) arite(nu)
  real ( kind = 8 ) f(nu)
  real ( kind = 8 ) ff
  real ( kind = 8 ) h(nsub)
  real ( kind = 8 ) he
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) ig
  integer ( kind = 4 ) il
  integer ( kind = 4 ) indx(0:nsub)
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iu
  integer ( kind = 4 ) jg
  integer ( kind = 4 ) jl
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) node(nl,nsub)
  integer ( kind = 4 ) nquad
  real ( kind = 8 ) phii
  real ( kind = 8 ) phiix
  real ( kind = 8 ) phij
  real ( kind = 8 ) phijx
  real ( kind = 8 ) pp
  real ( kind = 8 ) qq
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) x
  real ( kind = 8 ) xleft
  real ( kind = 8 ) xn(0:nsub)
  real ( kind = 8 ) xquad(nsub)
  real ( kind = 8 ) xquade
  real ( kind = 8 ) xrite
!
!  Zero out the arrays that hold the coefficients of the matrix
!  and the right hand side.
!
  f(1:nu) = 0.0D+00
  adiag(1:nu) = 0.0D+00
  aleft(1:nu) = 0.0D+00
  arite(1:nu) = 0.0D+00
!
!  For interval number IE,
!
  do ie = 1, nsub

    he = h(ie)
    xleft = xn(node(1,ie))
    xrite = xn(node(2,ie))
!
!  consider each quadrature point IQ,
!
    do iq = 1, nquad

      xquade = xquad(ie)
!
!  and evaluate the integrals associated with the basis functions
!  for the left, and for the right nodes.
!
      do il = 1, nl

        ig = node(il,ie)
        iu = indx(ig)

        if ( 0 < iu ) then

          call phi ( il, xquade, phii, phiix, xleft, xrite )

          f(iu) = f(iu) + he * ff ( xquade ) * phii
!
!  Take care of boundary nodes at which U' was specified.
!
          if ( ig == 0 ) then

            x = 0.0D+00
            f(iu) = f(iu) - pp ( x ) * ul

          else if ( ig == nsub ) then

            x = 1.0D+00
            f(iu) = f(iu) + pp ( x ) * ur

          end if
!
!  Evaluate the integrals that take a product of the basis
!  function times itself, or times the other basis function
!  that is nonzero in this interval.
!
          do jl = 1, nl

            jg = node(jl,ie)
            ju = indx(jg)

            call phi ( jl, xquade, phij, phijx, xleft, xrite )

            aij = he * ( pp ( xquade ) * phiix * phijx &
                       + qq ( xquade ) * phii  * phij   )
!
!  If there is no variable associated with the node, then it's
!  a specified boundary value, so we multiply the coefficient
!  times the specified boundary value and subtract it from the
!  right hand side.
!
            if ( ju <= 0 ) then

              if ( jg == 0 ) then

                f(iu) = f(iu) - aij * ul

              else if ( jg == nsub ) then

                f(iu) = f(iu) - aij * ur

              end if
!
!  Otherwise, we add the coefficient we've just computed to the
!  diagonal, or left or right entries of row IU of the matrix.
!
            else

              if ( iu == ju ) then
                adiag(iu) = adiag(iu) + aij
              else if ( ju < iu ) then
                aleft(iu) = aleft(iu) + aij
              else
                arite(iu) = arite(iu) + aij
              end if

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
!! FF evaluates the right hand side function.
!
!  Discussion:
!
!    This routine evaluates the function F(X) in the differential equation.
!
!      -d/dx (p du/dx) + q u  =  f
!
!    at the point X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FF, the value of the function.
!
  implicit none

  real ( kind = 8 ) ff
  real ( kind = 8 ) x

  ff = 0.0D+00

  return
end
subroutine geometry ( h, ibc, indx, nl, node, nsub, nu, xl, xn, xquad, xr )

!*****************************************************************************80
!
!! GEOMETRY sets up the geometry for the interval [XL,XR].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) H(NSUB)
!    H(I) is the length of subinterval I.  This code uses
!    equal spacing for all the subintervals.
!
!    Input, integer ( kind = 4 ) IBC.
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
  implicit none

  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nsub

  real ( kind = 8 ) h(nsub)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) indx(0:nsub)
  integer ( kind = 4 ) node(nl,nsub)
  integer ( kind = 4 ) nu
  real ( kind = 8 ) xl
  real ( kind = 8 ) xn(0:nsub)
  real ( kind = 8 ) xquad(nsub)
  real ( kind = 8 ) xr
!
!  Set the value of XN, the locations of the nodes.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node      Location'
  write ( *, '(a)' ) ' '
  do i = 0, nsub
    xn(i)  =  ( real ( nsub - i, kind = 8 ) * xl   &
              + real (        i, kind = 8 ) * xr ) &
              / real ( nsub,     kind = 8 )
    write ( *, '(2x,i8,2x,g14.6)' ), i, xn(i)
  end do
!
!  Set the lengths of each subinterval.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Subint    Length'
  write ( *, '(a)' ) ' '
  do i = 1, nsub
    h(i) = xn(i) - xn(i-1)
    write ( *, '(2x,i8,2x,g14.6)' ) i, h(i)
  end do
!
!  Set the quadrature points, each of which is the midpoint
!  of its subinterval.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Subint    Quadrature point'
  write ( *, '(a)' ) ' '
  do i = 1, nsub
    xquad(i) = 0.5D+00 * ( xn(i-1) + xn(i) )
    write ( *, '(2x,i8,2x,g14.6)' ) i, xquad(i)
  end do
!
!  Set the value of NODE, which records, for each interval,
!  the node numbers at the left and right.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Subint  Left Node  Right Node'
  write ( *, '(a)' ) ' '
  do i = 1, nsub
    node(1,i) = i-1
    node(2,i) = i
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, node(1,i), node(2,i)
  end do
!
!  Starting with node 0, see if an unknown is associated with
!  the node.  If so, give it an index.
!
  nu = 0
!
!  Handle first node.
!
  i = 0
  if ( ibc == 1 .or. ibc == 3 ) then
    indx(i) = -1
  else
    nu = nu+1
    indx(i) = nu
  end if
!
!  Handle nodes 1 through nsub-1
!
  do i = 1, nsub-1
    nu = nu + 1
    indx(i) = nu
  end do
!
!  Handle the last node.
!
  i = nsub

  if ( ibc == 2 .or. ibc == 3 ) then
    indx(i) = -1
  else
    nu = nu+1
    indx(i) = nu
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of unknowns NU = ', nu
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node  Unknown'
  write ( *, '(a)' ) ' '
  do i = 0, nsub
    write ( *, '(2x,i8,2x,i8)' )  i, indx(i)
  end do

  return
end
subroutine init ( ibc, nquad, ul, ur, xl, xr )

!*****************************************************************************80
!
!! INIT assigns values to variables which define the problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IBC.
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
!    Output, integer ( kind = 4 ) NQUAD.
!    The number of quadrature points used in a subinterval.
!    This code uses NQUAD = 1.
!
  implicit none

  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) nquad
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr
!
!  IBC declares what the boundary conditions are.
!
  ibc = 1
!
!  NQUAD is the number of quadrature points per subinterval.
!  The program as currently written cannot handle any value for
!  NQUAD except 1!
!
  nquad = 1
!
!  Set the values of U or U' at the endpoints.
!
  ul = 0.0D+00
  ur = 1.0D+00
!
!  Define the location of the endpoints of the interval.
!
  xl = 0.0D+00
  xr = 1.0D+00
!
!  Print out the values that have been set.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The equation is to be solved for'
  write ( *, '(a,g14.6)' ) 'X greater than XL = ', xl
  write ( *, '(a,g14.6)' ) ' and less than XR = ', xr
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
    write ( *, '(a,g14.6)' ), '  At X = XR, U''=', ur
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'Number of quadrature points per element is ', nquad

  return
end
subroutine output ( f, ibc, indx, nsub, nu, ul, ur, xn )

!*****************************************************************************80
!
!! OUTPUT prints out the computed solution.
!
!  Discussion:
!
!    We simply print out the solution vector F, except that, for
!    certain boundary conditions, we are going to have to get the
!    value of the solution at XL or XR by using the specified
!    boundary value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IBC.
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
  implicit none

  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) nu

  real ( kind = 8 ) f(nu)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) indx(0:nsub)
  real ( kind = 8 ) u
  real ( kind = 8 ) ul
  real ( kind = 8 ) ur
  real ( kind = 8 ) xn(0:nsub)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Computed solution coefficients:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Node    X(I)        U(X(I))'
  write ( *, '(a)' ) ' '

  do i = 0, nsub
!
!  If we're at the first node, check the boundary condition.
!
    if ( i == 0 ) then

      if ( ibc == 1 .or. ibc == 3 ) then
        u = ul
      else
        u = f(indx(i))
      end if
!
!  If we're at the last node, check the boundary condition.
!
    else if ( i == nsub ) then

      if ( ibc == 2 .or. ibc == 3 ) then
        u = ur
      else
        u = f(indx(i))
      end if
!
!  Any other node, we're sure the value is stored in F.
!
    else

      u = f(indx(i))

    end if

    write ( *, '(i8,f6.2,g14.6)' ) i, xn(i), u

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
!    The evaluation is done at a point X in an interval [XLEFT,XRITE].
!
!    In this interval, there are just two nonzero basis functions.
!    The first basis function is a line which is 1 at the left
!    endpoint and 0 at the right.  The second basis function is 0 at
!    the left endpoint and 1 at the right.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) il
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
!  If X is outside of the interval, just set everything to 0.
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
!    The function P appears in the differential equation as;
!
!      - d/dx (p du/dx) + q u  =  f
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) PP, the value of the function.
!
  implicit none

  real ( kind = 8 ) pp
  real ( kind = 8 ) x

  pp = 1.0D+00

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
!    01 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) nu

  real ( kind = 8 ) adiag(nu)
  real ( kind = 8 ) aleft(nu)
  real ( kind = 8 ) arite(nu)
  real ( kind = 8 ) f(nu)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Printout of tridiagonal linear system:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Equation  ALEFT  ADIAG  ARITE  RHS'
  write ( *, '(a)' ) ' '

  do i = 1, nu
    write ( *, '(2x,i8,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, aleft(i), adiag(i), arite(i), f(i)
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
!    The function Q appears in the differential equation as:
!
!      - d/dx (p du/dx) + q u  =  f
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) QQ, the value of the function.
!
  implicit none

  real ( kind = 8 ) qq
  real ( kind = 8 ) x

  qq = 0.0D+00

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
!    01 November 2006
!
!  Author:
!
!    John Burkardt
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
!    system to be solved.
!    On output, F contains the solution of the linear system.
!
!    Input, integer ( kind = 4 ) NU, the number of equations to be solved.
!
  implicit none

  integer ( kind = 4 ) nu

  real ( kind = 8 ) adiag(nu)
  real ( kind = 8 ) aleft(nu)
  real ( kind = 8 ) arite(nu-1)
  real ( kind = 8 ) f(nu)
  integer ( kind = 4 ) i
!
!  Carry out Gauss elimination on the matrix, saving information
!  needed for the backsolve.
!
  arite(1) = arite(1) / adiag(1)
  do i = 2, nu - 1
    adiag(i) = adiag(i) - aleft(i) * arite(i-1)
    arite(i) = arite(i) / adiag(i)
  end do
  adiag(nu) = adiag(nu) - aleft(nu) * arite(nu-1)
!
!  Carry out the same elimination steps on F that were done to the
!  matrix.
!
  f(1) = f(1) / adiag(1)
  do i = 2, nu
    f(i) = ( f(i) - aleft(i) * f(i-1) ) / adiag(i)
  end do
!
!  And now carry out the steps of "back substitution".
!
  do i = nu-1, 1, -1
    f(i) = f(i) - arite(i) * f(i+1)
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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  character ( len = 8 )  date
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  character ( len = 10 ) time
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y
  character ( len = 5 )  zone

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
